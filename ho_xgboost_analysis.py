"""Signed-weight HO diphoton training and independent cut comparison.

The CLI validates complete ROOT populations before calling this module. ML
packages are imported only when this optional workflow is requested.
"""
from __future__ import annotations

import csv
import hashlib
import html
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import tempfile

import numpy as np

from ho_xgboost_core import (
    class_balanced_abs_weights, grouped_sample_split, project_partition_weights,
    scan_thresholds, signed_yield_summary,
)

FEATURE_NAMES = ["m_gg", "pt_gamma1", "eta_gamma1", "pt_gamma2", "eta_gamma2",
                 "deltaR_gg", "deltaPhi_gg", "pt_gg", "y_gg"]
DEFAULT_MODEL_PARAMS = {
    "objective": "binary:logistic", "eval_metric": "logloss", "n_estimators": 300,
    "max_depth": 3, "learning_rate": 0.05, "subsample": 0.9,
    "colsample_bytree": 0.9, "n_jobs": 1,
}
SUMMARY_FIELDS = [
    "sample", "category", "analysis", "detector_response", "response_mode", "input_file",
    "raw_cross_section_pb", "cross_section_pb", "weight_scale", "entries_read",
    "selected_entries", "sum_weight", "sum_selected_weight", "analysis_efficiency",
    "selected_cross_section_pb", "expected_events", "mc_events_after_analysis",
    "normalization_kind", "normalization_sum_weight", "ihixs_record_sha256", "shower_quality",
    "evaluation_partition", "partition_fraction", "partition_entries", "selected_source_events",
    "mc_variance_events", "mc_error_events", "effective_count", "positive_events", "negative_events",
]


def validate_config(config):
    validation = float(config.get("validation_size", 0.2))
    test = float(config.get("test_size", 0.2))
    if not (math.isfinite(validation) and math.isfinite(test)
            and validation > 0 and test > 0 and validation + test < 1):
        raise ValueError("HO validation_size and test_size must be positive and sum to less than one")
    for key, default in (("systematics", 0.0), ("min_background_effective_count", 25.0)):
        value = float(config.get(key, default))
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"xgboost.{key} must be finite and nonnegative")
        if key == "min_background_effective_count" and value == 0:
            raise ValueError("xgboost.min_background_effective_count must be positive")
    seed = config.get("seed", 12345)
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("xgboost.seed must be a nonnegative integer")
    params = config.get("model_params", {}) or {}
    if not isinstance(params, dict):
        raise ValueError("xgboost.model_params must be a mapping (install PyYAML for nested YAML)")
    if params.get("objective", "binary:logistic") != "binary:logistic":
        raise ValueError("HO XGBoost requires objective binary:logistic")
    if params.get("early_stopping_rounds") is not None or params.get("callbacks"):
        raise ValueError("HO baseline uses a fixed training fit; early stopping/callbacks are not supported")


def _dependencies():
    try:
        import xgboost as xgb
        from sklearn.metrics import roc_auc_score, roc_curve
        os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("HO xgboost mode requires xgboost, scikit-learn, and matplotlib; "
                           "install requirements-analysis.txt in the PyROOT environment") from exc
    return xgb, plt, roc_auc_score, roc_curve


def _dataset(loaded, luminosity):
    features, eligible, weights, labels, names, sources, entries, baseline_rows = ([] for _ in range(8))
    samples = {}
    for sample, rows, raw, event_ids, row_ids in loaded:
        if sample.name in samples:
            raise ValueError(f"duplicate sample identity: {sample.name}")
        if not (len(rows) == len(raw) == len(event_ids) == len(row_ids)):
            raise ValueError(f"inconsistent HO metadata row lengths: {sample.name}")
        samples[sample.name] = sample
        for row, weight, source, entry in zip(rows, raw, event_ids, row_ids):
            accepted = row["n_selected_photons"] >= 2
            values = [row[key] for key in FEATURE_NAMES]
            if accepted and any(not math.isfinite(value) or value <= -900 for value in values):
                raise ValueError(f"invalid selected diphoton features: {sample.name}, entry {entry}")
            features.append(values)
            eligible.append(accepted)
            weights.append(weight)
            labels.append(int(sample.category == "Signal"))
            names.append(sample.name)
            sources.append(source)
            entries.append(entry)
            baseline_rows.append(row)
    data = {
        "features": np.asarray(features, dtype=float), "eligible": np.asarray(eligible, dtype=bool),
        "raw_weights": np.asarray(weights, dtype=float), "labels": np.asarray(labels, dtype=int),
        "sample_ids": np.asarray(names), "source_events": np.asarray(sources, dtype=np.int64),
        "entries": np.asarray(entries, dtype=np.int64), "samples": samples, "rows": baseline_rows,
    }
    physical = np.zeros(len(weights), dtype=float)
    for name, sample in samples.items():
        factor = luminosity * 1000 * sample.cross_section_pb * sample.weight_scale / sample.normalization_sum_weight
        if not math.isfinite(factor) or factor <= 0:
            raise ValueError(f"HO physical scale must be finite and positive: {name}")
        mask = data["sample_ids"] == name
        physical[mask] = data["raw_weights"][mask] * factor
    if not np.all(np.isfinite(physical)):
        raise ValueError("nonfinite HO physical event weights")
    data["physical_weights"] = physical
    return data


def _sample_summaries(data, selected, split, luminosity, partition="test"):
    projected = project_partition_weights(data["physical_weights"], data["sample_ids"], split, partition)
    rows = []
    for name, sample in data["samples"].items():
        in_sample = data["sample_ids"] == name
        in_partition = in_sample & (split.partitions == partition)
        accepted = in_partition & selected
        fraction = split.sample_fractions[name][partition]
        summary = signed_yield_summary(projected, data["sample_ids"], data["source_events"], selected=accepted)
        selected_sum = math.fsum(data["raw_weights"][accepted])
        efficiency = selected_sum / fraction / sample.normalization_sum_weight
        rows.append({
            "sample": name, "category": sample.category, "analysis": sample.analysis_name,
            "detector_response": sample.detector_response, "response_mode": sample.response_mode,
            "input_file": str(sample.var_file), "raw_cross_section_pb": sample.cross_section_pb,
            "cross_section_pb": sample.cross_section_pb * sample.weight_scale, "weight_scale": sample.weight_scale,
            "entries_read": int(np.sum(in_sample)), "selected_entries": summary["selected_entries"],
            "sum_weight": math.fsum(data["raw_weights"][in_sample]), "sum_selected_weight": selected_sum,
            "analysis_efficiency": efficiency,
            "selected_cross_section_pb": summary["sum_weight"] / (luminosity * 1000),
            "expected_events": summary["sum_weight"], "mc_events_after_analysis": summary["selected_entries"],
            "normalization_kind": sample.normalization_kind,
            "normalization_sum_weight": sample.normalization_sum_weight,
            "ihixs_record_sha256": sample.ihixs_record_sha256, "shower_quality": sample.shower_quality,
            "evaluation_partition": partition, "partition_fraction": fraction,
            "partition_entries": int(np.sum(in_partition)),
            "selected_source_events": summary["selected_source_events"],
            "mc_variance_events": summary["variance"], "mc_error_events": summary["mc_error"],
            "effective_count": summary["effective_count"], "positive_events": summary["positive_weight"],
            "negative_events": summary["negative_weight"],
        })
    return rows


def _totals(rows, systematics, minimum):
    signal = sum(row["expected_events"] for row in rows if row["category"] == "Signal")
    background = sum(row["expected_events"] for row in rows if row["category"] != "Signal")
    svar = sum(row["mc_variance_events"] for row in rows if row["category"] == "Signal")
    bvar = sum(row["mc_variance_events"] for row in rows if row["category"] != "Signal")
    effective = background * background / bvar if bvar > 0 else 0.0
    positive = signal > 0 and background > 0
    return {
        "signal_expected_events": signal, "background_expected_events": background,
        "signal_mc_error_events": math.sqrt(svar), "background_mc_error_events": math.sqrt(bvar),
        "background_effective_count": effective,
        "approx_significance_s_over_sqrt_b": signal / math.sqrt(background) if positive else None,
        "significance": signal / math.sqrt(background + (systematics * background) ** 2) if positive else None,
        "statistically_supported": bool(positive and (effective >= minimum or math.isclose(
            effective, minimum, rel_tol=1e-12, abs_tol=0.0))),
    }


def _csv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + "\n")


def _provenance(loaded, run_tag):
    output = []
    for sample, _, _, _, _ in loaded:
        stat = sample.var_file.stat()
        records = {}
        # Small sidecars are hashed; ROOT population integrity is checked by the
        # CLI, and exact group/row assignments are persisted separately.
        for path in (sample.dat_file, sample.sample_dir / "campaign.json",
                     sample.sample_dir / f"normalization-{run_tag}.json"):
            if path.is_file():
                records[path.name] = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        output.append({
            "sample": sample.name, "input_file": str(sample.var_file.resolve()),
            "input_size": stat.st_size, "input_mtime_ns": stat.st_mtime_ns,
            "tree_entries": sample.tree_entries, "sum_tree_weight": sample.sum_tree_weight,
            "normalization_kind": sample.normalization_kind,
            "normalization_sum_weight": sample.normalization_sum_weight,
            "cross_section_pb": sample.cross_section_pb, "weight_scale": sample.weight_scale,
            "ihixs_record_sha256": sample.ihixs_record_sha256, "shower_quality": sample.shower_quality,
            "records": records,
        })
    return output


def _plots(stage, data, split, scores, model, threshold, baseline_mask, plt, roc_curve, roc_auc_score):
    test = data["eligible"] & (split.partitions == "test")
    projected = project_partition_weights(data["physical_weights"], data["sample_ids"], split, "test")
    y, s, w = data["labels"][test], scores[test], np.abs(projected[test])
    diagnostics = {"auc_unweighted": None, "auc_absolute_weight": None}
    fig, ax = plt.subplots(figsize=(6, 5))
    if len(np.unique(y)) == 2 and all(np.sum(w[y == label]) > 0 for label in (0, 1)):
        diagnostics = {"auc_unweighted": float(roc_auc_score(y, s)),
                       "auc_absolute_weight": float(roc_auc_score(y, s, sample_weight=w))}
        for label, weight in (("Unweighted", None), ("Absolute physical weights", w)):
            fpr, tpr, _ = roc_curve(y, s, sample_weight=weight)
            ax.plot(fpr, tpr, label=label)
        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
        ax.legend()
    else:
        ax.text(.5, .5, "Insufficient test classes for ROC", ha="center", transform=ax.transAxes)
    ax.set(xlabel="Background efficiency (diagnostic)", ylabel="Signal efficiency (diagnostic)",
           title="Untouched test ROC — nonnegative diagnostic weights")
    fig.tight_layout(); fig.savefig(stage / "roc.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    importance = np.asarray(model.feature_importances_)
    order = np.argsort(importance)
    ax.barh(np.asarray(FEATURE_NAMES)[order], importance[order])
    ax.set(xlabel="Feature importance", title="HO XGBoost")
    fig.tight_layout(); fig.savefig(stage / "feature_importance.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(0, 1, 31)
    for partition, style in (("train", "--"), ("test", "-")):
        partition_weights = project_partition_weights(data["physical_weights"], data["sample_ids"], split, partition)
        for label, name, color in ((1, "Signal", "tab:blue"), (0, "Background", "tab:orange")):
            mask = data["eligible"] & (split.partitions == partition) & (data["labels"] == label)
            weights = np.abs(partition_weights[mask])
            if np.sum(weights) > 0:
                ax.hist(scores[mask], bins=bins, weights=weights / np.sum(weights), histtype="step",
                        linestyle=style, color=color, label=f"{name} {partition}")
    ax.set(xlabel="Classifier score", ylabel="Fraction of absolute physical weight",
           title="Training/test diagnostic — each curve normalized")
    ax.legend(); fig.tight_layout(); fig.savefig(stage / "score_distributions.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    mass = data["features"][:, 0]
    finite_mass = mass[test]
    if finite_mass.size:
        lo, hi = float(np.min(finite_mass)), float(np.max(finite_mass))
        if lo == hi:
            lo, hi = lo - 1, hi + 1
    else:
        lo, hi = 100., 150.
    selections = [("Two photons", data["eligible"]), ("Cut baseline", baseline_mask)]
    if threshold is not None:
        selections.append(("XGBoost", data["eligible"] & (scores >= threshold)))
    for row_axes, bins, view in ((axes[0], np.linspace(lo, hi, 51), "full test range"),
                                 (axes[1], np.linspace(100, 150, 51), "100–150 GeV view")):
        for label, name, ax in zip((1, 0), ("Signal", "Background"), row_axes):
            for selection, mask in selections:
                selected = mask & test & (data["labels"] == label)
                ax.hist(mass[selected], bins=bins, weights=projected[selected], histtype="step", label=selection)
            ax.axhline(0, color="black", linewidth=.5)
            ax.set(xlabel="Diphoton mass [GeV]", ylabel="Projected signed events / bin",
                   title=f"{name}: {view}")
            ax.legend()
    fig.suptitle("Untouched test events — viewing ranges do not change the selection")
    fig.tight_layout(); fig.savefig(stage / "mass_distributions.png", dpi=160); plt.close(fig)
    return diagnostics


def _html_report(stage, metadata, comparison, assets):
    def fmt(value):
        return "undefined" if value is None else f"{value:.6g}" if isinstance(value, float) else str(value)
    def table(fields, rows):
        return ('<div class="scroll"><table><thead><tr>' + ''.join(f'<th>{html.escape(field)}</th>' for field in fields)
                + '</tr></thead><tbody>' + ''.join('<tr>' + ''.join(
                    f'<td>{html.escape(fmt(row.get(field)))}</td>' for field in fields) + '</tr>' for row in rows)
                + '</tbody></table></div>')
    totals, samples = [], []
    for selection in ("baseline", "xgboost"):
        result = comparison[selection]
        if result is not None:
            totals.append({"selection": selection, **result["totals"]})
            samples.extend({"selection": selection, **row} for row in result["samples"])
    summary = table(["selection", "signal_expected_events", "signal_mc_error_events",
                     "background_expected_events", "background_mc_error_events", "background_effective_count",
                     "significance", "statistically_supported"], totals)
    detail = table(["selection", "sample", "selected_entries", "analysis_efficiency", "expected_events",
                    "mc_error_events", "positive_events", "negative_events"], samples)
    warnings = ''.join(f'<li>{html.escape(message)}</li>' for message in metadata["warnings"])
    pictures = ''.join(f'<figure><img src="{html.escape(name)}" alt="{html.escape(name)}"></figure>'
                       for name in assets if name.endswith('.png'))
    links = ''.join(f'<li><a href="{html.escape(name)}">{html.escape(name)}</a></li>' for name in assets)
    baseline = '; '.join(f"{item['variable']}: {item['min']} to {item['max']}" for item in metadata['baseline_cuts'])
    (stage / "index.html").write_text(f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{html.escape(metadata['name'])}</title>
<style>body{{font:16px system-ui,sans-serif;margin:2rem;color:#1f2933;line-height:1.5}}
.scroll{{overflow-x:auto}}table{{border-collapse:collapse;width:100%;font-size:.9rem}}
th,td{{text-align:right;padding:.55rem;border-bottom:1px solid #d9e2ec}}th{{background:#f0f4f8}}
th:first-child,td:first-child{{text-align:left}}img{{max-width:100%;height:auto}}figure{{margin:1.5rem 0}}
.notice{{background:#fff3cd;padding:1rem}}a{{color:#1f5eff}}</style></head><body>
<h1>{html.escape(metadata['name'])}</h1>
<p>HO diphoton XGBoost • campaign {html.escape(metadata['run_tag'])} •
{metadata['luminosity_fb']:.6g} fb<sup>−1</sup> • {html.escape(metadata['detector_response'])} response</p>
<p>Status: <strong>{html.escape(metadata['status'])}</strong>. Frozen score threshold:
<strong>{fmt(metadata['best_threshold'])}</strong>.</p>
<div class="notice"><ul>{warnings}</ul></div>
<h2>Independent test comparison</h2>
<p>Both selections use the same untouched test source events, projected by each sample's actual sampling fraction.
The classifier uses diphoton mass; no mass window is applied before training. Baseline: {html.escape(baseline)}.</p>
<p>Counting estimate: S / √(B + (δB)²), δ = {metadata['systematics']:.6g}.
MC errors combine response rows from each source event before squaring and condition on the saved normalization.
They do not include cross-section or detector systematic uncertainties.</p>{summary}
<h2>Per-sample signed yields</h2>{detail}
<p>Training uses class-balanced absolute physical weights. Signs are retained for all physics yields.
ROC and score overlays use nonnegative diagnostic weights; training scores never enter the test comparison.</p>
<h2>Diagnostics</h2>{pictures}<h2>Reproducibility and downloads</h2><ul>{links}</ul>
</body></html>''')


def run_ho_signal_background_analysis(loaded_samples, *, output_dir, metadata, baseline_cuts, config):
    """Train once; tune only on validation; publish independent test results."""
    validate_config(config)
    xgb, plt, roc_auc_score, roc_curve = _dependencies()
    output_dir = Path(output_dir)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    # Publish only after all calculations and report serialization succeed.
    with tempfile.TemporaryDirectory(prefix=".ho-xgboost-", dir=output_dir.parent) as temporary:
        stage = Path(temporary)
        result = _run(loaded_samples, stage, output_dir, dict(metadata), baseline_cuts, config,
                      xgb, plt, roc_auc_score, roc_curve)
        output_dir.mkdir(parents=True, exist_ok=True)
        _publish(stage, output_dir)
    return result


def _publish(stage, output_dir):
    """Replace generated artifacts, restoring the previous report on failure."""
    paths = sorted(stage.iterdir(), key=lambda path: path.name == "index.html")
    backup = stage / "previous"
    backup.mkdir()
    saved, published = [], []
    try:
        for path in paths:
            destination = output_dir / path.name
            if destination.exists():
                if not destination.is_file():
                    raise ValueError(f"report output is not a file: {destination}")
                os.replace(destination, backup / path.name)
                saved.append(path.name)
            os.replace(path, destination)
            published.append(path.name)
    except Exception:
        for name in reversed(published):
            (output_dir / name).unlink()
        for name in saved:
            os.replace(backup / name, output_dir / name)
        raise


def _run(loaded, stage, output_dir, metadata, baseline_cuts, config, xgb, plt, roc_auc_score, roc_curve):
    luminosity = metadata["luminosity_fb"]
    validation_size, test_size = float(config.get("validation_size", .2)), float(config.get("test_size", .2))
    seed, systematics = int(config.get("seed", 12345)), float(config.get("systematics", 0))
    minimum = float(config.get("min_background_effective_count", 25))
    data = _dataset(loaded, luminosity)
    provenance = _provenance(loaded, metadata["run_tag"])
    split = grouped_sample_split(data["sample_ids"], data["source_events"],
                                 fractions=(1-validation_size-test_size, validation_size, test_size), seed=seed)
    masks = {part: data["eligible"] & (split.partitions == part) for part in ("train", "validation", "test")}
    for part in ("train", "validation"):
        mask = masks[part]
        if any(np.sum(np.abs(data["physical_weights"][mask & (data["labels"] == label)])) <= 0 for label in (0, 1)):
            raise ValueError(f"HO {part} partition requires nonzero-weight eligible signal and background")
    params = {**DEFAULT_MODEL_PARAMS, "random_state": seed, **(config.get("model_params") or {})}
    model = xgb.XGBClassifier(**params)
    projected_train = project_partition_weights(data["physical_weights"], data["sample_ids"], split, "train")
    training_weights = class_balanced_abs_weights(data["labels"][masks["train"]], projected_train[masks["train"]])
    model.fit(data["features"][masks["train"]], data["labels"][masks["train"]], sample_weight=training_weights)
    scores = np.full(len(data["labels"]), np.nan)
    scores[data["eligible"]] = model.predict_proba(data["features"][data["eligible"]])[:, 1]
    if not np.all(np.isfinite(scores[data["eligible"]])):
        raise ValueError("classifier produced nonfinite scores")
    validation_weights = project_partition_weights(data["physical_weights"], data["sample_ids"], split, "validation")
    mask = masks["validation"]
    best, scan = scan_thresholds(scores[mask], data["labels"][mask], validation_weights[mask],
                                data["sample_ids"][mask], data["source_events"][mask],
                                systematics=systematics, min_background_effective_events=minimum)
    threshold = best["threshold"] if best is not None else None
    baseline_mask = np.asarray([all(cut.accepts(row) for cut in baseline_cuts) for row in data["rows"]])
    baseline_rows = _sample_summaries(data, baseline_mask, split, luminosity)
    summary_rows = [] if threshold is None else _sample_summaries(
        data, data["eligible"] & (scores >= threshold), split, luminosity)
    comparison = {
        "evaluation_partition": "test", "baseline": {
            "samples": baseline_rows, "totals": _totals(baseline_rows, systematics, minimum)},
        "xgboost": None if threshold is None else {
            "samples": summary_rows, "totals": _totals(summary_rows, systematics, minimum)},
    }
    warnings = ["Approximate counting sensitivity; finite Monte Carlo statistics and mass selection must be inspected."]
    if threshold is None:
        warnings.append("No validation threshold satisfies positive S/B and the minimum background effective count. "
                        "No XGBoost selection or test sensitivity is reported; increase simulation statistics or review the configured guard.")
    for selection in ("baseline", "xgboost"):
        if comparison[selection] is not None and not comparison[selection]["totals"]["statistically_supported"]:
            warnings.append(f"{selection}: test sensitivity has nonpositive yields or background effective count below {minimum:g}; "
                            "the frozen threshold is not retuned on test data.")
    metadata.update({
        "status": "insufficient_statistics" if threshold is None else "complete",
        "evaluation_partition": "test", "feature_names": FEATURE_NAMES, "best_threshold": threshold,
        "validation_best": best, "systematics": systematics, "min_background_effective_count": minimum,
        "seed": seed, "requested_split_fractions": {"train": 1-validation_size-test_size,
                                                     "validation": validation_size, "test": test_size},
        "sample_split_fractions": split.sample_fractions, "sample_group_counts": split.sample_group_counts,
        "eligible_rows_per_partition": {part: int(mask.sum()) for part, mask in masks.items()},
        "model_params": params, "training_weight_policy": "class_balanced_absolute_physical_weights",
        "mc_variance_policy": "sum of squared selected source-event contributions; fixed source normalization",
        "normalization_provenance": provenance, "warnings": warnings,
        "software_versions": {"python": platform.python_version(), **{
            name: importlib.metadata.version(name) for name in ("numpy", "xgboost", "scikit-learn", "matplotlib")}},
        "totals": None if threshold is None else comparison["xgboost"]["totals"],
        "comparison": comparison,
    })
    model.save_model(str(stage / "signal_background_xgboost.json"))
    metadata["diagnostics"] = _plots(stage, data, split, scores, model, threshold, baseline_mask, plt,
                                      roc_curve, roc_auc_score)
    _csv(stage / "summary.csv", SUMMARY_FIELDS, summary_rows)
    _csv(stage / "comparison.csv", ["selection", *SUMMARY_FIELDS],
         ({"selection": name, **row} for name, rows in (("baseline", baseline_rows), ("xgboost", summary_rows)) for row in rows))
    _json(stage / "comparison.json", comparison)
    _csv(stage / "threshold_scan.csv", list(scan[0]), scan)
    _csv(stage / "partitions.csv", ["sample", "sourceevent", "tree_entry", "partition", "eligible"], (
        {"sample": data["sample_ids"][i], "sourceevent": int(data["source_events"][i]),
         "tree_entry": int(data["entries"][i]), "partition": split.partitions[i], "eligible": bool(data["eligible"][i])}
        for i in range(len(scores))))
    _csv(stage / "scores.csv", ["sample", "source", "sourceevent", "tree_entry", "partition", "label", "score",
                                 "raw_signed_weight", "physical_weight", "partition_fraction", "projected_physical_weight"], (
        {"sample": data["sample_ids"][i], "source": str(data["samples"][data["sample_ids"][i]].var_file),
         "sourceevent": int(data["source_events"][i]), "tree_entry": int(data["entries"][i]),
         "partition": split.partitions[i], "label": int(data["labels"][i]), "score": float(scores[i]),
         "raw_signed_weight": float(data["raw_weights"][i]), "physical_weight": float(data["physical_weights"][i]),
         "partition_fraction": split.sample_fractions[data["sample_ids"][i]][split.partitions[i]],
         "projected_physical_weight": float(data["physical_weights"][i]) / split.sample_fractions[data["sample_ids"][i]][split.partitions[i]]}
        for i in np.flatnonzero(data["eligible"])))
    assets = ["summary.csv", "summary.json", "comparison.csv", "comparison.json", "metrics.json",
              "signal_background_xgboost.json", "model_metadata.json", "scores.csv", "partitions.csv",
              "threshold_scan.csv", "roc.png", "feature_importance.png", "score_distributions.png", "mass_distributions.png"]
    metadata["outputs"] = {name: str(output_dir / name) for name in assets}
    _json(stage / "metrics.json", metadata)
    _json(stage / "model_metadata.json", {key: value for key, value in metadata.items() if key not in ("comparison", "outputs")})
    _json(stage / "summary.json", {"metadata": metadata, "samples": summary_rows})
    _html_report(stage, metadata, comparison, assets)
    if _provenance(loaded, metadata["run_tag"]) != provenance:
        raise ValueError("HO input or normalization provenance changed during training; report not published")
    return {"metadata": metadata, "summary_rows": summary_rows,
            "assets": [output_dir / name for name in assets if name not in ("summary.csv", "summary.json")]}
