#!/usr/bin/env python3
"""Cut-based and XGBoost analyses of audited four-lepton manifests.

Only the user invocation trains a model. Signed yields use full source-event
normalization. Training uses absolute weights, with shared production events
and all their detector hypotheses assigned to the same partition.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
import numpy as np
from hfourlepton import make_fourlepton_report as report
from ho_xgboost_core import (grouped_sample_split, project_partition_weights,
    signed_yield_summary, scan_thresholds, class_balanced_abs_weights)

FEATURES = ("mZ1", "mZ2", "pt4l", "abs_y4l", "pt1", "pt2", "pt3", "pt4", "min_delta_r")
VARIABLES = (*FEATURES, "m4l")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_card(path):
    text = Path(path).read_text()
    if Path(path).suffix == ".json":
        card = json.loads(text)
    else:
        import yaml
        card = yaml.safe_load(text)
    if not isinstance(card, dict):
        raise ValueError("analysis card must be a mapping")
    for section in ("preselection", "cuts"):
        for name, limits in card.get(section, {}).items():
            if name not in VARIABLES or not isinstance(limits, dict) or set(limits) - {"min", "max"}:
                raise ValueError(f"invalid cut {section}.{name}")
            if any(not math.isfinite(float(v)) for v in limits.values()):
                raise ValueError("cut boundaries must be finite")
            if limits.get("min", -math.inf) >= limits.get("max", math.inf):
                raise ValueError("cut minimum must be below maximum")
    features = card.get("features", list(FEATURES))
    if not features or len(set(features)) != len(features) or set(features) - set(FEATURES):
        raise ValueError("features must be unique supported kinematics; m4l is excluded")
    return card


def cut_mask(values, cuts):
    mask = np.ones(len(values["m4l"]), dtype=bool)
    for name, limits in cuts.items():
        if "min" in limits:
            mask &= values[name] > limits["min"]
        if "max" in limits:
            mask &= values[name] < limits["max"]
    return mask


def audit_disjoint_slices(samples):
    ranges = {}
    for sample in samples:
        group = sample.raw.get("analysis_source_group", sample.production_group_id)
        if sample.category != "signal" or group == sample.production_group_id:
            continue
        proof = sample.raw.get("hard_event_slice", {})
        identity, artifact = proof.get("identity", {}), proof.get("artifact", {})
        if (identity.get("kind") != "disjoint_hard_event_slice" or identity.get("sample") != sample.name
                or artifact.get("sha256") != sample.raw.get("source_lhe_sha256")
                or artifact.get("path") != sample.raw.get("source_lhe") or group != sample.name):
            raise ValueError("independent decay channels require hard-event slice provenance")
        start, stop = identity["start"], identity["stop"]
        if not 0 <= start < stop or stop - start != sample.raw["nevents_requested"]:
            raise ValueError("invalid hard-event slice range")
        key = identity["source_sha256"]
        for other_start, other_stop in ranges.get(key, []):
            if start < other_stop and other_start < stop:
                raise ValueError("hard-event slices overlap across decay channels")
        ranges.setdefault(key, []).append((start, stop))


def candidate_row(row, name, group, category, scale, luminosity):
    pts = sorted((float(v) for v in row.get("lepton_pt", [])), reverse=True)
    if len(pts) != 4:
        raise ValueError("lepton_pt[4] is required; rerun the four-lepton analyzer")
    values = {key: float(row.get(key, 0.)) for key in ("m4l", "mZ1", "mZ2", "pt4l", "min_delta_r")}
    values.update(abs_y4l=abs(float(row.get("y4l", 0.))), **{f"pt{i+1}": pt for i, pt in enumerate(pts)})
    if not all(math.isfinite(v) for v in values.values()):
        raise ValueError("nonfinite analysis kinematics")
    return dict(sample=name, group=group, source=int(row["source_index"]), label=int(category == "signal"),
                weight=float(row["event_weight"]) * scale * luminosity * 1000., selected=True, **values)


def load_dataset(manifest, card, luminosity, *, for_ml=False):
    bundle = report.load_manifest(Path(manifest))
    errors = list(bundle.audit.errors)
    raw_manifest = json.loads(Path(manifest).read_text())
    if bundle.detector_profile != card["detector_response"]:
        raise ValueError("analysis card and manifest detector responses differ")
    if raw_manifest.get("selection_profile", "gem") != card["selection_profile"]:
        raise ValueError("selection profile differs; rerun reconstruction with the intended profile")
    names = {s.name for s in bundle.samples}
    from hfourlepton import pooled_controls
    requested = set(card["required_samples"])
    extra_names = set(pooled_controls.PROCESSES)
    use_published = extra_names <= requested and names.isdisjoint(extra_names)
    controls = (pooled_controls.load(manifest, bundle.cut_stages, published_defaults=use_published)
                if "pooled_backgrounds" in raw_manifest or use_published else [])
    if names.intersection(s["name"] for s, _, _ in controls):
        raise ValueError("a pooled process overlaps the manifest sample inventory")
    names.update(s["name"] for s, _, _ in controls)
    missing = requested - names
    if missing:
        raise ValueError(f"missing required samples: {sorted(missing)}")
    # A nominal card must say explicitly whether reducible controls belong in the sum.
    if names - requested:
        raise ValueError("manifest contains samples outside this card's declared sample inventory")
    audit_disjoint_slices(bundle.samples)
    versions = {s.raw.get("selection_version", "legacy_unversioned") for s in bundle.samples}
    if len(versions) != 1:
        raise ValueError("samples use different selection versions; reprocess them together")
    multiplicities = {}
    for sample in bundle.samples:
        group = sample.raw.get("analysis_source_group", sample.production_group_id) or sample.name
        multiplicities[group] = multiplicities.get(group, 0) + 1
    rows, metadata = [], []
    populations = {}
    for sample in bundle.samples:
        if sample.raw.get("status") != "complete":
            raise ValueError(f"incomplete sample: {sample.name}")
        for key, path in (("output_root_sha256", sample.output_root),
                          ("analysis_summary_sha256", Path(sample.raw["analysis_summary"]))):
            if not path or not path.is_file() or digest(path) != sample.raw.get(key):
                raise ValueError(f"missing or changed analysis product for {sample.name}: {key}")
        from hfourlepton.report_reader import load_report_records
        print(f"Auditing {sample.name}", flush=True)
        records = load_report_records(sample.output_root,
            selected_mask=report._final_stage(bundle.cut_stages).required_mask, progress=True)
        result = report.analyze_sample(sample, records, bundle.cut_stages, luminosity,
                                      (120., 130.), 1.e-9)
        errors.extend(f"{sample.name}: {e}" for e in result["audit"]["errors"])
        scale = result["normalization"].get("scale_pb_per_weight")
        if scale is None or scale <= 0:
            raise ValueError(f"invalid signed normalization: {sample.name}")
        group = sample.raw.get("analysis_source_group", sample.production_group_id) or sample.name
        # Without a saved LHE->HwSim event map, lost shower events can shift the
        # numeric source index independently between shared decay channels.
        if multiplicities[group] > 1 and sample.raw.get("shower_completion", {}).get("discarded_events", 0):
            raise ValueError("shared hard events with shower losses lack a common event map; use the HO campaign's disjoint decay inputs")
        sources = {int(s["source_index"]) for s in records.source_events}
        if group in populations and populations[group] != sources:
            raise ValueError("shared production channels have different source populations; cannot establish independent splits")
        populations[group] = sources
        meta = dict(name=sample.name, category=sample.category, group=group,
                    production_mode=sample.raw.get("production_mode", "ggf" if sample.category == "signal" else ""),
                    target_cross_section_pb=sample.target_cross_section_pb,
                    perturbative_order=sample.perturbative_order,
                    normalization=result["normalization"], root_sha256=digest(sample.output_root))
        metadata.append(meta)
        # Dummy zero-weight rows include every unselected source in split fractions.
        for source in sorted(sources):
            rows.append(dict(sample=sample.name, group=group, source=source,
                             label=int(sample.category == "signal"), weight=0., selected=False,
                             **{name: 0. for name in VARIABLES}))
        for row in records.four_lepton:
            if int(row["source_index"]) not in sources:
                raise ValueError("hypothesis has no SourceEvents entry")
            selected = report.passes_stage(row, bundle.cut_stages[-1])
            if not selected:
                # The zero-weight source row above already accounts for this
                # event in independent train/validation/test populations.
                continue
            rows.append(candidate_row(row, sample.name, group, sample.category, scale, luminosity))
        del records
    if errors:
        raise ValueError("four-lepton audit failed:\n" + "\n".join(errors))
    for sample, path, checksum in controls:
        name = sample["name"]
        norm = sample["normalization"]
        metadata.append(dict(name=name, category="reducible", group=name, production_mode="",
            target_cross_section_pb=sample["target_pb"], perturbative_order="LO", normalization=norm,
            pooled_input=str(path), pooled_input_sha256=checksum, selected_hypotheses=len(sample["selected"])))
        rows.extend(candidate_row(row, name, name, "reducible", norm["scale_pb_per_weight"], luminosity)
                    for row in sample["selected"])
    data = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    del rows
    # For ML every independent source, including unselected sources, contributes
    # to the exact per-sample split fractions. Store the dummy rows as NumPy
    # arrays to avoid constructing millions of identical Python dictionaries.
    if controls:
        blocks = {key: [col] for key, col in data.items()}
        for sample, _, _ in controls:
            name = sample["name"]
            ranges = sample["source_ranges"]
            sources = (np.concatenate([np.arange(low, high, dtype=np.int64) for low, high in ranges])
                       if for_ml else np.asarray([ranges[0][0]], dtype=np.int64))
            for key, col in data.items():
                if key in ("group", "sample"):
                    values = np.full(len(sources), name, dtype=f"<U{len(name)}")
                elif key == "source":
                    values = sources
                else:
                    values = np.zeros(len(sources), dtype=col.dtype)
                blocks[key].append(values)
        data = {key: np.concatenate(values) for key, values in blocks.items()}
    data["eligible"] = data["selected"] & cut_mask(data, card.get("preselection", {}))
    data["baseline"] = data["eligible"] & cut_mask(data, card.get("cuts", {}))
    return data, metadata, bundle


def selected_yield(data, weights, mask):
    # Only selected hypotheses contribute. Building groups for every generated
    # source for each individual process would repeat millions of empty groups.
    return signed_yield_summary(weights[mask], data["group"][mask], data["source"][mask])


def yields(data, weights, mask):
    return {name: selected_yield(data, weights, mask & select) for name, select in
            (("signal", data["label"] == 1), ("background", data["label"] == 0))}


def totals(data, weights, mask, card):
    result = yields(data, weights, mask)
    s, b = result["signal"]["sum_weight"], result["background"]["sum_weight"]
    effective = result["background"]["effective_count"]
    minimum = float(card.get("min_background_effective_events", 25.))
    systematics = float(card.get("background_fractional_uncertainty", 0.))
    if minimum <= 0 or systematics < 0 or not math.isfinite(minimum + systematics):
        raise ValueError("invalid significance statistics settings")
    supported = s > 0 and b > 0 and effective >= minimum
    result.update(statistically_supported=supported,
        counting_diagnostic=s / math.sqrt(b + (systematics*b)**2) if supported else None,
        background_fractional_uncertainty=systematics,
        interpretation="counting sensitivity only; no experimental likelihood or discovery claim")
    return result


def sample_yields(data, weights, mask, metadata):
    return [{"sample": item["name"], "category": item["category"],
             **selected_yield(data, weights, mask & (data["sample"] == item["name"]))} for item in metadata]


def write_plots(data, weights, selection, output, metadata, *, scores=None, split=None, threshold=None, importance=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from hfourlepton.observables import COLORS, component

    background_samples = {}
    for sample in metadata:
        if sample["category"] != "signal":
            background_samples.setdefault(component(sample), []).append(sample["name"])

    mass_edges = np.linspace(105., 140., 36)
    stages = [("Preselection", selection & data["eligible"])]
    if scores is None:
        stages.append(("Mass window", selection & data["baseline"]))
    elif threshold is not None:
        stages.append(("Preselection + frozen score (before mass window)", stages[0][1] & (scores >= threshold)))

    fig, axes = plt.subplots(1, len(stages), figsize=(6 * len(stages), 4.8), squeeze=False)
    for ax, (stage, mask) in zip(axes[0], stages):
        masses, names = data["m4l"][mask], data["sample"][mask]
        w, is_signal = weights[mask], data["label"][mask] == 1
        for label, subset, color in (("Signal", is_signal, COLORS["Signal"]),
                                     ("Total background", ~is_signal, "black")):
            ax.stairs(np.histogram(masses[subset], mass_edges, weights=w[subset])[0], mass_edges,
                      color=color, label=label, linewidth=2)
        for label, samples in background_samples.items():
            part = np.isin(names, samples)
            values = np.histogram(masses[part], mass_edges, weights=w[part])[0]
            ax.stairs(values, mass_edges, color=COLORS[label], linestyle="--",
                      label=label if part.any() else f"{label} (0 selected)")
        ax.set(xlabel=r"$m_{4\ell}$ [GeV]", ylabel="Signed expected events / GeV", title=stage)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output / "mass.svg")
    plt.close(fig)

    if scores is not None:
        edges = np.linspace(0., 1., 26)
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True, sharey=True)
        for ax, partition in zip(axes, ("train", "test")):
            mask = data["eligible"] & (split.partitions == partition)
            w = np.abs(data["weight"][mask])
            names, score, is_signal = data["sample"][mask], scores[mask], data["label"][mask] == 1
            for label, subset, color in (("Signal", is_signal, COLORS["Signal"]),
                                         ("Total background", ~is_signal, "black")):
                if w[subset].sum() > 0:
                    values = np.histogram(score[subset], edges, weights=w[subset])[0]
                    ax.stairs(values / w[subset].sum() / np.diff(edges), edges, color=color,
                              label=label, linewidth=2)
            # Each background component integrates to its share of the absolute background weight.
            total = w[~is_signal].sum()
            for label, samples in background_samples.items():
                part = np.isin(names, samples)
                if total > 0:
                    values = np.histogram(score[part], edges, weights=w[part])[0]
                    ax.stairs(values / total / np.diff(edges), edges, color=COLORS[label], linestyle="--",
                              label=label if part.any() else f"{label} (0 selected)")
                else:
                    ax.plot([], [], color=COLORS[label], label=f"{label} (no background weight)")
            ax.set(xlabel="Classifier score", ylabel="Absolute-weight density", title=partition.capitalize())
            ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(output / "scores.svg")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.barh(list(importance), list(importance.values()))
        ax.set(xlabel="XGBoost feature importance")
        fig.tight_layout(); fig.savefig(output / "features.svg"); plt.close(fig)


def run_analysis(data, metadata, card, mode, output):
    cuts = {"totals": totals(data, data["weight"], data["baseline"], card),
            "samples": sample_yields(data, data["weight"], data["baseline"], metadata)}
    result = {"cuts_full_sample": cuts}
    if mode == "cuts":
        write_plots(data, data["weight"], np.ones(len(data["weight"]), bool), output, metadata)
        return result
    from xgboost import XGBClassifier
    from sklearn.metrics import roc_auc_score
    split = grouped_sample_split(data["group"], data["source"],
        fractions=card.get("split_fractions", [.6, .2, .2]), seed=card.get("seed", 12345))
    features = card.get("features", list(FEATURES))
    matrix = np.column_stack([data[name] for name in features])
    train = data["eligible"] & (split.partitions == "train")
    train_weights = class_balanced_abs_weights(data["label"][train], data["weight"][train])
    params = {"n_estimators": 200, "max_depth": 3, "learning_rate": .05,
              "subsample": .8, "colsample_bytree": 1., "objective": "binary:logistic",
              "eval_metric": "logloss", "tree_method": "hist", "n_jobs": 1,
              "random_state": card.get("seed", 12345)}
    allowed = {"n_estimators", "max_depth", "learning_rate", "subsample", "colsample_bytree", "min_child_weight", "reg_lambda", "n_jobs"}
    if set(card.get("xgboost", {})) - allowed:
        raise ValueError("unsupported XGBoost parameter; the binary objective and split seed are fixed")
    params.update(card.get("xgboost", {}))
    model = XGBClassifier(**params)
    model.fit(matrix[train], data["label"][train], sample_weight=train_weights)
    scores = np.zeros(len(matrix))
    scores[data["eligible"]] = model.predict_proba(matrix[data["eligible"]])[:, 1]
    val = data["baseline"] & (split.partitions == "validation")
    vw = project_partition_weights(data["weight"], data["group"], split, "validation")
    best, scan = scan_thresholds(scores[val], data["label"][val], vw[val],
        data["group"][val], data["source"][val],
        systematics=card.get("background_fractional_uncertainty", 0.),
        min_background_effective_events=card.get("min_background_effective_events", 25.))
    tw = project_partition_weights(data["weight"], data["group"], split, "test")
    test = split.partitions == "test"
    result.update(features=features, split_group_counts=split.sample_group_counts,
        split_fractions=split.sample_fractions, threshold=best, validation_scan=scan,
        cuts_test=totals(data, tw, data["baseline"] & test, card), xgboost_test=None,
        training_weight_policy="class-balanced absolute physical weights; signed weights for yields")
    if best is not None:
        selected = data["baseline"] & test & (scores >= best["threshold"])
        result["xgboost_test"] = {"totals": totals(data, tw, selected, card),
                                  "samples": sample_yields(data, tw, selected, metadata)}
    result["auc_absolute_weight_diagnostic"] = {}
    for partition in ("train", "validation", "test"):
        mask = data["eligible"] & (split.partitions == partition)
        y = data["label"][mask]
        w = np.abs(data["weight"][mask])
        result["auc_absolute_weight_diagnostic"][partition] = (
            float(roc_auc_score(y, scores[mask], sample_weight=w))
            if set(y) == {0, 1} and all(w[y == label].sum() > 0 for label in (0, 1)) else None)
    model.save_model(output / "model.json")
    result["feature_importance"] = dict(zip(features, (float(x) for x in model.feature_importances_)))
    # Persist frozen assignments/scores for overtraining and mass-sculpting inspection.
    with (output / "scores.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "source_group", "source", "partition", "m4l", "weight", "score", "baseline"])
        for i in np.flatnonzero(data["eligible"]):
            writer.writerow([data["sample"][i], data["group"][i], int(data["source"][i]), split.partitions[i],
                             data["m4l"][i], data["weight"][i], scores[i], bool(data["baseline"][i])])
    # Signed mass spectra show sculpting even though the mass is not a feature.
    edges = np.linspace(105., 140., 36)
    spectra = {}
    for label in (0, 1):
        before = test & data["eligible"] & (data["label"] == label)
        after = before & (scores >= best["threshold"]) if best else np.zeros(len(test), bool)
        spectra[str(label)] = {"before": np.histogram(data["m4l"][before], edges, weights=tw[before])[0].tolist(),
                               "after": np.histogram(data["m4l"][after], edges, weights=tw[after])[0].tolist()}
    result["test_mass_spectra"] = {"edges_gev": edges.tolist(), "spectra": spectra}
    write_plots(data, tw, test, output, metadata, scores=scores, split=split,
                threshold=best["threshold"] if best else None, importance=result["feature_importance"])
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("cuts", "xgboost"))
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--luminosity-fb", type=float, default=100.)
    args = p.parse_args(argv)
    if not math.isfinite(args.luminosity_fb) or args.luminosity_fb <= 0:
        p.error("luminosity must be positive and finite")
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory is nonempty; choose a new output directory")
    card = read_card(args.config)
    data, metadata, bundle = load_dataset(args.manifest, card, args.luminosity_fb, for_ml=args.mode == "xgboost")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = run_analysis(data, metadata, card, args.mode, args.output_dir)
    result.update(schema_version=1, mode=args.mode, card=card, luminosity_fb=args.luminosity_fb,
                  manifest=str(args.manifest.resolve()), manifest_sha256=digest(args.manifest),
                  card_sha256=digest(args.config), samples=metadata)
    (args.output_dir / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    rows = result["cuts_full_sample"]["samples"]
    with (args.output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    brief = {k: v for k, v in result.items() if k not in ("validation_scan", "samples", "test_mass_spectra")}
    gallery = ''.join(f'<p><img src="{path.name}" alt="{path.stem}" style="max-width:100%"></p>' for path in sorted(args.output_dir.glob("*.svg")))
    (args.output_dir / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>Four-lepton analysis</title>'
        '<h1>Four-lepton ' + args.mode + '</h1><p>Signed physical yields. Counting sensitivity is not an experimental discovery test.</p>'
        '<p><a href="summary.json">Full results</a> · <a href="summary.csv">Sample yields</a></p>' + gallery + '<pre>' +
        html.escape(json.dumps(brief, indent=2)) + '</pre>')
    print(args.output_dir / "index.html")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error))
