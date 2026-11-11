#!/usr/bin/env python3
"""Extend a frozen classifier test with independent, completed qq samples.

No fitting, threshold scan, generation, or changes to the original analysis.
Independent cross-section estimates are averaged by their preselection hard
event exposure, so adding another file does not add a second physical process.
"""
from __future__ import annotations
import argparse
import copy
import csv
from datetime import datetime, timezone
import fcntl
import html
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import analyze_fourlepton as analysis, reprocess as revision
from hfourlepton import report_reader
from hfourlepton.replace_qq_backgrounds import PHYSICS_FIELDS, QQ
import ho_xgboost_core


def frozen_analysis(directory):
    completion = json.loads((directory / "completion.json").read_text())
    if completion.get("result") != "complete":
        raise ValueError("base XGBoost analysis is incomplete")
    for filename, checksum in completion["files"].items():
        path = directory / filename
        if not path.resolve().is_relative_to(directory.resolve()):
            raise ValueError("invalid completed analysis artifact")
        revision.verify(path, checksum)
    result = json.loads((directory / "summary.json").read_text())
    if result.get("mode") != "xgboost" or not result.get("threshold") or not result["threshold"]["valid"]:
        raise ValueError("a validation-selected frozen classifier threshold is required")
    revision.verify(result["manifest"], result["manifest_sha256"])
    base = json.loads(Path(result["manifest"]).read_text())
    if (base.get("status") != "complete" or result["card"]["detector_response"] != "ssc"
            or result["card"]["selection_profile"] != "ssc"):
        raise ValueError("this test extension requires a completed SSC analysis")
    return result, base


def compatible_subset(base, extra, summary):
    for key in ("collider", "selection_profile", "detector_profile", "detector", "order", "generation_cuts", "cut_stages"):
        if base.get(key) != extra.get(key):
            raise ValueError(f"independent test physics differs: {key}")
    old = {s["name"]: s for s in base["samples"]}
    new = {s["name"]: s for s in extra["samples"] if s["name"] in QQ}
    if len(old) != len(base["samples"]) or set(new) != QQ or sum(s["name"] in QQ for s in extra["samples"]) != len(QQ):
        raise ValueError("test extension requires exactly one sample for each qq channel")
    training = {s["name"]: s for s in summary["samples"]}
    for name in sorted(QQ):
        left, right = old[name], new[name]
        parent = base
        if isinstance(left.get("source_manifest"), dict):
            source = left["source_manifest"]
            revision.verify(source["path"], source["sha256"])
            parent = json.loads(Path(source["path"]).read_text())
        if (parent.get("sampling", {}).get("continuum_m4l_gev")
                != extra.get("sampling", {}).get("continuum_m4l_gev")):
            raise ValueError(f"continuum generation window changed: {name}")
        for key in PHYSICS_FIELDS:
            if left.get(key) != right.get(key):
                raise ValueError(f"independent sample changes {name} {key}")
        if right.get("status") != "complete":
            raise ValueError(f"incomplete additional sample: {name}")
        for key in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
            if not right.get(key) or right[key] == left.get(key):
                raise ValueError(f"additional sample overlaps the current model data: {name} {key}")
        if (not right.get("seeds", {}).get("production")
                or right["seeds"]["production"] == left.get("seeds", {}).get("production")):
            raise ValueError(f"additional sample reuses the production seed: {name}")
        if training[name]["root_sha256"] != left["output_root_sha256"]:
            raise ValueError(f"current training provenance differs from its manifest: {name}")
    subset = copy.deepcopy(extra)
    subset["samples"] = [copy.deepcopy(new[name]) for name in sorted(QQ)]
    return subset


def exposure_coefficients(summary, subset):
    original = {s["name"]: s for s in summary["samples"]}
    result = {}
    for sample in subset["samples"]:
        name = sample["name"]
        old = original[name]
        fraction = summary["split_fractions"][old["group"]]["test"]
        old_exposure = old["normalization"]["denominator_event_count"] * fraction
        new_exposure = sample["normalization"]["denominator_event_count"]
        if not 0 < fraction < 1 or not all(math.isfinite(n) and n > 0 for n in (old_exposure, new_exposure)):
            raise ValueError(f"invalid hard-event exposure: {name}")
        result[name] = dict(original_test_exposure=old_exposure, additional_test_exposure=new_exposure,
                            original_fraction=old_exposure / (old_exposure + new_exposure),
                            additional_fraction=new_exposure / (old_exposure + new_exposure))
    return result


def original_test(directory, summary):
    metadata = {s["name"]: s for s in summary["samples"]}
    rows = []
    with (directory / "scores.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["partition"] != "test":
                continue
            sample = row["sample"]
            group = row["source_group"]
            if metadata[sample]["group"] != group or row["baseline"] not in ("True", "False"):
                raise ValueError("invalid frozen score assignments")
            rows.append(dict(sample=sample, group="original/" + group, source=int(row["source"]),
                label=int(metadata[sample]["category"] == "signal"), m4l=float(row["m4l"]),
                weight=float(row["weight"]) / summary["split_fractions"][group]["test"],
                score=float(row["score"]), baseline=row["baseline"] == "True", eligible=True))
    data = {key: np.asarray([row[key] for row in rows]) for key in rows[0]}
    if any(not np.all(np.isfinite(data[key])) for key in ("weight", "score", "m4l")):
        raise ValueError("nonfinite frozen test values")
    selected = data["baseline"] & (data["score"] >= summary["threshold"]["threshold"])
    actual = analysis.totals(data, data["weight"], selected, summary["card"])
    for kind in ("signal", "background"):
        for key in ("sum_weight", "variance", "effective_count"):
            if not math.isclose(actual[kind][key], summary["xgboost_test"]["totals"][kind][key],
                                rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"frozen test scores do not reproduce the original {kind} {key}")
    return data


def pool(original, extra, coefficients):
    """One cross-section estimate per channel; disjoint source namespaces."""
    old = {key: value.copy() for key, value in original.items()}
    new = {key: value.copy() for key, value in extra.items()}
    if set(new["sample"]) - set(coefficients):
        raise ValueError("additional test contains an unexpected process")
    for name, values in coefficients.items():
        old["weight"][old["sample"] == name] *= values["original_fraction"]
        new["weight"][new["sample"] == name] *= values["additional_fraction"]
    if set(old["group"]) & set(new["group"]):
        raise ValueError("test source namespaces overlap")
    return {key: np.concatenate([old[key], new[key]]) for key in old}


def prepare(base_analysis, extra_manifest, output):
    if output.exists() and any(output.iterdir()):
        raise ValueError("test-extension output already exists; use evaluate to resume or choose a new directory")
    summary, base = frozen_analysis(base_analysis)
    extra = json.loads(extra_manifest.read_text())
    subset = compatible_subset(base, extra, summary)
    coefficients = exposure_coefficients(summary, subset)
    original_test(base_analysis, summary)
    plan = dict(schema_version=1, prepared_at_utc=datetime.now(timezone.utc).isoformat(),
        base_analysis=str(base_analysis), extra_manifest=revision.artifact(extra_manifest),
        base_manifest=revision.artifact(Path(summary["manifest"])),
        frozen_files={name: revision.artifact(base_analysis / name) for name in
                      ("completion.json", "summary.json", "model.json", "scores.csv")},
        threshold=summary["threshold"]["threshold"], coefficients=coefficients,
        additional_usage="all additional qq events are test-only; no fitting or threshold scan",
        prior_usage="the additional events may have appeared in an earlier development campaign; they are independent of the current model training and validation",
        software={name: revision.artifact(Path(module.__file__).resolve()) for name, module in
                  (("analysis", analysis), ("revision", revision), ("statistics", ho_xgboost_core),
                   ("report", analysis.report), ("reader", report_reader))})
    plan["software"]["extension"] = revision.artifact(Path(__file__).resolve())
    subset["analysis_revision"] = {"kind": "independent_frozen_test_subset", "parent": plan["extra_manifest"]}
    output.mkdir(parents=True, exist_ok=True)
    revision.save(output / "manifest_additional_ssc.json", subset)
    plan["additional_manifest"] = revision.artifact(output / "manifest_additional_ssc.json")
    revision.save(output / "plan.json", plan)
    print(f"Prepared independent test extension: {output / 'plan.json'}", flush=True)


def evaluate(output):
    with (output / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan = json.loads((output / "plan.json").read_text())
        for record in (plan["extra_manifest"], plan["base_manifest"], plan["additional_manifest"],
                       *plan["frozen_files"].values(), *plan["software"].values()):
            revision.verify(record["path"], record["sha256"])
        completed = output / "completion.json"
        if completed.exists():
            done = json.loads(completed.read_text())
            revision.verify(output / "plan.json", done["plan_sha256"])
            for filename, checksum in done["files"].items():
                revision.verify(output / filename, checksum)
            print("Already complete:", output / "index.html")
            return
        directory = Path(plan["base_analysis"])
        summary, base = frozen_analysis(directory)
        original = original_test(directory, summary)
        subset_path = Path(plan["additional_manifest"]["path"])
        subset = json.loads(subset_path.read_text())
        compatible_subset(base, subset, summary)
        if exposure_coefficients(summary, subset) != plan["coefficients"] or plan["threshold"] != summary["threshold"]["threshold"]:
            raise ValueError("fixed test plan changed")
        card = copy.deepcopy(summary["card"])
        card["required_samples"] = sorted(QQ)
        data, metadata, _ = analysis.load_dataset(subset_path, card, summary["luminosity_fb"], for_ml=True)
        from xgboost import XGBClassifier
        model = XGBClassifier()
        model.load_model(directory / "model.json")
        eligible = data["eligible"]
        features = summary["features"]
        if features != summary["card"]["features"]:
            raise ValueError("frozen model feature order differs from its card")
        scores = model.predict_proba(np.column_stack([data[name][eligible] for name in features]))[:, 1]
        extra = {key: data[key][eligible] for key in original if key not in ("group", "score")}
        extra["group"] = np.asarray(["additional/" + str(g) for g in data["group"][eligible]])
        extra["score"] = scores
        combined = pool(original, extra, plan["coefficients"])
        selected = combined["baseline"] & (combined["score"] >= plan["threshold"])
        result = dict(schema_version=1, threshold=plan["threshold"], model=plan["frozen_files"]["model.json"],
            coefficients=plan["coefficients"], original_test=summary["xgboost_test"],
            extended_test={"totals": analysis.totals(combined, combined["weight"], selected, card),
                "samples": analysis.sample_yields(combined, combined["weight"], selected, summary["samples"])},
            additional_test={"background": analysis.yields(extra, extra["weight"],
                extra["baseline"] & (scores >= plan["threshold"]))["background"], "normalization": metadata},
            note="The classifier and validation-selected cut are fixed. Additional qq events have no training or validation assignment. Original signal and gg test events are unchanged.")
        revision.save(output / "summary.json", result)
        with (output / "additional_scores.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(list(extra))
            writer.writerows(zip(*(extra[k] for k in extra)))
        # A dedicated diagnostic shows the fixed cut without retraining or a
        # misleading new train/test ROC comparison.
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        for ax, label, title in zip(axes, (1, 0), ("Signal", "Background")):
            mask = combined["label"] == label
            edges = np.linspace(105., 140., 36)
            for selection, text in ((mask, "Preselection"), (mask & selected, "Frozen score cut + mass window")):
                ax.stairs(np.histogram(combined["m4l"][selection], edges, weights=combined["weight"][selection])[0], edges, label=text)
            ax.set(xlabel=r"$m_{4\ell}$ [GeV]", ylabel="Signed expected events / GeV", title=title)
            ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output / "mass.svg")
        plt.close(fig)
        totals = result["extended_test"]["totals"]
        state = "supported" if totals["statistically_supported"] else "more_statistics_needed"
        (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>Extended frozen-model test</title>'
            '<h1>Extended frozen-model test</h1><p>' + html.escape(result["note"]) + '</p>'
            '<p>Earlier development-campaign qq events are independent of the current model training and validation.</p>'
            '<p>Independent normalized qq rate estimates are averaged by hard-event exposure; each process contributes one rate.</p>'
            f'<p>Statistical support: {state}</p><p><a href="summary.json">Full results</a> · <a href="plan.json">Frozen plan</a></p>'
            '<img src="mass.svg" style="max-width:100%" alt="Frozen classifier test mass spectra">'
            '<pre>' + html.escape(json.dumps(totals, indent=2)) + '</pre>')
        revision.save(completed, dict(state=state, plan_sha256=revision.digest(output / "plan.json"),
            files={name: revision.digest(output / name) for name in
                   ("summary.json", "additional_scores.csv", "mass.svg", "index.html")}))
        print(f"Extended test: {state}; background effective events = {totals['background']['effective_count']:.2f}", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--base-analysis", required=True, type=Path)
    p.add_argument("--additional-manifest", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p = sub.add_parser("evaluate")
    p.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.action == "prepare":
        prepare(args.base_analysis.resolve(), args.additional_manifest.resolve(), args.output_dir.resolve())
    else:
        evaluate(args.output_dir.resolve())
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(f"ERROR: {error}")
