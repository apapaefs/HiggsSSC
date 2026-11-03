#!/usr/bin/env python3
"""Replace the three qq samples with an independently generated CMS v3 set.

One cross section per channel is used. Inputs and completed products are never
edited. The assembled manifest, strict report, acceptance check, cuts and ML
live in a fresh directory with resumable, checksum-verified stage receipts.
"""
import argparse
import copy
import fcntl
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hfourlepton import reprocess as revision

QQ = frozenset({"bkg_qq_4e", "bkg_qq_4mu", "bkg_qq_2e2mu"})
PHYSICS_FIELDS = ("category", "channel", "model", "process", "matching", "perturbative_order",
                  "production_mode", "hard_model", "pdf_name", "pdf_lhaid", "higgs_mass_gev",
                  "normalization_rule", "matrix_element_generation_cuts")


def assemble(base, replacement, provenance):
    old = {s["name"]: s for s in base["samples"]}
    new = {s["name"]: s for s in replacement["samples"]}
    if len(old) != len(base["samples"]) or len(new) != len(replacement["samples"]):
        raise ValueError("duplicate sample names")
    if set(new) != QQ or not QQ <= old.keys():
        raise ValueError("replacement must contain exactly the three qq continuum channels")
    for field in ("collider", "selection_profile", "detector_profile", "order", "generation_cuts"):
        if base.get(field) != replacement.get(field):
            raise ValueError("incompatible campaigns: " + field)
    for s in [*old.values(), *new.values()]:
        if s.get("selection_version") != "cms_fiducial_v3" or s.get("ancestry_schema_version") != 1:
            raise ValueError("all samples must use ancestry-aware CMS v3")
    for name in QQ:
        for field in PHYSICS_FIELDS:
            if old[name].get(field) != new[name].get(field):
                raise ValueError(f"replacement changes {name} {field}")
        if old[name].get("input_lhe_sha256") == new[name].get("input_lhe_sha256"):
            raise ValueError("replacement is the same hard-event population: " + name)
        if new[name]["nevents_requested"] <= old[name]["nevents_requested"]:
            raise ValueError("replacement does not increase the event count: " + name)
    result = copy.deepcopy(base)
    result["samples"] = []
    for s in base["samples"]:
        entry = copy.deepcopy(new.get(s["name"], s))
        entry["source_manifest"] = provenance["replacement" if s["name"] in new else "base"]
        result["samples"].append(entry)
    result["analysis_revision"] = {"kind": "replace_qq_backgrounds", **provenance}
    result["configuration_fingerprint"] = revision.fingerprint(result["analysis_revision"])
    result["sampling"] = {"preset": "custom", "events_by_sample": {
        s["name"]: s["nevents_requested"] for s in result["samples"]},
        "description": "CMS v3; independent higher-statistics qq samples replace the original three channels"}
    result["software"] = {"source_campaigns": {
        "base": base.get("software", {}), "replacement": replacement.get("software", {})}}
    result["tracked_inputs"] = {"assembly": revision.artifact(Path(__file__).resolve())}
    result["status"] = "complete"
    return result


def statistical_support(cuts, ml):
    cut_totals = cuts["cuts_full_sample"]["totals"]
    test = ml.get("xgboost_test")
    ml_totals = test["totals"] if test else None
    return {"cuts_supported": bool(cut_totals["statistically_supported"]),
            "cuts_background_effective_events": cut_totals["background"]["effective_count"],
            "ml_supported": bool(ml.get("threshold") is not None and ml_totals
                                 and ml_totals["statistically_supported"]),
            "ml_test_background_effective_events": ml_totals["background"]["effective_count"] if ml_totals else None}


def execute(base_path, replacement_path, output, luminosity):
    base, _ = revision.audit_parent(base_path)
    replacement, _ = revision.audit_parent(replacement_path)
    if base_path.is_relative_to(output) or replacement_path.is_relative_to(output):
        raise ValueError("assembly directory must not contain an input campaign")
    sources = [Path(__file__).resolve(), REPO / "hfourlepton/reprocess.py",
               REPO / "hfourlepton/make_fourlepton_report.py", REPO / "hfourlepton/validate_lhc.py",
               REPO / "hfourlepton/analyze_fourlepton.py", REPO / "ho_xgboost_core.py",
               *[REPO / f"hfourlepton/analysis_cards/lhc13_{k}.json" for k in ("cuts", "xgboost")]]
    identity = {"base": revision.artifact(base_path), "replacement": revision.artifact(replacement_path),
                "luminosity_fb": luminosity, "source_hashes": {str(p.relative_to(REPO)): revision.digest(p) for p in sources}}
    combined = assemble(base, replacement, identity)
    if output.exists() and any(output.iterdir()) and not (output / "state.json").is_file():
        raise ValueError("nonempty output without assembly state")
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state_path = output / "state.json"
        state = json.loads(state_path.read_text()) if state_path.exists() else {"identity": identity, "state": "running"}
        if state["identity"] != identity:
            raise ValueError("assembly inputs or software changed; use a new output directory")
        revision.save(state_path, state)
        manifest = output / "manifest_perfect.json"
        if state.get("manifest_sha256"):
            revision.verify(manifest, state["manifest_sha256"])
        else:
            revision.save(manifest, combined)
            state["manifest_sha256"] = revision.digest(manifest)
            revision.save(state_path, state)
        stage_id = revision.fingerprint({"identity": identity, "manifest": state["manifest_sha256"]})
        py = sys.executable
        revision.stage(output, "report", stage_id, lambda work: [py, REPO / "hfourlepton/make_fourlepton_report.py",
                       manifest, "--output-dir", work, "--luminosity-fb", luminosity, "--strict"])
        validation = revision.stage(output, "lhc_validation", stage_id, lambda work: [py, REPO / "hfourlepton/validate_lhc.py",
                           "--manifest", manifest, "--output-dir", work], validation=True)
        for name, minimum in (("cuts", 400.), ("xgboost", 100.)):
            card = json.loads((REPO / f"hfourlepton/analysis_cards/lhc13_{name}.json").read_text())
            card["min_background_effective_events"] = minimum
            path = output / (name + "_card.json")
            if path.exists() and json.loads(path.read_text()) != card:
                raise ValueError("assembly analysis card changed")
            if not path.exists():
                revision.save(path, card)
            revision.stage(output, name, stage_id, lambda work: [py, REPO / "hfourlepton/analyze_fourlepton.py", name,
                   "--manifest", manifest, "--config", path, "--output-dir", work, "--luminosity-fb", luminosity])
        support = statistical_support(json.loads((output / "cuts/summary.json").read_text()),
                                      json.loads((output / "xgboost/summary.json").read_text()))
        statistics = "supported" if support["cuts_supported"] and support["ml_supported"] else "more_statistics_needed"
        state.update(state="complete", validation_status=validation["result"],
                     statistical_support=support, statistics_status=statistics)
        revision.save(state_path, state)
        links = "".join(f'<li><a href="{name}/index.html">{name}</a></li>'
                        for name in ("report", "lhc_validation", "cuts", "xgboost"))
        (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>CMS v3 validation</title>'
            f'<h1>CMS v3 with independent qq statistics</h1><p>{validation["result"]}</p>'
            f'<p>Cut and ML statistical support: {statistics}</p><ul>{links}</ul>')
        print("Completed:", output / "index.html")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-manifest", type=Path, required=True)
    p.add_argument("--qq-manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--luminosity-fb", type=float, default=100.)
    args = p.parse_args()
    if not 0 < args.luminosity_fb < float("inf"):
        p.error("luminosity must be positive and finite")
    execute(args.base_manifest.resolve(), args.qq_manifest.resolve(), args.output_dir.resolve(), args.luminosity_fb)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        sys.exit(f"Error: {exc}")
