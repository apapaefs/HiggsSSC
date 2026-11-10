#!/usr/bin/env python3
"""Replace the three qq samples with an independent higher-statistics set.

One cross section per channel is used. Inputs and completed products are never
edited. Supports 13 TeV CMS v3 and 40 TeV SSC analyses. The assembled manifest,
strict report, cuts and ML live in a fresh directory with resumable,
checksum-verified stage receipts. Only CMS runs include an LHC acceptance check.
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
                  "normalization_rule", "matrix_element_generation_cuts", "selection_version",
                  "ancestry_schema_version", "photon_dressing_dr", "alternative_pair_veto",
                  "prompt_decay_selection", "response_mode")


def profile_settings(manifest):
    energy = manifest.get("collider", {}).get("sqrt_s_tev")
    profile = manifest.get("detector_profile")
    selection = manifest.get("selection_profile")
    if manifest.get("order") == "ho":
        if energy == 13. and selection == "cms" and profile == "perfect":
            return "lhc13", "CMS v3", "cms_fiducial_v3", 1
        if energy == 40. and selection == "ssc" and profile in ("perfect", "ssc"):
            return "ssc40", "SSC 40 TeV", "ssc_v1", 0
    raise ValueError("assembly supports 13 TeV HO CMS/perfect or 40 TeV HO SSC profiles")


def assemble(base, replacement, provenance):
    _, label, selection_version, ancestry_version = profile_settings(base)
    old = {s["name"]: s for s in base["samples"]}
    new = {s["name"]: s for s in replacement["samples"]}
    if len(old) != len(base["samples"]) or len(new) != len(replacement["samples"]):
        raise ValueError("duplicate sample names")
    if set(new) != QQ or not QQ <= old.keys():
        raise ValueError("replacement must contain exactly the three qq continuum channels")
    for field in ("collider", "selection_profile", "detector_profile", "detector", "order",
                  "generation_cuts", "cut_stages"):
        if base.get(field) != replacement.get(field):
            raise ValueError("incompatible campaigns: " + field)
    if (base.get("sampling", {}).get("continuum_m4l_gev")
            != replacement.get("sampling", {}).get("continuum_m4l_gev")):
        raise ValueError("incompatible campaigns: continuum generation window")
    for s in [*old.values(), *new.values()]:
        if s.get("selection_version") != selection_version or s.get("ancestry_schema_version") != ancestry_version:
            raise ValueError(f"all samples must use {label} selection {selection_version}")
    for name in QQ:
        for field in PHYSICS_FIELDS:
            if old[name].get(field) != new[name].get(field):
                raise ValueError(f"replacement changes {name} {field}")
        if old[name].get("input_lhe_sha256") == new[name].get("input_lhe_sha256"):
            raise ValueError("replacement is the same hard-event population: " + name)
        old_seed = old[name].get("seeds", {}).get("production")
        if old_seed is not None and old_seed == new[name].get("seeds", {}).get("production"):
            raise ValueError("replacement reuses the production seed: " + name)
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
        "description": f"{label}; independent higher-statistics qq samples replace the original three channels"}
    result["software"] = {"source_campaigns": {
        "base": base.get("software", {}), "replacement": replacement.get("software", {})}}
    result["tracked_inputs"] = {"assembly": revision.artifact(Path(__file__).resolve())}
    result["status"] = "complete"
    return result


def audit_profile_pair(ssc, perfect):
    """The perfect/SSC comparison must use identical generated and showered events."""
    if profile_settings(perfect)[0] != "ssc40" or perfect["detector_profile"] != "perfect":
        raise ValueError("paired manifest must be a 40 TeV perfect-response SSC analysis")
    left = {s["name"]: s for s in ssc["samples"]}
    right = {s["name"]: s for s in perfect["samples"]}
    if left.keys() != right.keys() or len(right) != len(perfect["samples"]):
        raise ValueError("perfect/SSC sample inventories differ")
    for name in left:
        for field in ("input_lhe_sha256", "source_lhe_sha256", "seeds", "nevents_requested"):
            if not left[name].get(field) or left[name][field] != right[name].get(field):
                raise ValueError(f"perfect/SSC {field} differs for {name}")
        def roots(sample):
            return (sample.get("shower_completion") or {}).get("artifacts", {}).get("roots")
        if not roots(left[name]) or roots(left[name]) != roots(right[name]):
            raise ValueError(f"perfect/SSC showers differ for {name}")


def statistical_support(cuts, ml):
    cut_totals = cuts["cuts_full_sample"]["totals"]
    test = ml.get("xgboost_test")
    ml_totals = test["totals"] if test else None
    return {"cuts_supported": bool(cut_totals["statistically_supported"]),
            "cuts_background_effective_events": cut_totals["background"]["effective_count"],
            "ml_supported": bool(ml.get("threshold") is not None and ml_totals
                                 and ml_totals["statistically_supported"]),
            "ml_test_background_effective_events": ml_totals["background"]["effective_count"] if ml_totals else None}


def execute(base_path, replacement_path, output, luminosity, *, paired_perfect=False):
    base, _ = revision.audit_completed_manifest(base_path)
    replacement, _ = revision.audit_completed_manifest(replacement_path)
    card_prefix, label, _, _ = profile_settings(base)
    if card_prefix == "ssc40" and base["detector_profile"] != "ssc":
        raise ValueError("use manifest_ssc.json as the primary SSC input")
    if paired_perfect and card_prefix != "ssc40":
        raise ValueError("--paired-perfect applies only to SSC assembly")
    pairs = None
    input_paths = {"base": base_path, "replacement": replacement_path}
    if paired_perfect:
        input_paths.update(base_perfect=base_path.with_name("manifest_perfect.json"),
                           replacement_perfect=replacement_path.with_name("manifest_perfect.json"))
        pairs = [revision.audit_completed_manifest(input_paths[key])[0]
                 for key in ("base_perfect", "replacement_perfect")]
        audit_profile_pair(base, pairs[0])
        audit_profile_pair(replacement, pairs[1])
    if any(path.is_relative_to(output) for path in input_paths.values()):
        raise ValueError("assembly directory must not contain an input campaign")
    sources = [Path(__file__).resolve(), REPO / "hfourlepton/reprocess.py",
               REPO / "hfourlepton/make_fourlepton_report.py", REPO / "hfourlepton/validate_lhc.py",
               REPO / "hfourlepton/report_reader.py", REPO / "hfourlepton/observables.py", REPO / "hfourlepton/analyze_fourlepton.py",
               REPO / "ho_xgboost_core.py",
               *[REPO / f"hfourlepton/analysis_cards/{card_prefix}_{k}.json" for k in ("cuts", "xgboost")]]
    identity = {**{key: revision.artifact(path) for key, path in input_paths.items()},
                "luminosity_fb": luminosity, "source_hashes": {str(p.relative_to(REPO)): revision.digest(p) for p in sources}}
    combined = assemble(base, replacement, identity)
    combined_perfect = assemble(*pairs, {**identity, "base": identity["base_perfect"],
                                         "replacement": identity["replacement_perfect"]}) if pairs else None
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
        manifest = output / f"manifest_{base['detector_profile']}.json"
        if state.get("manifest_sha256"):
            revision.verify(manifest, state["manifest_sha256"])
        else:
            revision.save(manifest, combined)
            state["manifest_sha256"] = revision.digest(manifest)
            revision.save(state_path, state)
        manifests = [manifest]
        if combined_perfect is not None:
            perfect_path = output / "manifest_perfect.json"
            if state.get("perfect_manifest_sha256"):
                revision.verify(perfect_path, state["perfect_manifest_sha256"])
            else:
                revision.save(perfect_path, combined_perfect)
                state["perfect_manifest_sha256"] = revision.digest(perfect_path)
                revision.save(state_path, state)
            manifests.insert(0, perfect_path)
        stage_id = revision.fingerprint({"identity": identity, "manifest": state["manifest_sha256"],
                                         "perfect_manifest": state.get("perfect_manifest_sha256")})
        py = sys.executable
        revision.stage(output, "report", stage_id, lambda work: [py, REPO / "hfourlepton/make_fourlepton_report.py",
                       *manifests, "--output-dir", work, "--luminosity-fb", luminosity, "--strict", "--observables"])
        validation = {"result": "not_applicable_ssc"}
        stage_names = ["report", "cuts", "xgboost"]
        if card_prefix == "lhc13":
            validation = revision.stage(output, "lhc_validation", stage_id, lambda work: [py, REPO / "hfourlepton/validate_lhc.py",
                               "--manifest", manifest, "--output-dir", work], validation=True)
            stage_names.insert(1, "lhc_validation")
        for name, minimum in (("cuts", 400.), ("xgboost", 100.)):
            card = json.loads((REPO / f"hfourlepton/analysis_cards/{card_prefix}_{name}.json").read_text())
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
                        for name in stage_names)
        validation_text = f'<p>{validation["result"]}</p>' if card_prefix == "lhc13" else ""
        (output / "index.html").write_text(f'<!doctype html><meta charset="utf-8"><title>{label} analysis</title>'
            f'<h1>{label} with independent qq statistics</h1>{validation_text}'
            f'<p>Cut and ML statistical support: {statistics}</p><ul>{links}</ul>')
        print("Completed:", output / "index.html")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base-manifest", type=Path, required=True)
    p.add_argument("--qq-manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--luminosity-fb", type=float, default=100.)
    p.add_argument("--paired-perfect", action="store_true",
                   help="SSC: include sibling manifest_perfect.json inputs in the comparison report")
    args = p.parse_args()
    if not 0 < args.luminosity_fb < float("inf"):
        p.error("luminosity must be positive and finite")
    execute(args.base_manifest.resolve(), args.qq_manifest.resolve(), args.output_dir.resolve(), args.luminosity_fb,
            paired_perfect=args.paired_perfect)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError, KeyError) as exc:
        sys.exit(f"Error: {exc}")
