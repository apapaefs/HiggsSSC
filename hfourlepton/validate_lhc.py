#!/usr/bin/env python3
"""Compare 13 TeV CMS-like particle-level acceptances with CMS Table 8.

This is a deliberately approximate acceptance benchmark, not a reproduction
of CMS detector reconstruction, its likelihood, or ATLAS's matrix-element
quadruplet assignment. Insufficient statistics cannot pass this check.
"""
import argparse
import html
import json
import numpy as np
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hfourlepton.analyze_fourlepton import load_dataset, read_card, signed_yield_summary

REFERENCE = {"ggf": (.402, .001), "vbf": (.445, .002), "wh": (.329, .002),
             "zh": (.340, .003), "tth": (.315, .004)}
SOURCE = "https://arxiv.org/abs/2103.04956"


def compare(data, metadata, *, relative_tolerance=.10, min_effective_events=1000):
    rows = []
    for mode, (reference, ref_error) in REFERENCE.items():
        samples = [s for s in metadata if s["production_mode"] == mode]
        if len(samples) != (6 if mode in ("wh", "zh") else 3):
            raise ValueError(f"missing decay channels or production components for {mode}")
        names = {s["name"] for s in samples}
        selected = data["eligible"] & np.isin(data["sample"], list(names))
        stats = signed_yield_summary(data["weight"], data["group"], data["source"], selected=selected)
        # At L=1 fb^-1 the denominator is the inclusive sigma*BR in fb.
        denominator = 1000 * sum(s["target_cross_section_pb"] for s in samples)
        acceptance = stats["sum_weight"] / denominator
        error = stats["mc_error"] / denominator
        relative = acceptance / reference - 1
        supported = stats["effective_count"] >= min_effective_events and 0 < acceptance < 1
        status = "insufficient_statistics" if not supported else "pass" if abs(relative) <= relative_tolerance else "investigate"
        rows.append(dict(mode=mode, acceptance=acceptance, mc_error=error,
            reference=reference, reference_stat_error=ref_error, relative_difference=relative,
            selected_effective_events=stats["effective_count"], status=status))
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args(argv)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("validation output is nonempty; choose a new output directory")
    raw = json.loads(args.manifest.read_text())
    if raw.get("collider", {}).get("sqrt_s_tev") != 13. or raw.get("higgs", {}).get("mass_gev") != 125.:
        raise ValueError("the CMS acceptance benchmark requires sqrt(s)=13 TeV and mH=125 GeV")
    card = read_card(REPO / "hfourlepton/analysis_cards/lhc13_cuts.json")
    data, metadata, _ = load_dataset(args.manifest, card, 1.)
    rows = compare(data, metadata)
    output = {"status": "pass_approximate_acceptance_check" if all(r["status"] == "pass" for r in rows) else "not_validated",
        "reference": SOURCE, "reference_table": 8, "rows": rows,
        "relative_tolerance": .10, "minimum_selected_effective_events": 1000,
        "tolerance_meaning": "10% relative acceptance agreement is a broad implementation benchmark, not the CMS theory uncertainty",
        "limitations": ["HwSim objects do not carry prompt/hadron ancestry; dressing and isolation use stable particles",
                        "No CMS lepton ID, impact parameter, vertex or detector-resolution model",
                        "Associated-production combinatorial signal is not separately fitted",
                        "No ATLAS matrix-element pairing or lepton-jet overlap emulation",
                        "Native NLO rates for non-ggF modes; pure-HEFT ihixs normalization for ggF"],
        "detector_validation": "not_tested_by_this_particle_level_comparison"}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "comparison.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    table = "".join(f"<tr><td>{r['mode']}</td><td>{r['acceptance']:.4f} ± {r['mc_error']:.4f}</td>"
                    f"<td>{r['reference']:.3f}</td><td>{r['relative_difference']:+.1%}</td><td>{r['status']}</td></tr>" for r in rows)
    (args.output_dir / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>LHC validation</title>'
        '<h1>CMS-like acceptance benchmark</h1><p>' + output["status"] + '</p><table><tr><th>Mode</th><th>Simulation</th>'
        '<th>CMS Table 8</th><th>Difference</th><th>Status</th></tr>' + table + '</table><p><a href="' + SOURCE + '">CMS reference</a></p><ul>' +
        ''.join('<li>' + html.escape(x) + '</li>' for x in output["limitations"]) + '</ul>')
    print(args.output_dir / "comparison.json")
    return 0 if output["status"] == "pass_approximate_acceptance_check" else 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit(str(error))
