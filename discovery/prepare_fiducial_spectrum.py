"""Export the inclusive SSC four-lepton mass spectrum before the classifier cut."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def project(report):
    summary = json.loads((report / "report.json").read_text())
    samples = {r["name"]: r for r in summary["samples"]
               if r["detector_profile"] == "ssc" and r["included_in_physics_totals"]}
    edges = np.arange(70, 202, 2)
    groups = {"signal_fb": "signal_", "qq_continuum_fb": "bkg_qq_", "gg_continuum_fb": "bkg_gg_",
              "zbb_ee_fb": "bkg_zbb_ee", "zbb_mumu_fb": "bkg_zbb_mumu", "ttbar_fb": "bkg_ttbar"}
    sample_group = {}
    for sample in samples:
        key, = [name for name, prefix in groups.items() if sample.startswith(prefix)]
        sample_group[sample] = key
    expected = {"signal_fb": 21, "qq_continuum_fb": 3, "gg_continuum_fb": 3,
                "zbb_ee_fb": 1, "zbb_mumu_fb": 1, "ttbar_fb": 1}
    if any(sum(group == key for group in sample_group.values()) != count for key, count in expected.items()):
        raise ValueError("Incomplete set of signal, continuum or reducible samples")
    columns = {key: np.zeros(len(edges) - 1) for key in groups}
    seen = {sample: set() for sample in samples}
    totals = dict.fromkeys(samples, 0.)
    with (report / "observables/data/m4l.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            if (row["profile"], row["region"], row["channel"]) != ("ssc", "selected", "inclusive"):
                continue
            sample = row["sample"]
            rate = float(row["cross_section_pb"])
            totals[sample] += rate
            if row["bin_type"] != "regular":
                continue
            low, high = float(row["low"]), float(row["high"])
            index = np.searchsorted(edges, low)
            if not (index < len(edges) - 1 and (low, high) == tuple(edges[index:index + 2])
                    and index not in seen[sample]):
                raise ValueError(f"Invalid or duplicate fiducial bin for {sample} at {low}")
            seen[sample].add(index)
            columns[sample_group[sample]][index] += 1000 * rate
    for name, row in samples.items():
        if (len(seen[name]) != len(edges) - 1
                or not math.isclose(totals[name], row["selected"]["cross_section_pb"], rel_tol=1e-7, abs_tol=1e-12)):
            raise ValueError(f"Fiducial histogram fails inclusive cut-flow check for {name}")
    return dict(low_gev=edges[:-1], high_gev=edges[1:], **columns)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path, help="Four-lepton student report produced with --observables")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    columns = project(args.report)
    with args.output.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(zip(*columns.values()))
