"""Check electron-fake multiplicities and the resolved photon–quark generation boundary."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import itertools
import json
import math
from pathlib import Path


def fake_pairs(root):
    multiplicities = dict(all_masses=Counter(), higgs_window=Counter())
    attempted = saved = 0
    norm = 0.
    for receipt in sorted((root / "packets/tt_0a").glob("*/*/receipt.json")):
        r = json.loads(receipt.read_text())
        n = r["attempted_events"]
        attempted += n
        saved += r["saved_events"]
        norm += n * 1000 * r["cross_section_pb"]
        with receipt.with_name("packet.csv").open() as stream:
            for _, event in itertools.groupby(csv.DictReader(stream), key=lambda row: row["source_event"]):
                masses = [float(row["mass_gev"]) for row in event if row["component"] == "electron_electron"]
                multiplicities["all_masses"][len(masses)] += 1
                multiplicities["higgs_window"][sum(123.3 <= m < 126.7 for m in masses)] += 1
    if not attempted:
        return {}
    out = dict(attempted=attempted, saved=saved, eligible_pair_multiplicities=multiplicities)
    if not any(m for m in multiplicities["higgs_window"] if m):
        # Exact one-sided binomial limit on the probability of *any* eligible collision.
        # Display the dependence on maximum eligible unordered pairs, since the
        # experiment has no event-level upper limit of exactly two candidate electrons.
        upper = (norm / attempted) * (-math.expm1(math.log(.05) / attempted))
        p = .9 / 33.
        out["zero_count_95_upper_fb_by_maximum_pairs"] = {
            str(k): upper * .988 * min(1., k * p**2) for k in (1, 3, 6)
        }
        out["zero_count_95_upper_fb_arbitrary_multiplicity"] = upper * .988
    return out


def separations(lhe):
    out = []
    with gzip.open(lhe, "rt") as f:
        for line in f:
            if line.strip() != "<event>":
                continue
            n = int(next(f).split()[0])
            photons, quarks = [], []
            for _ in range(n):
                row = next(f).split()
                particle, status = abs(int(row[0])), int(row[1])
                if status != 1 or particle not in (1, 2, 3, 4, 5, 22):
                    continue
                px, py, pz = map(float, row[6:9])
                pt = math.hypot(px, py)
                if pt == 0.:
                    continue
                (photons if particle == 22 else quarks).append((math.asinh(pz/pt), math.atan2(py, px)))
            if not photons:
                raise ValueError("missing LHE photon")
            out.append(min((math.hypot(eta - qeta, math.remainder(phi - qphi, 2*math.pi))
                            for eta, phi in photons for qeta, qphi in quarks), default=math.inf))
    return out


def collinear(root):
    output = {}
    for hard_file in sorted(root.glob("tt_*dec*/round-*/hard.json")):
        meta = json.loads(hard_file.read_text())
        work = hard_file.parent
        sample = meta["sample"]
        distances = separations(work / "events.lhe.gz")
        bins = ((0., .05), (.05, .10), (.10, .20), (.20, .45), (.45, math.inf))
        acc = defaultdict(Counter)
        attempted = 0
        sources = root / "packets" / sample / work.name
        for receipt in sorted(sources.glob("*/receipt.json")):
            r = json.loads(receipt.read_text())
            if r["integration_id"] != meta["integration_id"]:
                raise ValueError("source mismatch")
            attempted += r["attempted_events"]
            with receipt.with_name("hard-photons.csv").open() as f:
                source = {int(row["saved_entry"]): r["source_first_event"] + int(row["lhe_entry"])
                          for row in csv.DictReader(f)}
            with receipt.with_name("packet.csv").open() as f:
                for row in csv.DictReader(f):
                    m = float(row["mass_gev"])
                    if not 80. <= m <= 140.:
                        continue
                    d = distances[source[int(row["source_event"])]]
                    b = next(i for i, (low, high) in enumerate(bins) if low <= d < high or d == high == math.inf)
                    p = float(row["probability"]) * float(row["evweight"])
                    acc["scan"][b] += p
                    if 123.3 <= m < 126.7:
                        acc["window"][b] += p
        if attempted:
            scale = 1000*meta["cross_section_pb"]/attempted
            output[meta["integration_id"]] = dict(
                events=attempted, dr_edges=[0., .05, .10, .20, .45, None],
                accepted_rates_fb={region: [counts[b]*scale for b in range(len(bins))]
                                   for region, counts in acc.items()})
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("top", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    payload = dict(tt_dielectron_fake=fake_pairs(args.top), radiative_photon_quark=collinear(args.top))
    tmp = args.output.with_name(args.output.name + ".partial")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(args.output)
