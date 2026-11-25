"""Check electron-fake multiplicities and the resolved photon–quark generation boundary."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path


def fake_pairs(receipts):
    multiplicities = dict(all_masses=Counter(), higgs_window=Counter())
    attempted = saved = 0
    norm = 0.
    sources = {}
    for receipt in sorted(receipts):
        r = json.loads(receipt.read_text())
        n = r["attempted_events"]
        attempted += n
        saved += r["saved_events"]
        norm += n * 1000 * r["cross_section_pb"]
        sources[r["integration_id"]] = dict(cross_section_pb=r["cross_section_pb"],
                                            integration_error_pb=r["cross_section_error_pb"],
                                            hard_seed=r["hard_seed"], source_lhe_sha256=r["source_lhe_sha256"])
        with receipt.with_name("packet.csv").open() as stream:
            for _, event in itertools.groupby(csv.DictReader(stream), key=lambda row: row["source_event"]):
                masses = [float(row["mass_gev"]) for row in event if row["component"] == "electron_electron"]
                multiplicities["all_masses"][len(masses)] += 1
                multiplicities["higgs_window"][sum(123.3 <= m < 126.7 for m in masses)] += 1
    if not attempted:
        return {}
    out = dict(attempted=attempted, saved=saved, eligible_pair_multiplicities=multiplicities,
               hard_integrations=sources)
    hits = sum(nevents for pairs, nevents in multiplicities["higgs_window"].items() if pairs)
    if not hits:
        # Exact one-sided binomial limit on the probability of *any* eligible collision.
        # Display the dependence on maximum eligible unordered pairs, since the
        # experiment has no event-level upper limit of exactly two candidate electrons.
        upper = (norm / attempted) * (-math.expm1(math.log(.05) / attempted))
        p = .9 / 33.
        out["zero_count_95_upper_fb_by_maximum_pairs"] = {
            str(k): upper * .988 * min(1., k * p**2) for k in (1, 3, 6)
        }
        out["zero_count_95_upper_fb_arbitrary_multiplicity"] = upper * .988
    from scipy.stats import beta
    upper_p = 1. if hits == attempted else float(beta.ppf(.95, hits+1, attempted-hits))
    # For a specified pair, the event probability of two successful IDs is at
    # most p_e^2 regardless of how many other electron candidates are present.
    out["specified_electron_pair_95_upper_fb"] = norm/attempted * .988 * (.9/33.)**2 * upper_p
    out["window_events_with_any_electron_pair_and_distinct_tag"] = hits
    return out


def photon_top(particles, index):
    pending = [index]
    visited, tops = set(), set()
    while pending:
        i = pending.pop()
        if i in visited:
            continue
        visited.add(i)
        row = particles[i]
        code = int(row[0])
        if abs(code) == 6:
            tops.add(code)
            continue
        first, last = int(row[2]), int(row[3])
        if first:
            if not 1 <= first <= last <= i:
                raise ValueError("unexpected LHE particle-parent ordering")
            pending.extend(range(first-1,last))
    if len(tops) > 1:
        raise ValueError("photon associated with both top poles")
    return next(iter(tops),0)


def separations(lhe, origins=None, identities=None, muon_distances=None):
    out = []
    with gzip.open(lhe, "rt") as f:
        for line in f:
            if line.strip() != "<event>":
                continue
            n = int(next(f).split()[0])
            photons, quarks, muons = [], [], []
            records = [next(f).split() for _ in range(n)]
            if identities is not None:
                signature = "\n".join(" ".join([row[0],row[1],*row[6:10]]) for row in records)
                identity = hashlib.sha256(signature.encode()).digest()
                if identity in identities:
                    raise ValueError(f"repeated hard collision across source batches of {lhe}")
                identities.add(identity)
            signs = []
            for i, row in enumerate(records):
                particle, status = abs(int(row[0])), int(row[1])
                if status != 1 or particle not in (1, 2, 3, 4, 5, 13, 22):
                    continue
                px, py, pz = map(float, row[6:9])
                pt = math.hypot(px, py)
                if pt == 0.:
                    continue
                position = (math.asinh(pz/pt), math.atan2(py, px))
                if particle == 22:
                    photons.append(position)
                elif particle == 13:
                    muons.append(position)
                else:
                    quarks.append(position)
                if particle == 22 and origins is not None:
                    signs.append(photon_top(records, i))
            if not photons:
                raise ValueError("missing LHE photon")
            if origins is not None:
                origins[",".join(map(str,sorted(signs)))] += 1
            out.append(min((math.hypot(eta - qeta, math.remainder(phi - qphi, 2*math.pi))
                            for eta, phi in photons for qeta, qphi in quarks), default=math.inf))
            if muon_distances is not None:
                muon_distances.append(min((math.hypot(eta-meta, math.remainder(phi-mphi, 2*math.pi))
                                           for eta,phi in photons for meta,mphi in muons), default=math.inf))
    return out


def collinear(root, dileptonic, selected=None):
    output = {}
    seen = defaultdict(set)
    for hard_file in sorted(root.glob("tt_*dec*/round-*/hard.json")):
        meta = json.loads(hard_file.read_text())
        work = hard_file.parent
        sample = meta["sample"]
        sources = root / "packets" / sample / work.name
        receipts = [p for p in sorted(sources.glob("*/receipt.json"))
                    if selected is None or p.resolve() in selected]
        if not receipts:
            continue
        origins = Counter()
        muon_distances = []
        distances = separations(work / "events.lhe.gz", origins, seen[sample], muon_distances)
        expected = {"tt_0a_1dec":{"-6","6"}, "tt_0a_1dec_ll":{"-6","6"},
                    "tt_1a_1dec":{"-6,0","0,6"}, "tt_0a_2dec_same":{"-6,-6","6,6"},
                    "tt_0a_2dec_split":{"-6,6"}}
        if set(origins) - expected[sample] or sum(origins.values()) != meta["generated_events"]:
            raise ValueError(f"hard photon top ancestry or event count disagrees with {work}")
        bins = ((0., .05), (.05, .10), (.10, .20), (.20, .45), (.45, math.inf))
        acc = defaultdict(Counter)
        squared = defaultdict(Counter)
        muon_window, muon2_window = Counter(), Counter()
        attempted = 0
        for receipt in receipts:
            r = json.loads(receipt.read_text())
            if r["integration_id"] != meta["integration_id"]:
                raise ValueError("source mismatch")
            if (r["signed_attempted"] != r["attempted_events"] or
                    r["sumw2_attempted"] != r["attempted_events"]):
                raise ValueError("This diagnostic expects unit positive hard events")
            attempted += r["attempted_events"]
            with receipt.with_name("hard-photons.csv").open() as f:
                source = {int(row["saved_entry"]): r["source_first_event"] + int(row["lhe_entry"])
                          for row in csv.DictReader(f)}
            with receipt.with_name("packet.csv").open() as f:
                for event, outcomes in itertools.groupby(csv.DictReader(f), key=lambda row: int(row["source_event"])):
                    entry = source[event]
                    d = distances[entry]
                    b = next(i for i, (low, high) in enumerate(bins) if low <= d < high or d == high == math.inf)
                    weights = Counter()
                    for row in outcomes:
                        comp = row["component"]
                        if dileptonic and ((sample == "tt_0a_1dec" and comp == "decay_electron_ll") or
                                           (sample == "tt_0a_1dec_ll" and comp != "decay_electron_ll")):
                            continue
                        m = float(row["mass_gev"])
                        if 80. <= m < 140.:
                            p = float(row["probability"]) * float(row["evweight"])
                            weights["scan"] += p
                            if 123.3 <= m < 126.7:
                                weights["window"] += p
                    for region, value in weights.items():
                        acc[region][b] += value
                        squared[region][b] += value**2
                    muon_bin = "no hard muon" if math.isinf(muon_distances[entry]) else (
                        "dR < 0.05" if muon_distances[entry] < .05 else "dR >= 0.05")
                    muon_window[muon_bin] += weights["window"]
                    muon2_window[muon_bin] += weights["window"]**2
        if attempted > 1:
            scale = 1000*meta["cross_section_pb"]/attempted
            output[f"{sample}/{work.name}"] = dict(
                integration_id=meta["integration_id"], source_lhe_sha256=meta["lhe_sha256"],
                events=attempted, dr_edges=[0., .05, .10, .20, .45, None],
                hard_photon_top_ancestors=dict(origins),
                unique_hard_collisions_checked=len(distances),
                window_by_nearest_hard_muon={key: dict(rate_fb=value*scale,
                    mc_error_fb=scale*math.sqrt(max(0.,attempted*muon2_window[key]-value**2)/(attempted-1)))
                    for key,value in muon_window.items()},
                accepted_rates_fb={region: [acc[region][b]*scale for b in range(len(bins))]
                                   for region in ("scan", "window")},
                window_mc_error_fb=[
                    scale*math.sqrt(max(0.,attempted*squared["window"][b]-acc["window"][b]**2)/(attempted-1))
                    for b in range(len(bins))])
    return output


def frozen_inputs(top, spectrum):
    opener = gzip.open if spectrum.suffix == ".gz" else open
    with opener(spectrum, "rt") as stream:
        model = json.load(stream)
    receipts = set()
    counts = Counter()
    for entry in model["shards"]:
        if not entry["sample"].startswith("tt_"):
            continue
        path = (top.parent / entry["receipt"]).resolve()
        if not path.is_relative_to((top/"packets").resolve()) or path in receipts:
            raise ValueError("unexpected or repeated top receipt in the frozen spectrum")
        if (hashlib.sha256(path.read_bytes()).hexdigest() != entry["receipt_sha256"] or
            hashlib.sha256(path.with_name("packet.csv").read_bytes()).hexdigest() != entry["packet_sha256"]):
            raise ValueError(f"frozen source checksum changed in {path}")
        receipts.add(path)
        counts[entry["sample"]] += json.loads(path.read_text())["attempted_events"]
    sources = {k:v for k,v in model["sources"].items() if k.startswith("tt_")}
    if counts != {k:v["attempted_events"] for k,v in sources.items()} or not sources:
        raise ValueError("frozen source exposure disagrees with detector receipts")
    summary = dict(spectrum_sha256=hashlib.sha256(spectrum.read_bytes()).hexdigest(),
                   selected_shards=len(receipts),
                   sources={k:dict(attempted_events=v["attempted_events"], **v["window_125"])
                            for k,v in sources.items()},
                   top_window_rate_fb=sum(v["window_125"]["rate_fb"] for v in sources.values()),
                   top_window_error_fb=math.sqrt(sum(v["window_125"]["error_fb"]**2 for v in sources.values())))
    return receipts, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("top", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--ee-control", type=Path, help="independent t tbar -> direct ee hard/control campaign")
    parser.add_argument("--quark-collar", type=Path, help="independent lower-cut one-decay control campaign")
    parser.add_argument("--frozen-spectrum", type=Path, help="only use nominal receipts verified against this spectrum")
    parser.add_argument("--inclusive-one-decay", action="store_true",
                        help="Use all one-decay electron/jet outcomes when checking a standalone generation-cut control")
    args = parser.parse_args()
    args.top = args.top.resolve()
    selected, frozen = frozen_inputs(args.top, args.frozen_spectrum) if args.frozen_spectrum else (None, None)
    inclusive = (args.top/"packets/tt_0a").glob("*/*/receipt.json")
    payload = dict(tt_dielectron_fake=fake_pairs(p for p in inclusive if selected is None or p.resolve() in selected),
                   dilepton_enriched_combination=not args.inclusive_one_decay,
                   radiative_photon_quark=collinear(args.top, dileptonic=not args.inclusive_one_decay, selected=selected))
    if frozen:
        payload["frozen_top_window"] = frozen
    if args.ee_control:
        payload["direct_electron_pair_control"] = fake_pairs((args.ee_control/"packets").glob("*/receipt.json"))
    if args.quark_collar:
        # The disjoint angular interval uses the inclusive one-photon hard process
        # and samples all its eligible jet and electron fakes.
        rates = collinear(args.quark_collar, dileptonic=False)
        for name, result in rates.items():
            work = args.quark_collar / name
            check = json.loads((work/"collar-check.json").read_text())
            meta = json.loads((work/"hard.json").read_text())
            if result["events"] != check["events"] or not 0.01999 <= check["minimum"] <= check["maximum"] < 0.05001:
                raise ValueError("incomplete shower coverage or wrong hard photon-quark interval")
            result.update(hard_geometry_check=check, cross_section_pb=meta["cross_section_pb"],
                          cross_section_error_pb=meta["cross_section_error_pb"], hard_seed=meta["hard_seed"])
        if not rates:
            raise ValueError("quark–photon control has no completed showers")
        payload["single_radiative_quark_collar"] = rates
    tmp = args.output.with_name(args.output.name + ".partial")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(args.output)
