"""Combine audited, disjoint tagged detector packets at the independent hard-event level."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SAMPLES = {"w_lnu_aa", "z_ll_aa", "z_ee_a", "tt_0a", "tt_1a", "tt_2a"}
RADIATIVE = {"tt_0a_1dec", "tt_1a_1dec", "tt_0a_2dec_same", "tt_0a_2dec_split"}
COMPONENTS = {
    "tt_0a":{"decay_decay","decay_electron","decay_jet","electron_electron","electron_jet","jet_jet"},
    "tt_1a":{"decay_hard","electron_hard","hard_jet"},
    "tt_2a":{"hard_hard"}, "tt_0a_1dec":{"decay_electron","decay_jet"},
    "tt_0a_1dec_ll":{"decay_electron_ll"},
    "tt_1a_1dec":{"decay_hard"}, "tt_0a_2dec_same":{"decay_decay"},
    "tt_0a_2dec_split":{"decay_decay"}, "w_lnu_aa":{"hard_hard"},
    "z_ll_aa":{"hard_hard"}, "z_ee_a":{"electron_hard"},
}
EDGES = np.sort(np.r_[np.arange(75., 146.), 123.3, 126.7])
SIZE = len(EDGES) - 1
WINDOW = (EDGES[:-1] >= 123.3) & (EDGES[1:] <= 126.7)


class Moments:
    def __init__(self):
        self.events = self.saved = self.shards = self.window_count = 0
        self.sumw = self.sumw2 = 0.
        self.s1 = np.zeros(SIZE)
        self.gluon = np.zeros(SIZE)
        self.weighted_s1 = np.zeros(SIZE)
        self.s2 = np.zeros((SIZE, SIZE))
        self.counts = Counter()
        self.component = Counter()
        self.component_sq = Counter()
        self.component_w = Counter()
        self.max_component = Counter()
        self.max_window = 0.
        self.components_seen = Counter()

    def add_packet(self, path, receipt, replace_shower_decays, dileptonic):
        saved = 0
        weights = weights2 = 0.
        with path.open() as packet:
            reader = csv.DictReader(packet)
            sample = receipt["sample"]
            split_dilepton = dileptonic and sample in {"tt_0a_1dec","tt_0a_1dec_ll"}
            if split_dilepton and "direct_top_dilepton" not in reader.fieldnames:
                raise ValueError(f"Missing direct-dilepton classification in {path}")
            grouped = itertools.groupby(reader, lambda row: int(row["source_event"]))
            for event, outcomes in grouped:
                if event != saved:
                    raise ValueError(f"Missing or duplicated shower event {saved} in {path}")
                saved += 1
                central, gluon = Counter(), Counter()
                components = Counter()
                total = total_gluon = 0.
                weight = None
                for row in outcomes:
                    w, mass, p, g = (float(row[k]) for k in ("evweight", "mass_gev", "probability",
                                                             "probability_gluon_jets"))
                    if (not all(map(math.isfinite, (w, mass, p, g))) or not 0. <= p <= 1.
                            or not 0. <= g <= 1. or (weight is not None and w != weight)):
                        raise ValueError(f"Invalid detector outcome in {path} event {event}")
                    weight = w
                    total += p
                    total_gluon += g
                    component = row["component"]
                    if replace_shower_decays and "decay" in component.split("_"):
                        continue
                    if split_dilepton and ((sample=="tt_0a_1dec" and component=="decay_electron_ll")
                                          or (sample=="tt_0a_1dec_ll" and component!="decay_electron_ll")):
                        continue
                    if not split_dilepton and component=="decay_electron_ll":
                        component="decay_electron"
                    if p == 0. or not EDGES[0] <= mass < EDGES[-1]:
                        continue
                    i = int(np.searchsorted(EDGES, mass, side="right")) - 1
                    central[i] += w * p
                    gluon[i] += w * g
                    self.components_seen[component] += 1
                    if WINDOW[i]:
                        components[component] += w * p
                if max(total, total_gluon) > 1. + 1.e-10:
                    raise ValueError(f"Non-exclusive outcomes in {path} event {event}")
                weights += weight
                weights2 += weight * weight
                for i, yi in central.items():
                    self.s1[i] += yi
                    self.gluon[i] += gluon[i]
                    self.weighted_s1[i] += weight * yi
                    for j, yj in central.items():
                        self.s2[i, j] += yi * yj
                if components:
                    self.window_count += 1
                    self.component.update(components)
                    self.counts.update(components.keys())
                    self.max_window = max(self.max_window, abs(sum(components.values())))
                    for name, value in components.items():
                        self.component_sq[name] += value*value
                        self.component_w[name] += value*weight
                        self.max_component[name] = max(self.max_component[name], abs(value))
        if (saved != receipt["saved_events"] or saved > receipt["attempted_events"]
                or not math.isclose(weights, receipt["signed_saved"], abs_tol=1.e-6, rel_tol=1.e-9)
                or not math.isclose(weights2, receipt["sumw2_saved"], abs_tol=1.e-6, rel_tol=1.e-9)
                or receipt["sumw2_attempted"] + 1.e-6 < weights2
                or not math.isfinite(receipt["signed_attempted"])
                or not math.isfinite(receipt["sumw2_attempted"])):
            raise ValueError(f"Attempted/saved shower weights disagree with {path}")
        self.sumw += receipt["signed_attempted"]
        self.sumw2 += receipt["sumw2_attempted"]
        self.saved += saved
        self.events += receipt["attempted_events"]
        self.shards += 1


def window_report(rate, mc, integration, count):
    r = float(np.sum(rate[WINDOW]))
    index = np.ix_(WINDOW, WINDOW)
    emc = math.sqrt(max(0., float(np.sum(mc[index]))))
    eint = math.sqrt(max(0., float(np.sum(integration[index]))))
    return dict(rate_fb=r, mc_error_fb=emc, integration_error_fb=eint, error_fb=math.hypot(emc, eint),
                contributing_events=count, effective_events=(r/emc)**2 if emc else None)


def make_report(rate, gluon, mc, integration, count):
    return dict(rates_fb=rate.tolist(), gluon_jet_rates_fb=gluon.tolist(),
                covariance_fb2=(mc+integration).tolist(), mc_covariance_fb2=mc.tolist(),
                integration_covariance_fb2=integration.tolist(), window_125=window_report(rate, mc, integration, count))


def read_campaign(root, expected, radiative, region, dileptonic):
    groups = {}
    provenance, hashes, seeds = [], set(), set()
    intervals = defaultdict(list)
    event_hashes = set()
    for path in sorted(root.rglob("receipt.json")):
        receipt = json.loads(path.read_text())
        if "sample" not in receipt or receipt["sample"] not in expected:
            continue
        if receipt.get("generation_region", "nominal") != region:
            continue
        sample, seed, sha = receipt["sample"], receipt["seed"], receipt["lhe_sha256"]
        ident = (sample, receipt["integration_id"])
        block = receipt["event_block_sha256"]
        start, stop = receipt["source_first_event"], receipt["source_stop_event"]
        if not 0 <= start < stop or stop-start != receipt["attempted_events"]:
            raise ValueError(f"Malformed LHE source-event interval in {path}")
        intervals[receipt["source_lhe_sha256"]].append((start, stop))
        if sha in hashes or block in event_hashes or (sample, seed) in seeds:
            raise ValueError(f"Repeated LHE shard or seed for {path}")
        hashes.add(sha)
        event_hashes.add(block)
        seeds.add((sample, seed))
        sigma = float(receipt["cross_section_pb"]) * 1000.
        error = float(receipt["cross_section_error_pb"]) * 1000.
        if not math.isfinite(sigma) or sigma <= 0. or not math.isfinite(error) or error < 0.:
            raise ValueError(f"Invalid hard-process cross section in {path}")
        if ident not in groups:
            groups[ident] = (sigma, error, Moments())
        xs, dxs, moments = groups[ident]
        if sigma != xs or error != dxs:
            raise ValueError(f"Cross section changed within integration {ident}")
        packet = path.with_name("packet.csv")
        packet_sha = hashlib.sha256(packet.read_bytes()).hexdigest()
        if receipt.get("packet_sha256", packet_sha) != packet_sha:
            raise ValueError(f"Packet contents changed in {path.parent}")
        moments.add_packet(packet, receipt, radiative and sample in {"tt_0a", "tt_1a"}, dileptonic)
        provenance.append(dict(receipt=str(path.relative_to(root)), packet_sha256=packet_sha,
                               receipt_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), **{
                                   k:receipt[k] for k in ("sample","seed","integration_id","lhe_sha256",
                                                         "source_lhe_sha256","source_first_event","source_stop_event",
                                                         "event_block_sha256","response_sha256")}))
    if not groups:
        raise ValueError(f"No completed tagged-background receipts under {root}")
    for parent, blocks in intervals.items():
        blocks.sort()
        if any(right[0] < left[1] for left, right in zip(blocks, blocks[1:])):
            raise ValueError(f"Overlapping LHE source-event ranges under {parent}")
    return groups, provenance


def summarize(root, radiative=False, region="nominal", dileptonic=False):
    radiative = radiative or dileptonic
    expected = SAMPLES | RADIATIVE if radiative else SAMPLES
    if dileptonic:
        expected = expected | {"tt_0a_1dec_ll"}
    if region == "outside_nominal":
        if radiative:
            raise ValueError("The complementary phase-space control applies to electroweak processes")
        expected = {"w_lnu_aa","z_ll_aa","z_ee_a"}
    groups, provenance = read_campaign(root, expected, radiative, region, dileptonic)
    sources = {}
    for sample in sorted({key[0] for key in groups}):
        batches = [(key[1], *values) for key, values in groups.items() if key[0] == sample]
        events = sum(m.events for _, _, _, m in batches)
        rate, gluon = np.zeros(SIZE), np.zeros(SIZE)
        mc, integration = np.zeros((SIZE, SIZE)), np.zeros((SIZE, SIZE))
        count = saved = shards = 0
        components, component_counts, in_range = Counter(), Counter(), Counter()
        component_variance, component_max = Counter(), Counter()
        max_packet = 0.
        details = []
        for group, sigma, error, m in batches:
            if m.events < 2 or m.sumw <= 0.:
                raise ValueError(f"Ill-defined normalization for {sample} {group}")
            # Equivalent integrations of one hard process are exposure-weighted
            # independent estimates. Sibling showers of the same integration
            # share its cross-section error and were pooled before this step.
            alpha = m.events / events
            mean = m.s1 / m.sumw
            rate += alpha * sigma * mean
            gluon += alpha * sigma * m.gluon / m.sumw
            cross = np.outer(m.weighted_s1, mean)
            residual = m.s2 - cross - cross.T + m.sumw2 * np.outer(mean, mean)
            mc += (alpha*sigma/m.sumw)**2 * (m.events/(m.events-1.)) * residual
            integration += (alpha*error)**2 * np.outer(mean, mean)
            count += m.window_count
            saved += m.saved
            shards += m.shards
            component_counts.update(m.counts)
            in_range.update(m.components_seen)
            for name, value in m.component.items():
                scale = alpha*sigma/m.sumw
                components[name] += scale * value
                avg = value/m.sumw
                component_variance[name] += scale*scale*m.events/(m.events-1.)*(
                    m.component_sq[name]-2*avg*m.component_w[name]+avg*avg*m.sumw2)
                component_max[name] = max(component_max[name],scale*m.max_component[name])
            max_packet = max(max_packet, alpha*sigma/m.sumw*m.max_window)
            details.append(dict(integration_id=group, generated=m.events, signed_attempted=m.sumw,
                                sumw2_attempted=m.sumw2, saved=m.saved, cross_section_fb=sigma,
                                cross_section_error_fb=error, shards=m.shards))
        expected_components = COMPONENTS[sample]
        if radiative and sample in {"tt_0a","tt_1a"}:
            expected_components = {k for k in expected_components if "decay" not in k}
        missing_components = sorted(expected_components-set(components))
        sources[sample] = dict(**make_report(rate, gluon, mc, integration, count), attempted_events=events,
                               saved_events=saved, shards=shards, integration_groups=details,
                               component_window_rates_fb=dict(components), component_window_events=dict(component_counts),
                               component_window_mc_error_fb={k:math.sqrt(max(0,v)) for k,v in component_variance.items()},
                               component_largest_window_event_fb=dict(component_max),largest_window_event_fb=max_packet,
                               zero_window_components=missing_components,
                               component_in_range_hypotheses=dict(in_range))
    sums = [np.sum([s[key] for s in sources.values()],axis=0) for key in
            ("rates_fb","gluon_jet_rates_fb","mc_covariance_fb2","integration_covariance_fb2")]
    count = sum(s["window_125"]["contributing_events"] for s in sources.values())
    missing = sorted(expected - set(sources))
    empty = sorted(name for name, source in sources.items() if not source["window_125"]["contributing_events"])
    upper = 0.
    for name in empty:
        source = sources[name]
        batches = [(sigma, m) for (sample, _), (sigma, _, m) in groups.items() if sample == name]
        # For independent positive unit-weight draws with no nonzero response,
        # P(no hits)=(1-p)^N. Every in-window event contributes at most its
        # trigger probability. If both specific ME photons must pass, both
        # .90 x .96 identification factors also bound the event response.
        if any(m.sumw != m.events or m.sumw2 != m.events for _,m in batches):
            raise ValueError("Zero-pass binomial limits require positive unit-weight collisions")
        # The dedicated matcher also enforces at most one stable representative
        # per LHE photon, even for an unusual branching shower record.
        two_me = name in {"tt_1a_1dec","tt_0a_2dec_same","tt_0a_2dec_split"}
        max_response = .988 * ((.90*.96)**2 if two_me else 1.)
        confidence_alpha = .05 / (2*len(expected))
        p_upper = -math.expm1(math.log(confidence_alpha) / source["attempted_events"])
        bound = max(sigma for sigma,_ in batches) * max_response * p_upper
        source["window_125"]["zero_pass_upper_fb_95"] = bound
        source["zero_pass_limit"] = dict(alpha=confidence_alpha, max_detector_probability=max_response,
                                         max_group_cross_section_fb=max(sigma for sigma,_ in batches),
                                         conditional_on_hard_cross_sections=True,
                                         entire_simulated_mass_range_empty=not source["component_in_range_hypotheses"])
        if not source["component_in_range_hypotheses"]:
            source["zero_pass_limit"]["entire_range_upper_fb_95"] = bound
        upper += bound
    full = make_report(*sums, count)
    full["window_125"]["zero_strata_upper_fb_95"] = upper
    return dict(schema=1, edges_gev=EDGES.tolist(), top_decay_model="matrix_element" if radiative else "shower",
                dileptonic_electron_enrichment=dileptonic,
                **full, complete_sources=region=="nominal" and not missing, generation_region=region,
                complete_region=not missing, missing_sources=missing, zero_window_sources=empty,
                sources=sources, shards=provenance)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--top-radiative", action="store_true", help="Replace shower decay photons by four disjoint ME strata")
    parser.add_argument("--generation-region", choices=("nominal","outside_nominal"), default="nominal")
    parser.add_argument("--top-dileptonic", action="store_true", help="Use exclusive dileptonic e-fake enrichment with the decay ME samples")
    args = parser.parse_args()
    report = summarize(args.campaign, args.top_radiative, args.generation_region, args.top_dileptonic)
    temp = args.output.with_name(args.output.name+".partial")
    temp.write_text(json.dumps(report, separators=(",",":"), allow_nan=False)+"\n")
    temp.replace(args.output)
    print(json.dumps({k:report[k] for k in ("window_125", "missing_sources", "zero_window_sources")}, indent=2))
