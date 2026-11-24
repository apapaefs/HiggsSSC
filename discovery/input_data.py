"""Load the signed-weight spectra in campaign order and pool independent hard events."""
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile

import numpy as np
from scipy.optimize import brentq

from likelihood import EDGES

DATA = Path(__file__).resolve().parent / "data"
SIGNAL = "signal_gg_h_aa"
BACKGROUNDS = ("bkg_prompt_aa", "bkg_loop_gg_aa", "bkg_gamma_j", "bkg_dijet_double_fake", "bkg_dy_ee")
SCENARIOS = {"central": "observed_heavy_flavour", "upper": "upper_heavy_flavour"}


def load_archive():
    path = DATA / "spectra.tar.gz"
    expected = (DATA / "SHA256SUMS").read_text().split()[0]
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError("Input archive checksum mismatch")
    with tarfile.open(path) as tar:
        names = tar.getnames()
        manifest = json.load(tar.extractfile("manifest.json"))
        hashes = manifest["files"]
        if (manifest["schema"] != 1 or len(set(names)) != len(names)
                or set(names) != {"manifest.json", *hashes}):
            raise ValueError("Unexpected archive contents")
        objects = {}
        for name, digest in hashes.items():
            raw = tar.extractfile(name).read()
            if hashlib.sha256(raw).hexdigest() != digest:
                raise ValueError(f"Corrupt archive member {name}")
            objects[name] = json.loads(raw)
    fit, reach = objects["records/final_fit.json"], objects["records/final_reach.json"]
    packets = sorted(name for name in hashes if name.startswith("packets/"))
    packets.sort(key=lambda name: int(name.split("/")[2].split("-")[0]))
    if (hashes["records/final_fit.json"] != reach["fit_sha256"]
            or hashes["baseline.json"] != fit["reference_sha256"]
            or [hashes[p] for p in packets] != [r["sha256"] for r in fit["independent_new_inputs"]]):
        raise ValueError("Archive does not reproduce the signed-off production manifest")
    return objects, packets


def source_keys(packet, sample):
    provenance = packet.get("provenance", {})
    keys = set()
    for identifier in (packet["independent_id"], packet.get("source_lhe_sha256"),
                       provenance.get("sha256") if provenance.get("parent_lhe") else None):
        if identifier is not None:
            is_hash = re.fullmatch("[0-9a-fA-F]{64}", identifier)
            keys.add(("hard", identifier.lower()) if is_hash else (sample, identifier))
    original = dict(packet)
    if "hard_seed" in provenance:
        original["hard_seed"] = provenance["hard_seed"]
    for part in [original, *packet.get("shards", [])]:
        seeds = list({int(part[k]) for k in ("seed", "hard_seed", "generator_seed") if k in part})
        for k in ("auxiliary_generator_seeds", "continuation_generator_seeds"):
            seeds.extend(int(seed) for seed in part.get(k, []))
        for seed in seeds:
            key = (sample, seed)
            if key in keys:
                raise ValueError(f"Repeated seed in {packet['independent_id']}")
            keys.add(key)
    return keys


def pool(objects, packet_names):
    base = objects["baseline.json"]["datasets"]
    allowed = {(d["sample"], d["stratum"]) for d in base}
    members, used, added = defaultdict(list), set(), defaultdict(int)
    for name, packet in [(None, d) for d in base] + [(p, objects[p]) for p in packet_names]:
        sample, stratum = packet["sample"], packet.get("stratum", "full")
        if sample == "bkg_dijet_double_fake" and stratum in {
                "inclusive_minbias", "inclusive_minimum_bias", "minimum_bias", "inclusive minimum bias"}:
            stratum = "full"
        group = (sample, stratum)
        keys = source_keys(packet, sample)
        if group not in allowed or used & keys or "first_source" in packet or "last_source_exclusive" in packet:
            raise ValueError(f"Unexpected or overlapping exposure in {name or sample}")
        used.update(keys)
        n = int(packet["generated_events"])
        edges = packet.get("edges_gev", packet.get("mass_bin_edges_GeV"))
        rates = np.asarray(packet.get("rate_fb", packet.get("rates_fb")), dtype=float)
        cov = np.asarray(packet["covariance_fb2"], dtype=float)
        if (n <= 0 or n != packet["generated_events"] or not np.array_equal(edges, EDGES)
                or rates.shape != (55,) or cov.shape != (55, 55)
                or not np.isfinite(rates).all() or not np.isfinite(cov).all()
                or not np.allclose(cov, cov.T)
                or np.linalg.eigvalsh(cov)[0] < -1.e-9 * max(1., cov.max())):
            raise ValueError(f"Invalid spectrum in {name or sample}")
        members[group].append((n, rates, cov))
        if name and name.startswith("packets/final/"):
            added["/".join(group)] += n
    rates, covs, total = {}, {}, {}
    for group, chunks in members.items():
        sample, _ = group
        n = sum(chunk[0] for chunk in chunks)
        rates[sample] = rates.get(sample, 0.) + sum(h * (size / n) for size, h, _ in chunks)
        covs[sample] = covs.get(sample, 0.) + sum(c * (size / n)**2 for size, _, c in chunks)
        total["/".join(group)] = n
    if set(rates) != {*BACKGROUNDS, SIGNAL}:
        raise ValueError("Missing process in the published diphoton spectrum")
    return rates, covs, total, dict(added)


def fourlepton_rates(extended, reducible):
    conditions = reducible["discovery_conditions"]
    if extended["threshold"] != conditions["frozen_model_threshold"]:
        raise ValueError("The four-lepton evaluation and the zero-count bound use different selections")
    samples = extended["extended_test"]["samples"]
    def rate(group):
        return sum(s["sum_weight"] for s in samples if s["category"] == group) / conditions["reference_test_luminosity_fb_inverse"]
    signal, central = rate("signal"), rate("irreducible")
    inputs = reducible["old_zero_count_inputs"]
    if any(row["selected_positive_response_events"] for row in inputs):
        raise ValueError("Joint zero-count bound requires zero observed responses")
    # Maximize the all-zero probability at fixed sum of the three unknown rates.
    ordered = sorted(inputs, key=lambda r: r["raw_source_events"] / r["conditional_rate_fb"])
    def log_zero(total_rate):
        size = scale = 0.
        for i, row in enumerate(ordered, 1):
            size += row["raw_source_events"]
            scale += row["conditional_rate_fb"]
            multiplier = size / (scale - total_rate)
            if i < len(ordered) and multiplier > ordered[i]["mc_exposure_fb_inverse"]:
                continue
            if i == 1:
                return size * math.log1p(-total_rate / scale)
            return sum(p["raw_source_events"] * math.log(p["raw_source_events"] /
                       (multiplier * p["conditional_rate_fb"])) for p in ordered[:i])
        raise ValueError("Unbounded heavy-flavour limit")
    upper = brentq(lambda rate: log_zero(rate) - math.log(.05), 0.,
                   -math.log(.05) / ordered[0]["mc_exposure_fb_inverse"])
    return signal, {"central": central, "upper": central + upper}
