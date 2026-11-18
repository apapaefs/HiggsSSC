"""Portable, source-grouped analysis products for disjoint heavy-flavour control shards.

The full detector and shower trees remain available in the shared release. The
compact form contains the complete finally selected hypotheses (all kinematics),
and raw signed moments for each earlier cut before cross-section normalization.
"""
from __future__ import annotations

from collections import defaultdict
import copy
from dataclasses import fields
import json
import math
from pathlib import Path

from hfourlepton import make_fourlepton_report as report
from hfourlepton import shared_events as shared

PROCESSES = ("bkg_zbb_ee", "bkg_zbb_mumu", "bkg_ttbar")
STATS = tuple(f.name for f in fields(report.WeightStats))


def raw_stats(row):
    return {key: row[key] for key in STATS}


def sum_stats(rows):
    result = report.WeightStats()
    for row in rows:
        result.merge(report.WeightStats(**raw_stats(row)))
    return result


def compact_shard(path):
    from hfourlepton.report_reader import load_report_records
    bundle = report.load_manifest(path)
    sample, = bundle.samples
    stage = report._final_stage(bundle.cut_stages)
    records = load_report_records(sample.output_root, selected_mask=stage.required_mask)
    result = report.analyze_sample(sample, records, bundle.cut_stages, 1., (120., 130.), 1.e-9)
    if bundle.audit.errors or result["audit"]["errors"] or not result["included_in_physics_totals"]:
        raise ValueError(f"{path}: {bundle.audit.errors + result['audit']['errors']}")
    if (sample.name not in PROCESSES or sample.category != "reducible" or sample.decay_bias_enabled
            or sample.raw["selection_version"] != "ssc_v1" or sample.perturbative_order != "LO"):
        raise ValueError("compact publication requires the unbiased SSC heavy-flavour controls")
    selected = [r for r in records.four_lepton if report.passes_stage(r, stage)]
    src = sorted(int(r["source_index"]) for r in records.source_events)
    partition = json.loads(path.read_text())["reducible_control_partition"]
    if not src or src != list(range(src[0], src[0] + len(src))) or len(src) > partition["events"]:
        raise ValueError(f"noncontiguous saved event indices in {path}")
    if any(int(row["source_index"]) < src[0] or int(row["source_index"]) > src[-1] for row in selected):
        raise ValueError("a selected hypothesis lacks its source")
    if result["selected"]["raw"]["entries"] != len({r["source_index"] for r in selected}):
        raise ValueError("selected population differs from the full report")
    for row in selected:
        row["source_index"] += partition["start"]
    cutflows = {"all": result["cutflow"], **{k: v["cutflow"] for k,v in result["channels"].items()}}
    return dict(name=sample.name, profile=sample.detector_profile, label=sample.label,
        order=sample.perturbative_order, source_kind=sample.source_kind, category=sample.category,
        source_range=[partition["start"] + src[0], partition["start"] + src[-1] + 1],
        partition={k: partition[k] for k in ("index","source_sha256","start","stop","events","sha256")},
        cuts=[dict(name=s.name, bit=s.bit, required_mask=s.required_mask) for s in bundle.cut_stages],
        generation_seeds=sample.raw["seeds"], native_pb=sample.native_cross_section_pb,
        target_pb=sample.target_cross_section_pb, normalization=result["normalization"],
        generated=result["generated"], selected=selected, selected_weight=result["selected"]["raw"],
        cutflows={k: [raw_stats(row) for row in values] for k,values in cutflows.items()},
        audit=result["audit"], detector_audit=result["detector_audit"],
        closure=result["response_closure"], optional=result["optional_weights"],
        manifest=str(path), manifest_sha256=shared.digest(path), root_sha256=sample.raw["output_root_sha256"],
        shower_fingerprint=sample.raw["shower_completion"]["fingerprint"])


def pool(shards, total):
    if not shards:
        raise ValueError("empty pool")
    first = shards[0]
    selected, source_ranges, inventory = [], [], []
    previous = 0
    grouped = defaultdict(list)
    for i, shard in enumerate(shards):
        part = shard["partition"]
        if (part["index"] != i or part["start"] != previous or part["stop"]-previous != part["events"]
                or part["events"] != shard["normalization"]["denominator_event_count"]
                or part["source_sha256"] != first["partition"]["source_sha256"]
                or any(shard[k] != first[k] for k in ("name","profile","cuts","native_pb","target_pb","order",
                                                        "source_kind","category","label"))):
            raise ValueError("incompatible shards or incomplete/overlapping parent partition")
        previous = part["stop"]
        if source_ranges and shard["source_range"][0] < source_ranges[-1][1]:
            raise ValueError("pooled independent source ids overlap")
        source_ranges.append(shard["source_range"])
        selected.extend(shard["selected"])
        inventory.append({k:shard[k] for k in ("partition","manifest","manifest_sha256","root_sha256",
                                                "shower_fingerprint","generation_seeds")})
        for key, rows in shard["cutflows"].items():
            grouped[key].append(rows)
    if previous != total:
        raise ValueError("pooled exposure differs from complete hard parent")
    cutflows = {channel:[sum_stats(rows).as_dict() for rows in zip(*stages)] for channel,stages in grouped.items()}
    norm = dict(first["normalization"])
    for key in ("denominator_sumw","denominator_sumw2","denominator_sumabsw","denominator_event_count",
                "denominator_npositive","denominator_nnegative","denominator_nzero"):
        norm[key] = sum(s["normalization"][key] for s in shards)
    scale = norm["scale_pb_per_weight"] = first["native_pb"] / norm["denominator_sumw"]
    selected_stats = report.grouped_weight_stats(selected)
    check = sum_stats(s["selected_weight"] for s in shards)
    if check.entries != selected_stats.entries or not math.isclose(check.sumw2, selected_stats.sumw2, rel_tol=1e-12):
        raise ValueError("selected source weight variance changed while pooling")
    norm.update(selected_cross_section_pb=selected_stats.sumw * scale,
                selected_cross_section_sumw2_pb2=selected_stats.sumw2 * scale * scale)
    closure = {key:sum(s["closure"][key] for s in shards) for key in
               ("source_events","evaluated","missing","failed","direct_delta_disagreements")}
    closure.update(ok=all(s["closure"]["ok"] for s in shards), tolerance=first["closure"]["tolerance"],
        max_abs_delta=max(s["closure"]["max_abs_delta"] for s in shards),
        mean_abs_delta=sum(s["closure"]["mean_abs_delta"] * s["closure"]["evaluated"] for s in shards)/closure["evaluated"])
    optional = copy.deepcopy(first["optional"])
    if any(s["optional"]["variation_names"] != optional["variation_names"] for s in shards):
        raise ValueError("optional-weight definitions differ between shards")
    for key in ("source_events","events_with_optional_branches","events_with_variations","synthesized_name_events"):
        if optional[key] is not None:
            optional[key] = sum(s["optional"][key] for s in shards)
    for i, item in enumerate(optional["variations"]):
        for key in ("generated_raw","selected_raw"):
            item[key] = sum_stats(s["optional"]["variations"][i][key] for s in shards).as_dict()
    warnings = sorted({w for shard in shards for w in shard["audit"]["warnings"]})
    return dict(schema_version=1, kind="pooled_fourlepton_controls", name=first["name"],
        detector_profile=first["profile"], label=first["label"], source_kind=first["source_kind"],
        native_pb=first["native_pb"], target_pb=first["target_pb"], generated=sum_stats(s["generated"] for s in shards).as_dict(),
        normalization=norm, cutflows=cutflows, cuts=first["cuts"], selected=selected, source_ranges=source_ranges,
        source_sha256=first["partition"]["source_sha256"], partitions=inventory, optional_weights=optional,
        response_closure=closure, audit=dict(ok=True, errors=[], warnings=warnings))


def load(manifest, stages, *, published_defaults=False):
    import gzip
    import hashlib
    path = Path(manifest)
    raw = json.loads(path.read_text())
    refs = raw.get("pooled_backgrounds")
    root = path.parent
    packed = refs is None and published_defaults
    if packed:
        from hfourlepton import student_backgrounds
        root = student_backgrounds.DIRECTORY
        refs = student_backgrounds.controls_for(raw, root)
    if refs is None:
        return []
    if set(refs) != set(PROCESSES):
        raise ValueError("incomplete pooled background inventory")
    results = []
    expected = [dict(name=s.name,bit=s.bit,required_mask=s.required_mask) for s in stages]
    for name in PROCESSES:
        info = refs[name]
        data_path = shared.safe_path(root, info["path"])
        data = data_path.read_bytes()
        if hashlib.sha256(data).hexdigest() != info["sha256"]:
            raise ValueError(f"changed pooled background {name}")
        if packed:
            data = gzip.decompress(data)
            if hashlib.sha256(data).hexdigest() != info["public_sha256"]:
                raise ValueError(f"packed control differs from the public release: {name}")
        sample = json.loads(data)
        if (sample["schema_version"] != 1 or sample["kind"] != "pooled_fourlepton_controls"
                or sample["name"] != name or sample["detector_profile"] != raw["detector_profile"]
                or sample["cuts"] != expected or sample["audit"]["errors"]):
            raise ValueError(f"incompatible pooled background {name}")
        results.append((sample, data_path, info["sha256"]))
    return results


def sample_report(sample, path, stages, luminosity, signal_region, tolerance, mass_range, mass_width, observables):
    norm = sample["normalization"]
    rows = sample["selected"]
    audit = copy.deepcopy(sample["audit"])
    if sample["response_closure"]["max_abs_delta"] > tolerance:
        audit["errors"].append("pooled control exceeds the requested response-closure tolerance")
        audit["ok"] = False
    def yields(channel, subset):
        low, high = signal_region
        flow = []
        for stage, raw in zip(stages, sample["cutflows"][channel]):
            value = report._yield_dict(report.WeightStats(**raw_stats(raw)), norm, luminosity)
            flow.append(dict(**raw, name=stage.name, label=stage.label, bit=stage.bit,
                cross_section_pb=value["cross_section_pb"], cross_section_sumw2_pb2=value["cross_section_sumw2_pb2"],
                yield_events=value["events"], yield_sumw2=value["events_sumw2"],
                yield_stat_uncertainty=value["stat_uncertainty_events"], normalization_status=value["status"]))
        return dict(cutflow=flow, selected=report._yield_dict(report.grouped_weight_stats(subset),norm,luminosity),
            signal_region=report._yield_dict(report.grouped_weight_stats(subset,lambda r: low<r["m4l"]<high),norm,luminosity),
            mass_spectrum=report.mass_spectrum(subset,norm,luminosity,mass_range,mass_width))
    identity = dict(input_lhe_sha256=sample["source_sha256"],source_lhe_sha256=sample["source_sha256"],
        production_group_id=sample["name"],nevents_requested=int(norm["denominator_event_count"]),
        source_events_read=sample["generated"]["entries"],seeds={key:[p["generation_seeds"][key] for p in sample["partitions"]]
            for key in ("production","herwig","analysis")})
    result = dict(name=sample["name"],label=sample["label"],category="reducible",channel="inclusive",
        detector_profile=sample["detector_profile"],perturbative_order="LO",source_kind=sample["source_kind"],
        production_group_id=sample["name"],comparison_identity=identity,decay_bias={"enabled":False},
        included_in_physics_totals=True, exclusion_reason=None,exclusion_reasons=[],manifest=str(path),
        root_file=None, shard_manifests=[p["manifest"] for p in sample["partitions"]],native_cross_section_pb=sample["native_pb"],
        target_cross_section_pb=sample["target_pb"],weight_scale=1.,generated=sample["generated"],normalization=norm,
        channels={ch:yields(ch,[r for r in rows if report.normalize_channel(r["channel"])==ch]) for ch in sample["cutflows"] if ch!="all"},
        response_closure=sample["response_closure"],optional_weights=sample["optional_weights"],
        detector_audit=dict(ok=True,shards=len(sample["partitions"]),manifest_profile=sample["detector_profile"]),
        mass_resolution=None,audit=audit,**yields("all",rows))
    if observables:
        from hfourlepton.observables import sample_observables
        result["observables"] = sample_observables(rows,norm,signal_region)
    return result
