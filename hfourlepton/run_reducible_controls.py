"""Generate unbiased reducible controls, then shower disjoint LHE shards.

Example (runtime options after -- are passed to the ordinary four-lepton CLI):
  python hfourlepton/run_reducible_controls.py --run-tag controls_pilot \
    --cores 80 --workers 30 --events-per-shard 100 --nevents 1000 -- \
    --mg5-dir /path/to/MG5 --herwig-module herwig/stable

The index is a campaign inventory, not a nominal physics manifest. Shards of
one process share its inclusive cross section: pool their weighted responses
and source denominators; never add their inclusive cross sections.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import asdict, replace
import fcntl
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time

if str(Path(__file__).resolve().parents[1]) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import higher_order as higher
from hfourlepton.minimal_runtime import _encode, _decode
from hgammagamma import ho_shower_completion as completion
from hgammagamma import run_gammagamma_ho_campaign as diphoton

SAMPLES = ("bkg_zbb_ee", "bkg_zbb_mumu", "bkg_ttbar")
LIGHT_INCOMING = "g u c d s u~ c~ d~ s~"
STANDARD_INCOMING = "g u c d s b u~ c~ d~ s~ b~"
INCLUSIVE_HEAVY_FLAVOUR_CUTS = {
    "ptj": 0.0, "etaj": -1.0, "drjj": 0.0, "drjl": 0.0,
    "ptb": 0.0, "etab": -1.0, "drbb": 0.0, "drbl": 0.0,
}
PROMPT_MODES = (
    "/Herwig/Particles/W+/W+->nu_e,e+;",
    "/Herwig/Particles/W+/W+->nu_mu,mu+;",
    "/Herwig/Particles/W-/W-->nu_ebar,e-;",
    "/Herwig/Particles/W-/W-->nu_mubar,mu-;",
)


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, default=_encode, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def lhe_parts(source):
    """Stream complete events and the shared header/footer, rejecting grouping."""
    header, event, tail = [], None, []
    started = False
    with four.open_maybe_gzip(Path(source)) as stream:
        for line in stream:
            if re.search(r"<eventgroup\b", line):
                raise ValueError("grouped LHE events are unsupported")
            if re.search(r"<event\b", line):
                if event is not None:
                    raise ValueError("nested LHE event")
                if not started:
                    yield "header", "".join(header)
                    started = True
                elif any(x.strip() for x in tail):
                    raise ValueError("nonempty content between LHE events")
                tail = []
                event = [line]
            elif event is not None:
                event.append(line)
                if "</event>" in line:
                    yield "event", "".join(event)
                    event = None
            elif not started:
                header.append(line)
            else:
                tail.append(line)
    if event is not None or not started or "</LesHouchesEvents>" not in "".join(tail):
        raise ValueError("incomplete or empty LHE file")
    yield "footer", "".join(tail)


def split_lhe(source, destination, events_per_shard):
    """Partition all source records once; publish the completed directory atomically."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if events_per_shard < 1:
        raise ValueError("events_per_shard must be positive")
    identity = {"source_lhe": str(source), "source_sha256": four.sha256_file(source),
                "events_per_shard": events_per_shard}
    marker = destination / "partition.json"
    if destination.exists():
        if not marker.is_file():
            raise ValueError(f"incomplete shard directory: {destination}")
        saved = json.loads(marker.read_text())
        if saved["identity"] != identity:
            raise ValueError("shard source or partition changed; use a new tag")
        for shard in saved["shards"]:
            if four.sha256_file(Path(shard["lhe"])) != shard["sha256"]:
                raise ValueError("shard artifact changed")
        return saved
    count, header, footer = 0, None, None
    for kind, data in lhe_parts(source):
        if kind == "event":
            count += 1
        elif kind == "header":
            header = data
        else:
            footer = data
    # This adapter pools source populations directly and therefore requires
    # positive unit handler weights, as produced by the unweighted LO runs.
    weights = four.predecay_lhe_weight_summary(source, count)
    if weights["denominator_npositive"] != count or abs(weights["denominator_sumw"] - count) > 1.e-8 * count:
        raise ValueError("reducible shard pooling requires positive unit LHE handler weights")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".partition-", dir=destination.parent) as tmp:
        temporary = Path(tmp)
        output, shards, start, seen = None, [], 0, 0
        try:
            for kind, data in lhe_parts(source):
                if kind != "event":
                    continue
                if output is None:
                    filename = f"part-{len(shards):05d}.lhe.gz"
                    output = gzip.open(temporary / filename, "wt")
                    output.write(header)
                    start = seen
                output.write(data)
                seen += 1
                if seen - start == events_per_shard or seen == count:
                    output.write(footer)
                    output.close()
                    output = None
                    shards.append({"index": len(shards), "start": start, "stop": seen,
                                   "events": seen - start, "lhe": str(destination / filename),
                                   "sha256": four.sha256_file(temporary / filename)})
        finally:
            if output is not None:
                output.close()
        if seen != count or four.sha256_file(source) != identity["source_sha256"]:
            raise ValueError("source LHE changed during partition")
        record = {"identity": identity, "source_events": count,
                  "source_weight_summary": weights, "shards": shards}
        write_json(temporary / "partition.json", record)
        temporary.rename(destination)
    return record


def finite_card(text):
    """Prevent replay, preserving the LO shower and forced-decay settings."""
    line = "set /Herwig/EventHandlers/theLHReader:AllowedToReOpen No\n"
    if "AllowedToReOpen" in text:
        raise ValueError("LO template already specifies reader reopening; inspect before adapting")
    if text.count("saverun ") != 1:
        raise ValueError("expected one LO saverun command")
    text = text.replace("create ThePEG::LesHouchesFileReader theLHReader LesHouches.so",
                        "create ThePEG::LesHouchesFileReader theLHReader")
    text = text.replace("insert LesHouchesHandler:LesHouchesReaders[0] theLHReader",
                        "insert LesHouchesHandler:LesHouchesReaders 0 theLHReader")
    line += "set /Herwig/EventHandlers/theLHReader:Cuts /Herwig/Cuts/NoCuts\n"
    # The finite source ends with an expected EOF exception. Disable repository
    # diagnostic dumps on that path; ROOT finalization and the EOF audit remain.
    line += "set /Herwig/Generators/theGenerator:DumpPeriod -1\n"
    return text.replace("saverun ", line + "saverun ", 1)


def prompt_evidence(read_log, card_path, destination):
    """Parse independent runtime BR queries, preserving their printed precision."""
    pattern = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
    tokens = [line.strip() for line in Path(read_log).read_text().splitlines()
              if re.fullmatch(pattern, line.strip())]
    if len(tokens) != len(PROMPT_MODES):
        raise ValueError(f"expected exactly four numeric Herwig BR query outputs in {read_log}")
    probabilities = [float(token.replace("D", "e").replace("d", "e")) for token in tokens]
    if any(not 0 < value < 1 for value in probabilities):
        raise ValueError("invalid queried W branching probabilities")
    data = {"schema_version": 1, "profile": "w_emu",
            "branching_ratio_tokens": dict(zip(PROMPT_MODES, tokens)),
            "branching_ratios": dict(zip(PROMPT_MODES, probabilities)),
            "selected_branching_sums": {"24": sum(probabilities[:2]), "-24": sum(probabilities[2:])},
            "query_log": completion._artifact(read_log), "card": completion._artifact(card_path)}
    write_json(destination, data)
    return completion.prompt_branching_transform(destination)


def run_finite_shower(sample, lhe, cfg, *, native_decay_table=False):
    """Use the existing LO physics card with the validated finite-LHE audit."""
    work = four.shower_work_dir(sample, cfg)
    directory = work / "herwig"
    (directory / "events").mkdir(parents=True, exist_ok=True)
    library = None
    if sample.prompt_decay_profile or sample.is_hf_biased:
        library = directory / cfg.decay_reweighter_plugin.name
        if library.exists():
            if four.sha256_file(library) != four.sha256_file(cfg.decay_reweighter_plugin):
                raise ValueError("saved decay reweighter changed")
        else:
            shutil.copy2(cfg.decay_reweighter_plugin, library)
    card_path = directory / f"{sample.name}.in"
    temporary = directory / f"{sample.name}.in.prepared"
    metadata = four.write_herwig_input(sample, lhe, four.seed_for(f"{sample.name}:herwig", cfg),
                                     cfg.events_for(sample), temporary, cfg, reweighter_library=library)
    text = finite_card(temporary.read_text())
    temporary.unlink()
    if native_decay_table and not sample.is_hf_biased:
        if cfg.herwig_b_decays is None:
            raise ValueError("matched native control requires a native decay table")
        text = text.replace("saverun ", "cd /Herwig/Particles\n"
                            f"read {cfg.herwig_b_decays}\ncd /Herwig/Generators\nsaverun ", 1)
    if sample.prompt_decay_profile:
        if (sample.bias_of or sample.name) != "bkg_ttbar" or sample.prompt_decay_profile != "w_emu":
            raise ValueError("finite LO adapter supports only the validated ttbar W->e/mu selection")
        queries = "".join(f"get {mode}:BranchingRatio\n" for mode in PROMPT_MODES)
        text = text.replace("saverun ", queries + "saverun ", 1)
    if card_path.exists():
        if card_path.read_text() != text:
            raise ValueError("saved finite-LO shower card changed; use a new tag")
    else:
        card_path.write_text(text)
    marker = directory / "finite_lo_completion.json"
    source_sha = four.sha256_file(lhe)
    if marker.exists():
        saved = json.loads(marker.read_text())
        if saved.get("input_lhe_sha256") != source_sha:
            raise ValueError("saved shower source changed")
        completion.validate_completion(saved, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe": {"events": saved["lhe_events"]}}, directory=work)
    else:
        if list((directory / "events").glob("**/*.root")):
            raise ValueError(f"partial shower products at {directory}; inspect before recovery")
        count = sum(kind == "event" for kind, _ in lhe_parts(lhe))
        execution = completion.finite_shower_plan(sample.name, cfg.events_for(sample), count)
        write_json(directory / "shower_execution.json", execution)
        if higher.logged_command([cfg.herwig, "read", card_path], cfg, directory, directory / "read.log"):
            raise ValueError(f"Herwig read failed: {directory / 'read.log'}")
        transform = (prompt_evidence(directory / "read.log", card_path, directory / "prompt_decay_weights.json")
                     if sample.prompt_decay_profile else None)
        code = higher.logged_command([cfg.herwig, "run", f"{sample.name}.run", f"-N{execution['execution_events']}"],
                                     cfg, directory, directory / "run.log")
        roots = sorted((directory / "events").glob("**/*.root"))
        if higher.logged_command([cfg.mg5_python or sys.executable, "-c", diphoton.ROOT_METADATA_SCRIPT, *roots],
                                 cfg, directory, directory / "root-metadata.log"):
            raise ValueError("cannot audit shower ROOT output")
        root_info = json.loads((directory / "root-metadata.log").read_text().splitlines()[-1])
        hf_audit = None
        if sample.is_hf_biased:
            from hfourlepton.hf_decay_audit import build_audit
            hf_audit = build_audit(directory / "run.log", cfg.herwig_b_decays, card_path,
                                   cfg.decay_reweighter_source, library)
        # Finite-source execution must consume the entire disjoint input once.
        # Raw source weights remain the normalization denominator; the separate
        # transform audits the physical prompt-decay weights in the ROOT output.
        saved = completion.build_completion(directory, sample.name, cfg.events_for(sample), count, root_info,
                                            "source_exhausted" if code else "requested_events",
                                            source_weights=completion.source_weight_summary(lhe, count),
                                            source_weight_transform=transform,
                                            heavy_flavour_weight_audit=hf_audit)
        saved["input_lhe_sha256"] = source_sha
        saved["adapter_sha256"] = four.sha256_file(Path(__file__))
        saved["fingerprint"] = completion.completion_fingerprint(saved)
        completion.validate_completion(saved, {"sample": sample.name,
            "nevents_requested": cfg.events_for(sample), "lhe": {"events": count}})
        write_json(marker, saved)
    root_input = work / f"{sample.name}_hwsim_roots.input"
    four.write_root_input_list(directory / "events", root_input, cfg)
    metadata["shower_completion"] = saved
    return root_input, metadata


def worker(job):
    payload = json.loads(Path(job).read_text(), object_hook=_decode)
    if payload["adapter_sha256"] != four.sha256_file(Path(__file__)):
        raise ValueError("worker adapter changed after preparation")
    cfg = four.Config(**payload["config"])
    samples = four.selected_samples(cfg)
    four.require_inputs(cfg, samples)
    planned = four.build_manifest(cfg, samples)
    existing = four.load_resumable_manifest(cfg)
    if existing and existing["configuration_fingerprint"] != planned["configuration_fingerprint"]:
        raise ValueError("worker configuration changed")
    manifest = four.merge_manifest_samples(existing, planned)
    manifest["reducible_control_partition"] = payload["partition"]
    manifest["reducible_control_adapter_sha256"] = payload["adapter_sha256"]
    manifest["status"] = "running"
    four.write_manifest(manifest, cfg)
    original = four.run_herwig
    four.run_herwig = run_finite_shower
    started = time.monotonic()
    try:
        for sample in samples:
            four.run_sample(sample, samples, manifest, cfg)
        manifest["status"] = "complete"
    except BaseException as exc:
        manifest.update(status="failed", error=str(exc))
        raise
    finally:
        four.run_herwig = original
        manifest["last_attempt_wall_seconds"] = time.monotonic() - started
        four.write_manifest(manifest, cfg)


def generation_policy(inclusive, sample_name=None):
    if not inclusive:
        return {"name": "legacy_generator_defaults", "cut_overrides": {}}
    definitions = {name: LIGHT_INCOMING if name.startswith("bkg_zbb_") else STANDARD_INCOMING
                   for name in SAMPLES}
    policy = {"name": "inclusive_heavy_flavour_v2",
              "cut_overrides": dict(INCLUSIVE_HEAVY_FLAVOUR_CUTS),
              "incoming_partons_by_sample": definitions,
              "zbb_convention": "massive flavour-creation amplitudes with no incoming bottom; existing five-flavour LO PDFs as an approximation, not a complete four-flavour scheme",
              "pdf_and_model": "existing sm model, bottom mass and NNPDF40_lo_as_01180 unchanged"}
    if sample_name is not None:
        policy["incoming_partons"] = definitions[sample_name]
    return policy


def flavour_creation_process_card(text, sample_name):
    if not sample_name.startswith("bkg_zbb_"):
        return text
    rewritten, count = re.subn(r"(?m)^define p\s*=\s*[^\n]+$", "define p = " + LIGHT_INCOMING, text)
    if count != 1:
        raise ValueError("expected one proton definition for massive Zbb flavour creation")
    return rewritten


def rewrite_inclusive_run_card(text):
    """Set known MG5 cuts explicitly, including keys hidden in default cards."""
    found, output = set(), []
    for line in text.splitlines(keepends=True):
        match = re.match(r"^(\s*)([^=#]+?)(\s*=\s*)(\w+)\b(.*)$", line.rstrip("\r\n"))
        if match and match[4].lower() in INCLUSIVE_HEAVY_FLAVOUR_CUTS:
            key = match[4].lower()
            if key in found:
                raise ValueError(f"duplicate generation cut {key}")
            found.add(key)
            line = f"{match[1]}{INCLUSIVE_HEAVY_FLAVOUR_CUTS[key]:.1f}{match[3]}{match[4]}{match[5]}\n"
        output.append(line)
    missing = set(INCLUSIVE_HEAVY_FLAVOUR_CUTS) - found
    if missing:
        output.append("\n# Explicit inclusive heavy-flavour generation cuts (HiggsSSC).\n")
        for key in INCLUSIVE_HEAVY_FLAVOUR_CUTS:
            if key in missing:
                output.append(f"{INCLUSIVE_HEAVY_FLAVOUR_CUTS[key]:.1f} = {key}\n")
    return "".join(output)


def actual_generation_cuts(path):
    keys = {*INCLUSIVE_HEAVY_FLAVOUR_CUTS, "ptl", "etal", "drll", "mmll", "maxjetflavor"}
    result = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"\s*([^=#]+?)\s*=\s*(\w+)\b", line)
        if match and match[2].lower() in keys:
            key = match[2].lower()
            if key in result:
                raise ValueError(f"duplicate generation cut {key}")
            result[key] = match[1].strip()
    return result


def effective_generation_cuts(path):
    """Read executable Fortran assignments, independently of the user card."""
    keys = {*INCLUSIVE_HEAVY_FLAVOUR_CUTS, "ptl", "etal", "drll", "mmll", "maxjetflavor"}
    result = {}
    for line in Path(path).read_text().splitlines():
        match = re.match(r"\s*(\w+)\s*=\s*([^!]+)", line)
        if match and match[1].lower() in keys:
            key = match[1].lower()
            if key in result:
                raise ValueError(f"duplicate executable generation cut {key}")
            result[key] = float(match[2].strip().replace("D", "e").replace("d", "e"))
    return result


def audit_inclusive_generation(process_dir):
    process_dir = Path(process_dir)
    card, include = process_dir / "Cards/run_card.dat", process_dir / "Source/run_card.inc"
    explicit, effective = actual_generation_cuts(card), effective_generation_cuts(include)
    for key, expected in INCLUSIVE_HEAVY_FLAVOUR_CUTS.items():
        value = float(explicit[key].replace("D", "e").replace("d", "e")) if key in explicit else None
        if value != expected or effective.get(key) != expected:
            raise ValueError(f"inclusive heavy-flavour cut {key} differs between policy/card/executable")
    return {"include": str(include), "include_sha256": four.sha256_file(include),
            "effective_generation_cuts": effective}


def audit_incoming_partons(source, process_card, sample_name):
    expected = LIGHT_INCOMING if sample_name.startswith("bkg_zbb_") else STANDARD_INCOMING
    definitions = re.findall(r"(?m)^define p\s*=\s*([^\n]+)$", Path(process_card).read_text())
    if len(definitions) != 1 or definitions[0].split() != expected.split():
        raise ValueError("exported proton definition differs from inclusive generation policy")
    pairs, count = {}, 0
    for kind, block in lhe_parts(source):
        if kind != "event":
            continue
        rows = [line.split() for line in block.splitlines()[1:]
                if line.strip() and not line.lstrip().startswith(("<", "#"))]
        if len(rows) < 3 or any(len(row) != 13 or int(row[1]) != -1 for row in rows[1:3]):
            raise ValueError("expected two incoming LHE particles")
        ids = tuple(int(row[0]) for row in rows[1:3])
        if sample_name.startswith("bkg_zbb_") and any(abs(pid) == 5 for pid in ids):
            raise ValueError("incoming bottom found in massive Zbb flavour-creation source")
        if any(pid != 21 and abs(pid) not in (1, 2, 3, 4, 5) for pid in ids):
            raise ValueError("unexpected incoming parton in reducible source")
        key = ",".join(map(str, ids))
        pairs[key] = pairs.get(key, 0) + 1
        count += 1
    return {"process_card": str(process_card), "process_card_sha256": four.sha256_file(Path(process_card)),
            "incoming_partons": expected, "incoming_parton_audit": {"events": count, "pairs": pairs,
            "bottom_incoming_allowed": not sample_name.startswith("bkg_zbb_")}}


@contextmanager
def apply_generation_policy(inclusive):
    """Scope the policy to this coordinator; leave the ordinary runner untouched."""
    if not inclusive:
        yield
        return
    original_patch = four.patch_run_card
    original_include = four.generate_and_patch_run_card_include
    original_identity = four.production_stage_identity
    original_process_card = four.process_card_text
    def patch(card, **kwargs):
        original_patch(card, **kwargs)
        Path(card).write_text(rewrite_inclusive_run_card(Path(card).read_text()))
    def include(process_dir, cfg):
        original_include(process_dir, cfg)
        audit_inclusive_generation(process_dir)
    def identity(*args, **kwargs):
        sample = args[0] if args else kwargs["sample"]
        return {**original_identity(*args, **kwargs), "control_generation_policy": generation_policy(True, sample.name)}
    def process_card(sample, *args, **kwargs):
        return flavour_creation_process_card(original_process_card(sample, *args, **kwargs), sample.name)
    four.patch_run_card = patch
    four.generate_and_patch_run_card_include = include
    four.production_stage_identity = identity
    four.process_card_text = process_card
    try:
        yield
    finally:
        four.patch_run_card = original_patch
        four.generate_and_patch_run_card_include = original_include
        four.production_stage_identity = original_identity
        four.process_card_text = original_process_card


def run_campaign(args, forwarded):
    if args.cores < 1 or args.workers < 1 or args.workers > args.cores or args.events_per_shard < 1:
        raise ValueError("positive core/shard counts required and workers must not exceed cores")
    base = ["--order", "lo", "--sample-set", "reducible", "--run-samples", ",".join(SAMPLES),
            "--selection-profile", "ssc", "--detector-response", "ssc", "--ebeam", "20000",
            "--run-tag", args.run_tag, "--nb-core", str(args.cores), "--nevents", str(args.nevents)]
    if not any(x == "--seed-base" or x.startswith("--seed-base=") for x in forwarded):
        base += ["--seed-base", str(1 + four.stable_offset(args.run_tag))]
    cfg = four.parse_config(base + forwarded)
    if (cfg.order != "lo" or cfg.preset != "standard" or cfg.sample_set != "reducible"
        or cfg.detector_response != "ssc" or cfg.selection_profile != "ssc" or cfg.ebeam != 20000
        or cfg.heavy_flavour_bias != "unbiased" or cfg.force or cfg.dry_run
        or cfg.source_backend != "mg5" or cfg.run_tag != args.run_tag or cfg.nb_core != args.cores):
        raise ValueError("adapter requires unbiased standard LO SSC 40 TeV MG5 controls and no force/dry-run")
    samples = four.selected_samples(cfg)
    if tuple(s.name for s in samples) != SAMPLES:
        raise ValueError("adapter requires the two Zbb samples and ttbar only")
    cfg.run_dir.mkdir(parents=True, exist_ok=True)
    index_path = cfg.run_dir / "reducible_controls.json"
    identity = {"config": asdict(cfg), "events_per_shard": args.events_per_shard,
                "runner_sha256": four.sha256_file(Path(four.__file__))}
    # Preserve legacy pilot identities while pinning any explicit new policy.
    if args.inclusive_heavy_flavour:
        identity["generation_policy"] = generation_policy(True)
    adapter_sha = four.sha256_file(Path(__file__))
    with (cfg.run_dir / ".reducible_controls.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if index_path.exists():
            index = json.loads(index_path.read_text(), object_hook=_decode)
            if index["identity"] != identity:
                raise ValueError("campaign identity changed; use a new tag")
        else:
            index = {"schema_version": 1, "identity": identity, "status": "prepared", "samples": {},
                     "normalization": "Pool source-normalized responses once per process; do not sum shard rates."}
            write_json(index_path, index)
        if args.stage == "prepare":
            print(index_path)
            return
        four.require_inputs(cfg, samples)
        if args.stage in {"generate", "all"}:
            for sample in samples:
                print(f"Generating {sample.name} with {cfg.nb_core} cores", flush=True)
                with apply_generation_policy(args.inclusive_heavy_flavour):
                    source = four.obtain_production_lhe(sample, samples, cfg)
                process_dir = four.production_dir(sample, cfg) / "mg5_process"
                include_audit = (audit_inclusive_generation(process_dir)
                                 if args.inclusive_heavy_flavour else {})
                if args.inclusive_heavy_flavour:
                    include_audit.update(audit_incoming_partons(source,
                        four.production_dir(sample, cfg) / "cards/process_card.dat", sample.name))
                partition = split_lhe(source, cfg.run_dir / "partitions" / sample.name, args.events_per_shard)
                run_card = four.production_dir(sample, cfg) / "mg5_process/Cards/run_card.dat"
                old_shards = {s["index"]: s for s in index["samples"].get(sample.name, {}).get("shards", [])}
                record = {**partition["identity"], "source_events": partition["source_events"],
                          "production_seed": four.seed_for(sample.production_group_id, cfg),
                          "generation_policy": generation_policy(args.inclusive_heavy_flavour, sample.name),
                          **include_audit,
                          "cross_section_pb": four.parse_lhe_init_cross_section_pb(source),
                          "generation_run_card": str(run_card), "generation_run_card_sha256": four.sha256_file(run_card),
                          "actual_generation_cuts": actual_generation_cuts(run_card), "shards": []}
                for shard in partition["shards"]:
                    tag = f"{cfg.run_tag}-{sample.name}-{shard['index']:05d}"
                    worker_cfg = replace(cfg, run_tag=tag, nb_core=1, nevents=shard["events"],
                        sample_nevents={sample.name: shard["events"]}, run_samples=sample.name,
                        source_backend="external-lhe", external_lhe={sample.name: Path(shard["lhe"])},
                        external_perturbative_order="LO", seed_base=1 + four.stable_offset(tag), postprocess_only=True)
                    job = cfg.run_dir / "jobs" / f"{sample.name}-{shard['index']:05d}.json"
                    partition_info = {**partition["identity"], **shard,
                                      "generation_policy": generation_policy(args.inclusive_heavy_flavour, sample.name)}
                    if not job.exists():
                        write_json(job, {"config": asdict(worker_cfg), "partition": partition_info,
                                         "adapter_sha256": adapter_sha})
                    record["shards"].append({**shard, "run_tag": tag, "manifest": str(worker_cfg.manifest_path),
                        "job": str(job), "log": str(job.with_suffix(".log")),
                        "status": old_shards.get(shard["index"], {}).get("status", "prepared")})
                index["samples"][sample.name] = record
                write_json(index_path, index)
            index["status"] = "generated"
            write_json(index_path, index)
        if args.stage == "generate":
            return
        if set(index["samples"]) != set(SAMPLES):
            raise ValueError("generate all three hard samples before showering")
        # Build shared read-only executables once before launching workers.
        four.run_runtime_command(["make", "-C", cfg.analysis_build_dir, cfg.analysis_target], cfg)
        four.build_decay_reweighter(cfg)
        jobs = [s for record in index["samples"].values() for s in record["shards"]]
        if index.get("worker_adapter_sha256") != adapter_sha:
            if any(Path(shard["manifest"]).exists() for shard in jobs):
                raise ValueError("worker adapter changed after showers started; inspect existing artifacts")
            index.setdefault("worker_adapter_history", []).append({"sha256": adapter_sha,
                "reason": "pin current adapter before first shower", "time_unix": time.time()})
            index["worker_adapter_sha256"] = adapter_sha
            for shard in jobs:
                payload = json.loads(Path(shard["job"]).read_text(), object_hook=_decode)
                payload["adapter_sha256"] = adapter_sha
                write_json(shard["job"], payload)
        index.update(status="showering", workers=args.workers)
        write_json(index_path, index)
        stop = threading.Event()
        def launch(shard):
            if stop.is_set():
                return None
            with Path(shard["log"]).open("a") as out:
                result = subprocess.run([cfg.mg5_python or sys.executable, Path(__file__).resolve(), "--worker", shard["job"]],
                    stdout=out, stderr=subprocess.STDOUT,
                    env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                         "PYTHONUNBUFFERED": "1"})
            if result.returncode:
                stop.set()
            return result.returncode
        failures = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(launch, shard): shard for shard in jobs}
            for future in as_completed(futures):
                shard = futures[future]
                code = future.result()
                shard["status"] = "prepared" if code is None else "complete" if code == 0 else "failed"
                if code:
                    failures.append(shard["log"])
                write_json(index_path, index)
                print(f"{shard['run_tag']}: {shard['status']}", flush=True)
        index["status"] = "failed" if failures else "complete"
        write_json(index_path, index)
        if failures:
            raise RuntimeError(f"{len(failures)} shard(s) failed; inspect {failures[0]}")
        print(f"Completed control index: {index_path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--run-tag")
    parser.add_argument("--stage", choices=("prepare", "generate", "shower", "all"), default="all")
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--events-per-shard", type=int, default=5000)
    parser.add_argument("--nevents", type=int, default=10000)
    parser.add_argument("--inclusive-heavy-flavour", action="store_true",
                        help="remove hard jet cuts; use massive Zbb flavour creation without incoming bottom (existing 5F PDFs); audit source and executable cuts")
    args, forwarded = parser.parse_known_args(argv)
    if args.worker:
        if forwarded:
            parser.error("worker does not accept extra arguments")
        worker(args.worker)
    else:
        if not args.run_tag:
            parser.error("--run-tag is required")
        run_campaign(args, forwarded[1:] if forwarded[:1] == ["--"] else forwarded)


if __name__ == "__main__":
    main()
