"""Publish the SSC heavy-flavour controls and reconstruct the same showers with perfect response."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hfourlepton import shared_events as shared

PROCESSES = ("bkg_zbb_ee", "bkg_zbb_mumu", "bkg_ttbar")


def read(path):
    return json.loads(Path(path).read_text())


def copy_shard(job):
    destination, name, shard, receipt, detector = job
    if shared.digest(shard["manifest"]) != receipt["manifest_sha256"]:
        raise ValueError(f"changed control manifest: {shard['manifest']}")
    original = read(shard["manifest"])
    src, = original["samples"]
    if (src["name"] != name or src["category"] != "reducible" or src["decay_bias"]["enabled"]
            or original["detector_profile"] != "ssc" or src["selection_version"] != "ssc_v1"
            or src["input_lhe_sha256"] != shard["sha256"]
            or src["output_root_sha256"] != receipt["output_root_sha256"]):
        raise ValueError(f"incorrect unbiased control: {name}/{shard['index']}")
    directory = destination / "controls" / name / f"{shard['index']:05d}"
    with shared.Publisher(directory, dict(source=receipt["manifest_sha256"])) as pub:
        ssc = pub.manifest(shard["manifest"], "ssc.json", "raw", expected_order="lo")
        paired = read(ssc)
        sample, = paired["samples"]
        ssc_summary = read(sample["analysis_summary"])
        out = directory / "raw/analysis/perfect" / name
        out.mkdir(parents=True, exist_ok=True)
        command = src["analyzer_cli"].copy()
        for key, value in (("--response-profile", "perfect"), ("--input-list", ssc_summary["input"]),
                           ("--output-dir", str(out)), ("--optional-weight-names-file",
                            sample["optional_weight_names_file"]["path"])):
            command[command.index(key) + 1] = value
        summary_path = out / f"{name}_{paired['run_tag']}_perfect.summary.json"
        if not summary_path.is_file():
            with (out / "analyzer.log").open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        summary = read(summary_path)
        for key in ("events_read", "sample", "category", "requested_channel", "seed",
                    "weight_scale", "selection_version", "generated_raw", "generated_scaled"):
            if summary[key] != ssc_summary[key]:
                raise ValueError(f"perfect reconstruction changed {name}/{key}")
        if summary["response_profile"] != "perfect":
            raise ValueError("wrong detector response on paired output")
        root = Path(summary["output_root"])
        sample.update(root_file=str(root), output_root=str(root), output_root_sha256=shared.digest(root),
                      analysis_summary=str(summary_path), analysis_summary_sha256=shared.digest(summary_path),
                      weight_sums={k: summary[k] for k in ("generated_raw", "generated_scaled","selected_scaled")})
        paired.update(detector_profile="perfect", response_mode="perfect", detector=detector)
        pub.json("perfect.json", paired)
        for path in (root, summary_path, out / "analyzer.log"):
            pub.register(path, shared.digest(path))
        shared.save(directory / "files.json", pub.files)
    print(f"Copied and paired {name}/{shard['index']:05d}, {summary['selected_rows']} selected", flush=True)


def audit_shard(job):
    from hfourlepton.pooled_controls import compact_shard
    destination, name, shard, receipt, _ = job
    directory = destination / "controls" / name / f"{shard['index']:05d}"
    if not (directory / "files.json").exists():
        raise ValueError(f"missing completed public copy: {directory}; finish the copy phase first")
    for profile in ("ssc", "perfect"):
        manifest = directory / f"{profile}.json"
        path = directory / f"{profile}.compact.json"
        if path.is_file() and read(path)["manifest_sha256"] == shared.digest(manifest):
            continue
        summary = compact_shard(manifest)
        shared.save(path, summary)
    print(f"Audited both profiles for {name}/{shard['index']:05d}", flush=True)


def jobs(args):
    receipt = read(args.receipt)
    if shared.digest(args.index) != receipt["index_sha256"]:
        raise ValueError("generation index differs from the completed unbiased control receipt")
    by_path = {r["manifest"]: r for r in receipt["inputs"]}
    index = read(args.index)
    base = read(args.base / "dataset.json")
    detector = read(args.base / base["manifests"]["perfect"])["detector"]
    inventory = []
    if set(index["samples"]) != set(PROCESSES):
        raise ValueError("incorrect control process inventory")
    for name in PROCESSES:
        parent = index["samples"][name]
        offset = 0
        for i, shard in enumerate(parent["shards"]):
            if (shard["index"] != i or shard["start"] != offset or shard["stop"]-offset != shard["events"]
                    or shard["status"] != "complete"):
                raise ValueError("noncontiguous or incomplete control production")
            offset = shard["stop"]
            inventory.append((args.destination, name, shard, by_path.pop(shard["manifest"]), detector))
        if offset != parent["source_events"]:
            raise ValueError("partition inventory does not cover the whole parent")
    if by_path:
        raise ValueError("unaccounted control manifests")
    return inventory


def seal(args, items):
    from datetime import datetime, timezone
    import os
    from hfourlepton import pooled_controls

    root = args.destination.resolve()
    base_path = args.base.resolve()
    if root.parent != base_path.parent or root == base_path:
        raise ValueError("the base and extension must be sibling version directories")
    base = shared.load_dataset(base_path)
    parents = read(args.index)["samples"]
    files = {}
    refs = {"ssc": {}, "perfect": {}}
    added = []
    for name in PROCESSES:
        halves = {profile:[] for profile in refs}
        for _, current, shard, receipt, _ in items:
            if current != name:
                continue
            directory = root / "controls" / name / f"{shard['index']:05d}"
            subfiles = read(directory / "files.json")
            for relative, info in subfiles.items():
                path = shared.safe_path(directory, relative)
                if path.stat().st_size != info["size"]:
                    raise ValueError(f"changed input after publication: {path}")
                files[str(path.relative_to(root))] = info
            for profile in refs:
                data = read(directory / f"{profile}.compact.json")
                if data["manifest_sha256"] != shared.digest(directory / f"{profile}.json"):
                    raise ValueError("full-tree audit no longer matches the public manifest")
                part = data["partition"]
                if (part["source_sha256"] != parents[name]["source_sha256"]
                        or any(part[k] != shard[k] for k in ("index","start","stop","events","sha256"))):
                    raise ValueError("audited sample is not the indexed parent hard-event slice")
                if data["root_sha256"] != read(data["manifest"])["samples"][0]["output_root_sha256"]:
                    raise ValueError("compact input differs from the public tree")
                if profile == "ssc" and data["root_sha256"] != receipt["output_root_sha256"]:
                    raise ValueError("SSC tree differs from the original control")
                halves[profile].append(data)
        for ssc, perfect in zip(halves["ssc"], halves["perfect"]):
            if any(ssc[key] != perfect[key] for key in ("shower_fingerprint","partition","source_range","generation_seeds","generated")):
                raise ValueError("SSC and perfect controls have different source populations")
        for profile, pieces in halves.items():
            result = pooled_controls.pool(pieces, parents[name]["source_events"])
            relative = f"controls/{name}/{profile}.json"
            shared.save(root / relative, result)
            refs[profile][name] = dict(path=relative, sha256=shared.digest(root / relative))
        added.append(dict(name=name, category="reducible", channel="inclusive",
                          requested=parents[name]["source_events"], saved=result["generated"]["entries"]))
    # Preserve all original hard parents, generation cards, indices and versioned
    # input evidence, alongside the 300 independently relocatable event shards.
    proof_dir = root / "controls/provenance"
    with shared.Publisher(proof_dir, dict(index_sha256=shared.digest(args.index))) as pub:
        pub.copy_file(args.index,"reducible_controls.json")
        pub.copy_file(args.receipt,"full_receipt.json")
        original = read(items[0][2]["manifest"])
        executable = Path(original["samples"][0]["analyzer_cli"][0])
        archived_executable = pub.copy_file(executable, Path("analyzer") / executable.name)
        for key in ("analyzer_source","analyzer_ancestry_header","analyzer_makefile"):
            artifact = original["tracked_inputs"][key]
            pub.copy_file(artifact["path"], Path("analyzer") / Path(artifact["path"]).name, artifact["sha256"])
        pub.json("perfect_reconstruction.json", dict(
            analyzer=pub.artifact(archived_executable), software=original["software"],
            output_profile="perfect", selection_profile="ssc", selection_version="ssc_v1",
            command_source="Original per-shard samples[0].analyzer_cli",
            argument_changes=["--response-profile perfect", "--input-list points to the published same-shower input list",
                              "--optional-weight-names-file points to the published mapping",
                              "--output-dir points to the shard's perfect reconstruction folder"],
            shared_settings=["analysis seed","cut mask","weight scale","muon resolution scale","sample","category","channel"]))
        pub.copy_file(__file__, "analyzer/publish_reducible_events.py")
        pub.copy_file(pooled_controls.__file__, "analyzer/pooled_controls.py")
        for name in PROCESSES:
            parent = parents[name]
            for key, sha in (("source_lhe","source_sha256"),("process_card","process_card_sha256"),
                             ("generation_run_card","generation_run_card_sha256"),("include","include_sha256")):
                pub.copy_file(parent[key], Path(name) / Path(parent[key]).name, parent[sha])
        for relative, info in pub.files.items():
            files[str((proof_dir / relative).relative_to(root))] = info

    for profile, relative in base["manifests"].items():
        raw = read(base_path / relative)
        if {s["name"] for s in raw["samples"]} & set(PROCESSES):
            raise ValueError("a new background already appears in the nominal manifest")
        raw["pooled_backgrounds"] = refs[profile]
        shared.save(root / relative, raw)
    for key, relative in base["cards"].items():
        card = read(base_path / relative)
        card["required_samples"].extend(PROCESSES)
        shared.save(root / relative, card)
    # Operational paths for the unchanged events point at their sealed public
    # base release. The links expose familiar event directories to students.
    for relative in ("nominal","supplemental"):
        link = root / relative
        if not link.exists():
            link.symlink_to(Path("..") / base_path.name / relative, target_is_directory=True)
    reference = root / "reference"
    reference.mkdir(exist_ok=True)
    old = reference / base_path.name
    if not old.exists():
        old.symlink_to(Path("../..") / base_path.name / "reference", target_is_directory=True)
    updated_plots = root.parent.parent / "student-results" / root.name / "fourlepton"
    (reference / "README.md").write_text(
        f"The archived instructor plots in {base_path.name}/ use the original 27 samples. "
        f"The updated instructor plots are in {updated_plots}/. "
        "Run student_events.py report, cuts or xgboost to create results with the three added heavy-flavour backgrounds.\n")
    for path in root.rglob("*"):
        if path.is_file() and not path.is_symlink() and str(path.relative_to(root)) not in files:
            files[str(path.relative_to(root))] = dict(sha256=shared.digest(path),size=path.stat().st_size)
    index = dict(schema_version=1, kind="fourlepton_shared_events", status="complete",
        published_utc=datetime.now(timezone.utc).isoformat(), base_release=dict(name=base_path.name,
            sha256=shared.digest(base_path / "dataset.json")), manifests=base["manifests"], cards=base["cards"],
        archived_reference=str(old.relative_to(root)), supplemental_manifests={},
        samples=[*base["samples"], *added], files=files,
        binary_bytes=base["binary_bytes"] + sum(v["size"] for k,v in files.items() if k.endswith((".root",".lhe.gz",".lhe"))))
    # Seal the inventories before publishing the completion marker.
    for path in root.rglob("*"):
        if not path.is_symlink():
            os.chmod(path, 0o555 if path.is_dir() else 0o444)
    shared.save(root / "dataset.json",index)
    (root / "dataset.json").chmod(0o444)
    root.chmod(0o555)
    print(f"Published {len(index['samples'])} four-lepton samples at {root}", flush=True)


def check_destination(args):
    root = args.destination.resolve()
    base = args.base.resolve()
    if root.parent != base.parent or root == base:
        raise ValueError("choose a new sibling directory of the immutable base release")
    root.mkdir(mode=0o700, parents=True, exist_ok=True)
    if (root / "dataset.json").exists():
        raise ValueError("release already published")
    checkpoint = root / ".reducible_publication.json"
    identity = dict(kind="fourlepton_reducible_release",
        index_sha256=shared.digest(args.index), receipt_sha256=shared.digest(args.receipt),
        base_sha256=shared.digest(base / "dataset.json"), base=str(base))
    if checkpoint.is_file():
        if read(checkpoint) != identity:
            raise ValueError("unfinished publication has different inputs")
    elif any(root.iterdir()):
        raise ValueError("refusing to publish into a nonempty, unmanaged directory")
    else:
        shared.save(checkpoint, identity)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("index", "receipt", "base", "destination"):
        p.add_argument("--"+key, type=Path, required=True)
    p.add_argument("--workers", type=int, default=12)
    p.add_argument("--phase", choices=("copy","audit","seal"), default="copy")
    args = p.parse_args()
    check_destination(args)
    items = jobs(args)
    if args.phase == "seal":
        seal(args, items)
    else:
        with ProcessPoolExecutor(args.workers) as executor:
            list(executor.map(copy_shard if args.phase == "copy" else audit_shard, items))


if __name__ == "__main__":
    main()
