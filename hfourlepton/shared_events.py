"""Publish completed four-lepton events for analysis without a private home.

Binary events are copied byte for byte (reflinks on Linux when supported).
Only paths in operational JSON, input lists and shower evidence are relocated.
The original manifests and evidence remain available alongside new checksums;
the ordinary analysis audits still apply, with no special trust bypass.
"""
from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hgammagamma import ho_shower_completion as completion

DEFAULT_DATASET = Path("/home/shared/HiggsSSC/fourlepton/ssc40-20261112")


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def safe_path(root, relative):
    root = Path(root).resolve()
    path = root / relative
    if Path(relative).is_absolute() or not path.resolve().is_relative_to(root):
        raise ValueError(f"path escapes dataset/workspace: {relative}")
    return path


def name(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) or value in (".", ".."):
        raise ValueError(f"invalid name: {value}")
    return value


class Publisher:
    def __init__(self, destination, identity):
        self.root = Path(destination).resolve()
        if (self.root / "dataset.json").exists():
            raise ValueError("dataset is already published; choose a new version directory")
        self.root.mkdir(mode=0o700, parents=True, exist_ok=True)
        self.lock = (self.root / ".publication.lock").open("a")
        try:
            fcntl.flock(self.lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            self.lock.close()
            raise ValueError("another publisher is using this release directory") from error
        try:
            checkpoint = self.root / ".publication.json"
            if checkpoint.exists():
                old = json.loads(checkpoint.read_text())
                if old["identity"] != identity:
                    raise ValueError("unfinished publication has different inputs; choose a new directory")
            elif any(p.name != ".publication.lock" for p in self.root.iterdir()):
                raise ValueError("publication destination is nonempty and is not a resumable publication")
            save(checkpoint, {"identity": identity})
        except BaseException:
            self.lock.close()
            raise
        self.identity = identity
        self.files = {}
        self.copied = {}
        self.showers = {}

    def __enter__(self):
        return self

    def __exit__(self, *error):
        self.lock.close()

    def register(self, path, checksum, **extra):
        relative = str(path.relative_to(self.root))
        self.files[relative] = dict(sha256=checksum, size=path.stat().st_size, **extra)
        return path

    def copy_file(self, source, relative, expected=None):
        source = Path(source).resolve()
        if not source.is_file():
            raise ValueError(f"missing completed input: {source}")
        if source in self.copied:
            path = self.copied[source]
            if expected and self.files[str(path.relative_to(self.root))]["sha256"] != expected:
                raise ValueError(f"conflicting input checksum: {source}")
            return path
        target = safe_path(self.root, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        before = source.stat()
        if not target.exists():
            temporary = target.with_name(target.name + ".copying")
            if sys.platform.startswith("linux"):
                subprocess.run(["cp", "--reflink=auto", "--preserve=timestamps", "--",
                                str(source), str(temporary)], check=True)
            else:
                shutil.copy2(source, temporary)
            temporary.replace(target)
        checksum = digest(target)
        after = source.stat()
        expected = expected or digest(source)
        if (checksum != expected or before.st_size != after.st_size
                or before.st_mtime_ns != after.st_mtime_ns
                or target.stat().st_size != before.st_size):
            raise ValueError(f"input changed or published copy differs: {source}")
        # Never chmod a source inode: reflinks/copies are separate inodes.
        if os.path.samestat(source.stat(), target.stat()):
            raise ValueError("publication must not use hardlinks to original events")
        self.copied[source] = target
        return self.register(target, checksum, original_path=str(source), original_sha256=expected)

    def write(self, relative, text):
        target = safe_path(self.root, relative)
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if target.read_text() != text:
                raise ValueError(f"unfinished publication product changed: {target}")
        else:
            target.write_text(text)
        return self.register(target, digest(target))

    def json(self, relative, value):
        return self.write(relative, json.dumps(value, indent=2, allow_nan=False) + "\n")

    def archive(self, source, expected=None):
        source = Path(source)
        checksum = expected or digest(source)
        return self.copy_file(source, Path("provenance") / (checksum + "-" + source.name), checksum)

    def artifact(self, path):
        info = path.stat()
        return dict(path=str(path), size=info.st_size, mtime_ns=info.st_mtime_ns,
                    sha256=self.files[str(path.relative_to(self.root))]["sha256"])

    def shower(self, original, namespace):
        fingerprint = original["fingerprint"]
        if fingerprint in self.showers:
            return copy.deepcopy(self.showers[fingerprint])
        sample = name(original["sample"])
        identity = dict(sample=sample, nevents_requested=original["requested_events"],
                        lhe_file=original["reader_lhe_path"], lhe={"events": original["lhe_events"]})
        completion.validate_completion(original, identity)
        folder = Path(namespace) / "events" / sample
        herwig = self.root / folder / "herwig"
        translated = {original["herwig_dir"]: str(herwig)}
        record = copy.deepcopy(original)
        source = record["source_weights"]["artifact"]
        lhe = self.copy_file(source["path"], folder / Path(source["path"]).name, source["sha256"])
        translated[original["reader_lhe_path"]] = str(lhe)
        record["source_weights"]["artifact"] = self.artifact(lhe)
        record["source_weights"]["fingerprint"] = completion.completion_fingerprint(record["source_weights"])
        record["reader_lhe_path"] = str(lhe)
        record["herwig_dir"] = str(herwig)
        roots = []
        for entry, metadata in zip(original["artifacts"]["roots"], record["root_metadata"]):
            path = self.copy_file(entry["path"], folder / "herwig/events" / Path(entry["path"]).name,
                                  entry["sha256"])
            translated[entry["path"]] = str(path)
            artifact = self.artifact(path)
            roots.append(artifact)
            metadata.update({key: artifact[key] for key in ("path", "size", "mtime_ns")})
        record["artifacts"]["roots"] = roots
        for key, entry in original["artifacts"].items():
            if key == "roots":
                continue
            archived = self.archive(entry["path"], entry["sha256"])
            text = archived.read_text()
            for old, new in sorted(translated.items(), key=lambda pair: -len(pair[0])):
                text = text.replace(old, new)
            path = self.write(folder / "herwig" / Path(entry["path"]).name, text)
            record["artifacts"][key] = self.artifact(path)
        record["fingerprint"] = completion.completion_fingerprint(record)
        self.json(folder / "herwig/matched_completion.json", record)
        self.showers[fingerprint] = record
        return copy.deepcopy(record)

    def manifest(self, source, relative, namespace="nominal"):
        source = Path(source).resolve()
        raw = json.loads(source.read_text())
        if raw.get("status") != "complete" or raw.get("order") != "ho":
            raise ValueError(f"publication requires a completed HO manifest: {source}")
        profile = raw["detector_profile"]
        if profile not in ("perfect", "ssc"):
            raise ValueError("unsupported detector response")
        self.archive(source)
        names = [name(sample["name"]) for sample in raw["samples"]]
        if not names or len(set(names)) != len(names):
            raise ValueError("empty or duplicate sample inventory")
        for sample in raw["samples"]:
            print(f"Publishing {namespace}/{profile}/{sample['name']}", flush=True)
            if sample["status"] != "complete":
                raise ValueError(f"incomplete sample: {sample['name']}")
            original = copy.deepcopy(sample)
            if (sample["input_lhe"] != sample["shower_completion"]["reader_lhe_path"]
                    or sample["input_lhe_sha256"] != sample["shower_completion"]["source_weights"]["artifact"]["sha256"]):
                raise ValueError("original sample input differs from its shower weight proof")
            folder = Path(namespace) / "analysis" / profile / sample["name"]
            sample["shower_completion"] = self.shower(sample["shower_completion"], namespace)
            completed = sample["shower_completion"]
            sample["input_lhe"] = completed["reader_lhe_path"]
            sample["output_root"] = sample["root_file"] = str(self.copy_file(
                sample["output_root"], folder / Path(sample["output_root"]).name,
                sample["output_root_sha256"]))
            summary_source = self.archive(sample["analysis_summary"], sample["analysis_summary_sha256"])
            summary = json.loads(summary_source.read_text())
            if summary["output_root"] != original["output_root"]:
                raise ValueError("original analyzer summary points to a different ROOT file")
            input_original = Path(summary["input"])
            listed = [str(Path(line.strip()).resolve()) for line in input_original.read_text().splitlines()
                      if line.strip() and not line.lstrip().startswith("#")]
            if listed != [entry["path"] for entry in original["shower_completion"]["artifacts"]["roots"]]:
                raise ValueError("original analyzer input list differs from completed shower")
            self.archive(input_original)
            summary["input"] = str(self.write(folder / "showers.input", "".join(
                entry["path"] + "\n" for entry in completed["artifacts"]["roots"])))
            summary["output_root"] = sample["output_root"]
            optional = sample.get("optional_weight_names_file")
            if optional:
                if summary.get("optional_weight_names_file") != optional["path"]:
                    raise ValueError("original optional-weight input differs from manifest")
                optional["path"] = str(self.copy_file(optional["path"], folder / "lhe_optional_weight_names.txt",
                                                      optional["sha256"]))
                summary["optional_weight_names_file"] = optional["path"]
            new_summary = self.json(folder / Path(sample["analysis_summary"]).name, summary)
            sample["analysis_summary"] = str(new_summary)
            sample["analysis_summary_sha256"] = digest(new_summary)
            sample["source_lhe"] = str(self.copy_file(sample["source_lhe"],
                Path(namespace) / "events" / sample["name"] / ("hard-" + Path(sample["source_lhe"]).name),
                sample["source_lhe_sha256"]))
            if sample.get("hard_event_slice"):
                proof = sample["hard_event_slice"]
                proof["artifact"]["path"] = sample["source_lhe"]
                proof["identity"]["source"] = str(self.copy_file(proof["identity"]["source"],
                    Path(namespace) / "production" / name(sample["production_group_id"]) /
                    Path(proof["identity"]["source"]).name, proof["identity"]["source_sha256"]))
            # Historical source manifests and normalization records keep their
            # original bytes/fingerprints; analysis never executes their paths.
            if sample.get("source_manifest"):
                provenance = sample["source_manifest"]
                provenance["path"] = str(self.archive(provenance["path"], provenance["sha256"]))
            sample["publication"] = dict(original_manifest_sha256=digest(source),
                original_shower_fingerprint=original["shower_completion"]["fingerprint"],
                binary_events_unchanged=True, original_analyzer_cli=sample.pop("analyzer_cli", None))
        raw["shared_dataset"] = dict(schema_version=1, role=namespace, original_manifest_sha256=digest(source))
        return self.json(relative, raw)

    def reference(self, source, relative):
        source = Path(source)
        if not source.is_dir():
            raise ValueError(f"missing reference directory: {source}")
        for path in sorted(source.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"reference directory contains a symlink: {path}")
            if path.is_file():
                self.copy_file(path, Path(relative) / path.relative_to(source))

    def finish(self, manifests, cards, extra):
        from hfourlepton.reprocess import audit_completed_manifest
        for relative in [*manifests.values(), *extra.values()]:
            print(f"Auditing relocated {relative}", flush=True)
            audit_completed_manifest(self.root / relative)
        samples = json.loads((self.root / manifests["ssc"]).read_text())["samples"]
        dataset = dict(schema_version=1, kind="fourlepton_shared_events", status="complete",
            published_utc=datetime.now(timezone.utc).isoformat(), manifests=manifests, cards=cards,
            publisher_sha256=digest(__file__),
            supplemental_manifests=extra, samples=[dict(name=s["name"], requested=s["nevents_requested"],
                saved=s["shower_completion"]["saved_events"], category=s["category"], channel=s["channel"])
                for s in samples], source=self.identity, files=self.files,
            binary_bytes=sum(entry["size"] for key, entry in self.files.items()
                             if key.endswith((".root", ".lhe", ".lhe.gz"))))
        save(self.root / "dataset.json", dataset)
        # The release is readable by every Timur account and writable by none
        # through normal file permissions, including the publishing owner.
        for path in self.root.rglob("*"):
            path.chmod(0o555 if path.is_dir() else 0o444)
        self.root.chmod(0o555)
        print(f"Published {len(samples)} samples: {self.root}", flush=True)


def load_dataset(path, *, full=False):
    path = Path(path).resolve()
    index = json.loads((path / "dataset.json").read_text())
    if (index.get("schema_version") != 1 or index.get("kind") != "fourlepton_shared_events"
            or index.get("status") != "complete"):
        raise ValueError("shared dataset is incomplete or unsupported")
    for relative, entry in index["files"].items():
        file = safe_path(path, relative)
        if not file.is_file() or file.stat().st_size != entry["size"] or not os.access(file, os.R_OK):
            raise ValueError(f"missing, unreadable or changed shared file: {file}")
        # Metadata is always checked; full verification also hashes large events.
        if (full or not relative.endswith((".root", ".lhe", ".lhe.gz"))) and digest(file) != entry["sha256"]:
            raise ValueError(f"shared checksum changed: {file}")
    for relative in [*index["manifests"].values(), *index["cards"].values(),
                     *index.get("supplemental_manifests", {}).values()]:
        if relative not in index["files"]:
            raise ValueError(f"unsealed dataset input: {relative}")
    return index


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--perfect", required=True, type=Path)
    parser.add_argument("--ssc", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--reference-dir", required=True, type=Path,
                        help="Completed revision with report/, cuts/ and xgboost/.")
    parser.add_argument("--extra-ssc", type=Path, help="Independent qq test manifest, kept separate from nominal samples.")
    parser.add_argument("--extended-reference", type=Path)
    args = parser.parse_args(argv)
    sources = {key: dict(path=str(value.resolve()), sha256=digest(value)) for key, value in
               (("perfect", args.perfect), ("ssc", args.ssc))}
    if args.extra_ssc:
        sources["extra_ssc"] = dict(path=str(args.extra_ssc.resolve()), sha256=digest(args.extra_ssc))
    sources["reference_dir"] = str(args.reference_dir.resolve())
    sources["extended_reference"] = str(args.extended_reference.resolve()) if args.extended_reference else None
    with Publisher(args.destination, sources) as publisher:
        manifests, cards, extra = {}, {}, {}
        for profile in ("perfect", "ssc"):
            relative = f"manifest_{profile}.json"
            publisher.manifest(getattr(args, profile), relative)
            manifests[profile] = relative
            for mode in ("cuts", "xgboost"):
                card = json.loads((args.reference_dir / mode / "summary.json").read_text())["card"]
                card["detector_response"] = profile
                card_path = f"cards/{mode}_{profile}.json"
                publisher.json(card_path, card)
                cards[f"{mode}_{profile}"] = card_path
        if args.extra_ssc:
            extra["qq_test_ssc"] = "supplemental/manifest_ssc.json"
            publisher.manifest(args.extra_ssc, extra["qq_test_ssc"], "supplemental")
        for mode in ("report", "cuts", "xgboost"):
            publisher.reference(args.reference_dir / mode, f"reference/{mode}")
        if args.extended_reference:
            publisher.reference(args.extended_reference, "reference/extended_test")
        publisher.finish(manifests, cards, extra)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        sys.exit(f"Error: {error}")
