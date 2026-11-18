"""Publish and verify the prepared 40 TeV diphoton analysis inputs."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from hfourlepton.shared_events import digest, safe_path, save

DEFAULT_DATASET = Path("/home/shared/HiggsSSC/diphoton/ssc40-20261118")
NOMINAL = {"signal_gg_h_aa", "bkg_prompt_aa", "bkg_gamma_j", "bkg_dy_ee"}
MANIFEST = "diphoton_samples.json"


def load_dataset(root, *, full=False):
    root = Path(root)
    index = json.loads((root / "dataset.json").read_text())
    if (index.get("schema_version"), index.get("kind"), index.get("status")) != (1, "diphoton_shared_events", "complete"):
        raise ValueError("invalid shared diphoton dataset")
    files = index["files"]
    if not {MANIFEST, *index["cards"].values()} <= files.keys():
        raise ValueError("unsealed diphoton manifest or cards")
    for relative, record in files.items():
        path = safe_path(root, relative)
        if not path.is_file() or not os.access(path, os.R_OK) or path.stat().st_size != record["size"]:
            raise ValueError(f"missing or unreadable diphoton input: {path}")
        if (full or path.suffix != ".root") and digest(path) != record["sha256"]:
            raise ValueError(f"changed diphoton input: {path}")
    return index


def load_prepared(root, run_tag):
    """Return sealed nominal records, or None for a normal generator campaign."""
    root = Path(root)
    source = root / MANIFEST
    if not source.is_file():
        if (root / "dataset.json").exists():
            raise ValueError(f"shared diphoton manifest missing: {source}")
        return None
    index = load_dataset(root, full=True)
    manifest = json.loads(source.read_text())
    if (manifest["energy_tev"] != 40 or manifest["run_tag"] != run_tag
            or index["run_tag"] != run_tag or set(manifest["samples"]) != NOMINAL):
        raise ValueError(f"unexpected prepared diphoton sample set, collision energy or run tag: {source}")
    for name, record in manifest["samples"].items():
        for key in ("var_file", "dat_file", "top_file"):
            if record[key] not in index["files"]:
                raise ValueError(f"unsealed {key} for {name}")
        if record["category"] != ("Signal" if name == "signal_gg_h_aa" else "Backgrounds"):
            raise ValueError(f"invalid process category for {name}")
    return manifest["samples"]


def publish(source, destination, tag, reference):
    """Audit the original production and publish only the analysis-level inputs."""
    from analyze_lo_varfiles import discover_samples
    from hgammagamma.make_gammagamma_report import parse_top_file
    import numpy as np
    import uproot

    source, destination = Path(source).resolve(), Path(destination).resolve()
    if destination.exists():
        raise ValueError(f"choose an unused release directory: {destination}")
    if destination.is_relative_to(source) or source.is_relative_to(destination):
        raise ValueError("source and release must be separate directories")
    reference = json.loads(Path(reference).read_text())["source_provenance"]["root_files"]
    if set(reference) != NOMINAL:
        raise ValueError("expected the four ROOT hashes from the original classifier replay")
    # This uses the full private generator/shower audit during publication.
    nominal = discover_samples(source, tag, NOMINAL)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".publishing-") as temporary:
        root = Path(temporary) / "release"
        root.mkdir()
        records = {}
        for sample in nominal:
            name = sample.name
            source_dir = sample.sample_dir
            raw = json.loads((source_dir / "campaign.json").read_text())
            if (raw["ebeam_gev"] != 20000 or sample.detector_response != "ssc"
                    or not sample.requires_full_sample or sample.shower_quality is None):
                raise ValueError(f"expected audited 40 TeV SSC production for {name}")
            folder = Path(sample.category) / "events" / name
            (root / folder).mkdir(parents=True)
            paths = {}
            for field, path in (("var_file", sample.var_file), ("dat_file", sample.dat_file),
                                ("top_file", sample.dat_file.with_suffix(".top"))):
                relative = folder / path.name
                shutil.copyfile(path, root / relative)
                if digest(root / relative) != digest(path):
                    raise ValueError(f"source changed during publication: {path}")
                paths[field] = str(relative)
            tree_path = root / paths["var_file"]
            if digest(tree_path) != reference[name]["sha256"]:
                raise ValueError(f"ROOT input differs from the classifier replay: {name}")
            with uproot.open(tree_path) as f:
                tree = f["Data2"]
                weights = tree["eventweight"].array(library="np").reshape(-1)
                if (len(weights) != sample.tree_entries or not np.isfinite(weights).all() or
                    not math.isclose(math.fsum(weights), sample.sum_tree_weight,
                                     abs_tol=1.e-8 * sample.sum_abs_weight)):
                    raise ValueError(f"ROOT response population differs from completed analysis: {name}")
            if not parse_top_file(root / paths["top_file"]):
                raise ValueError(f"no histogram distributions for {name}")
            provenance = root / "provenance" / name
            provenance.mkdir(parents=True)
            for path in (source_dir / "campaign.json", source_dir / f"normalization-{tag}.json"):
                shutil.copyfile(path, provenance / path.name)
            fields = ("category", "cross_section_pb", "cross_section_error_pb", "weight_scale",
                      "normalization_sum_weight", "normalization_kind", "ihixs_record_sha256", "shower_quality",
                      "tree_entries", "sum_tree_weight", "sum_abs_weight", "events_read", "sum_weight",
                      "analysis_name", "detector_response", "response_mode")
            records[name] = dict({field: getattr(sample, field) for field in fields}, **paths)
        save(root / MANIFEST, dict(run_tag=tag, energy_tev=40, samples=records))
        cards = {}
        for mode, filename in (("cuts", "ho_baseline_cuts.yaml"), ("xgboost", "ho_xgboost_baseline.yaml")):
            relative = f"cards/{mode}_ssc.yaml"
            (root / relative).parent.mkdir(exist_ok=True)
            shutil.copyfile(REPO / "hgammagamma" / "analysis_cards" / filename, root / relative)
            cards[mode + "_ssc"] = relative
        files = {str(p.relative_to(root)): dict(size=p.stat().st_size, sha256=digest(p))
                 for p in root.rglob("*") if p.is_file()}
        save(root / "dataset.json", dict(schema_version=1, kind="diphoton_shared_events", status="complete",
             run_tag=tag, samples=sorted(records), cards=cards, files=files))
        load_prepared(root, tag)
        # Seal metadata and data before a single same-filesystem rename.
        for path in root.rglob("*"):
            path.chmod(0o555 if path.is_dir() else 0o444)
        root.rename(destination)
        destination.chmod(0o555)
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--classifier-bundle-manifest", type=Path, required=True)
    args = parser.parse_args()
    print(publish(args.source, args.destination, args.run_tag, args.classifier_bundle_manifest))
