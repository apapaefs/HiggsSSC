"""Read gamma-gamma variable ROOT files produced by the shared LO/HO analysis."""

from __future__ import annotations

import math
import operator
from pathlib import Path
from typing import Iterable

try:
    import ROOT
except ImportError as exc:  # pragma: no cover - depends on the local HEP stack.
    raise RuntimeError("PyROOT is required to read gamma-gamma variable ROOT files") from exc

ROOT.gROOT.SetBatch(True)

FEATURE_NAMES = [
    "m_gg",
    "pt_gamma1",
    "eta_gamma1",
    "pt_gamma2",
    "eta_gamma2",
    "deltaR_gg",
    "deltaPhi_gg",
    "pt_gg",
    "y_gg",
    "n_selected_photons",
]

VARIABLE_COUNT = len(FEATURE_NAMES)
TREE_NAME = "Data2"


def _as_scalar(value) -> float:
    try:
        return float(value[0])
    except TypeError:
        return float(value)


def _finite(values: Iterable[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def _open_tree(filename: str | Path, strict: bool = False):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"ROOT variable file does not exist: {path}")

    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        if root_file:
            root_file.Close()
        raise OSError(f"Failed to open ROOT variable file: {path}")
    if strict and root_file.TestBit(ROOT.TFile.kRecovered):
        root_file.Close()
        raise OSError(f"Recovered ROOT variable file is not accepted for HO yields: {path}")

    tree = root_file.Get(TREE_NAME)
    if not tree:
        root_file.Close()
        raise KeyError(f"{path} does not contain a {TREE_NAME} tree")
    if strict and not tree.InheritsFrom("TTree"):
        root_file.Close()
        raise TypeError(f"{path}: {TREE_NAME} is not a TTree")
    if not tree.GetBranch("variables"):
        root_file.Close()
        raise KeyError(f"{path}: {TREE_NAME} tree does not contain a variables branch")
    if not tree.GetBranch("eventweight"):
        root_file.Close()
        raise KeyError(f"{path}: {TREE_NAME} tree does not contain an eventweight branch")
    if strict:
        for name, length in (("variables", VARIABLE_COUNT), ("eventweight", 1)):
            leaf = tree.GetLeaf(name)
            if not leaf or leaf.GetLeafCount() or int(leaf.GetLenStatic()) != length:
                root_file.Close()
                raise ValueError(f"{path}: invalid {name} branch length for HO yields")
    return path, root_file, tree


def read_ROOT_varfile(
    filename,
    sample_id,
    xsec=1.0,
    max_events=None,
    include_weight_feature=False,
    selected_only=False,
    strict=False,
):
    """Return feature rows, labels, and weighted event weights from a gamma-gamma tree."""

    if max_events is not None and int(max_events) <= 0:
        raise ValueError("max_events must be positive")

    before = Path(filename).stat() if strict else None
    path, root_file, tree = _open_tree(filename, strict=strict)
    try:
        n_entries = int(tree.GetEntries())
        if max_events is not None and not selected_only:
            n_entries = min(n_entries, int(max_events))

        features = []
        labels = []
        weights = []

        for entry in range(n_entries):
            entry_bytes = tree.GetEntry(entry)
            if strict and entry_bytes <= 0:
                raise OSError(f"Unreadable HO response tree entry {entry}: {path}")
            if strict:
                for name in ("variables", "eventweight"):
                    if tree.GetBranch(name).GetEntry(entry) <= 0:
                        raise OSError(f"Unreadable HO {name} branch entry {entry}: {path}")
            values = [float(tree.variables[index]) for index in range(VARIABLE_COUNT)]
            weight = _as_scalar(tree.eventweight)
            if not math.isfinite(weight) or not _finite(values):
                if strict:
                    raise ValueError(f"Non-finite HO response tree entry {entry}: {path}")
                continue
            if selected_only and values[9] < 2.0:
                continue
            row = [weight, *values] if include_weight_feature else values
            features.append(row)
            labels.append(sample_id)
            weights.append(weight * float(xsec))
            if selected_only and max_events is not None and len(features) >= int(max_events):
                break

        if strict:
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise OSError(f"ROOT variable file changed during inspection: {path}")
        return features, labels, weights
    finally:
        root_file.Close()


def read_named_ROOT_varfile(filename, max_events=None, strict=False):
    """Return named gamma-gamma rows and event weights from a ROOT variable file."""

    features, _, weights = read_ROOT_varfile(filename, sample_id=0, xsec=1.0,
                                          max_events=max_events, strict=strict)
    rows = [dict(zip(FEATURE_NAMES, row)) for row in features]
    return rows, weights


def read_ho_ROOT_varfile(filename):
    """Return every named row, signed weight, source-event ID and tree entry.

    Source-event IDs identify correlated detector-response hypotheses within
    one input file; callers must pair them with the sample identity when
    constructing groups across files. Metadata is kept separate from features.
    """

    before = Path(filename).stat()
    path, root_file, tree = _open_tree(filename, strict=True)
    try:
        source_branch = tree.GetBranch("sourceevent")
        if not source_branch:
            raise KeyError(f"{path}: {TREE_NAME} tree does not contain a sourceevent branch")
        source_leaf = tree.GetLeaf("sourceevent")
        if (not source_leaf or source_leaf.GetLeafCount()
                or int(source_leaf.GetLenStatic()) != 1
                or source_leaf.GetTypeName() != "Long64_t"):
            raise ValueError(f"{path}: invalid sourceevent branch; HO ML requires sourceevent[1]/L")

        branches = [(name, tree.GetBranch(name))
                    for name in ("variables", "eventweight", "sourceevent")]
        rows, weights, source_events, entries = [], [], [], []
        for entry in range(int(tree.GetEntries())):
            if tree.GetEntry(entry) <= 0:
                raise OSError(f"Unreadable HO response tree entry {entry}: {path}")
            for name, branch in branches:
                if branch.GetEntry(entry) <= 0:
                    raise OSError(f"Unreadable HO {name} branch entry {entry}: {path}")
            values = [float(tree.variables[index]) for index in range(VARIABLE_COUNT)]
            weight = _as_scalar(tree.eventweight)
            if not math.isfinite(weight) or not _finite(values):
                raise ValueError(f"Non-finite HO response tree entry {entry}: {path}")

            source_value = tree.sourceevent
            try:
                source_value = source_value[0]
            except TypeError:
                pass
            try:
                source_event = operator.index(source_value)
            except TypeError as exc:
                raise ValueError(f"Invalid HO sourceevent at entry {entry}: {path}") from exc
            if isinstance(source_value, bool) or not 0 <= source_event <= 2**63 - 1:
                raise ValueError(f"Invalid HO sourceevent at entry {entry}: {path}")

            rows.append(dict(zip(FEATURE_NAMES, values)))
            weights.append(weight)
            source_events.append(source_event)
            entries.append(entry)

        after = path.stat()
        stat_fields = ("st_dev", "st_ino", "st_size", "st_mtime_ns")
        if any(getattr(before, field) != getattr(after, field) for field in stat_fields):
            raise OSError(f"ROOT variable file changed during inspection: {path}")
        return rows, weights, source_events, entries
    finally:
        root_file.Close()


def sum_ROOT_varfile_weights(filename) -> float:
    """Return the finite weight sum over every detector-response tree row."""

    _, root_file, tree = _open_tree(filename)
    try:
        total = 0.0
        for entry in range(int(tree.GetEntries())):
            tree.GetEntry(entry)
            weight = _as_scalar(tree.eventweight)
            if math.isfinite(weight):
                total += weight
        return total
    finally:
        root_file.Close()
