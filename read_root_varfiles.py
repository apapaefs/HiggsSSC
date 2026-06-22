"""Read gamma-gamma variable ROOT files produced by the LO analysis."""

from __future__ import annotations

import math
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


def _open_tree(filename: str | Path):
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"ROOT variable file does not exist: {path}")

    root_file = ROOT.TFile.Open(str(path))
    if not root_file or root_file.IsZombie():
        raise OSError(f"Failed to open ROOT variable file: {path}")

    tree = root_file.Get(TREE_NAME)
    if not tree:
        root_file.Close()
        raise KeyError(f"{path} does not contain a {TREE_NAME} tree")
    if not tree.GetBranch("variables"):
        root_file.Close()
        raise KeyError(f"{path}: {TREE_NAME} tree does not contain a variables branch")
    if not tree.GetBranch("eventweight"):
        root_file.Close()
        raise KeyError(f"{path}: {TREE_NAME} tree does not contain an eventweight branch")
    return path, root_file, tree


def read_ROOT_varfile(filename, sample_id, xsec=1.0, max_events=None, include_weight_feature=False):
    """Return feature rows, labels, and weighted event weights from a gamma-gamma tree."""

    _, root_file, tree = _open_tree(filename)
    try:
        n_entries = int(tree.GetEntries())
        if max_events is not None:
            n_entries = min(n_entries, int(max_events))

        features = []
        labels = []
        weights = []

        for entry in range(n_entries):
            tree.GetEntry(entry)
            values = [float(tree.variables[index]) for index in range(VARIABLE_COUNT)]
            weight = _as_scalar(tree.eventweight)
            if not math.isfinite(weight) or not _finite(values):
                continue
            row = [weight, *values] if include_weight_feature else values
            features.append(row)
            labels.append(sample_id)
            weights.append(weight * float(xsec))

        return features, labels, weights
    finally:
        root_file.Close()


def read_named_ROOT_varfile(filename, max_events=None):
    """Return named gamma-gamma rows and event weights from a ROOT variable file."""

    features, _, weights = read_ROOT_varfile(filename, sample_id=0, xsec=1.0, max_events=max_events)
    rows = [dict(zip(FEATURE_NAMES, row)) for row in features]
    return rows, weights
