"""Frozen pair scoring for the SSC double-fake calculation."""
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from hgammagamma.frozen_gammagamma_classifier import (
    FEATURE_NAMES, load_frozen_classifier, validate_metadata,
)


def load_model(directory):
    if directory is None:
        return None
    directory = Path(directory)
    model_path = directory / "signal_background_xgboost.json"
    metadata_path = directory / "model_metadata.json"
    if (directory / "bundle.json").is_file():
        classifier, metadata, _, replay = load_frozen_classifier(directory)
    else:
        import xgboost as xgb
        replay = None
        metadata = json.loads(metadata_path.read_text())
        validate_metadata(metadata)
        classifier = xgb.XGBClassifier(n_jobs=1)
        classifier.load_model(model_path)
        classifier.set_params(n_jobs=1)
    threshold = float(metadata["best_threshold"])
    if classifier.n_features_in_ != len(FEATURE_NAMES):
        raise ValueError("feature count differs from the original classifier")
    provenance = dict(threshold=threshold, features=FEATURE_NAMES, replay=replay,
                      model_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
                      metadata_sha256=hashlib.sha256(metadata_path.read_bytes()).hexdigest())
    return classifier, provenance


def moments(pairs, mask, n, scale, scenario):
    selected = pairs[mask]
    a = np.minimum(1., scale * selected["f1"])
    b = np.minimum(1., scale * selected["f2"])
    q = np.minimum(1., scale * selected["fq"])
    g = np.minimum(1., scale * selected["fg"])
    weight = {"two_eligible_only": a*b*(selected["eligible_jets"] == 2),
              "independent_pair_sum": a*b, "all_quark_independent": q*q,
              "all_quark_any_correlation": q,
              "all_gluon_two_eligible": g*g*(selected["eligible_jets"] == 2)}[scenario]
    by_event = np.minimum(1., np.bincount(selected["source"].astype(int), weights=weight, minlength=n))
    if len(by_event) != n:
        raise ValueError("pair source ID outside the ROOT population")
    return float(by_event.sum()), float(by_event @ by_event), int(np.count_nonzero(by_event))


def response_rows(analyzer, root, seed, model=None):
    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp)/"moments.csv"
        features = Path(tmp)/"pairs.csv"
        args = [analyzer, "--dijet-bound", root, output, str(seed)]
        subprocess.run(args + ([features] if model else []), check=True)
        with output.open() as f:
            rows = list(csv.DictReader(f))
        if not model:
            return rows
        pairs = np.genfromtxt(features, names=True, delimiter=",", ndmin=1)
    classifier, metadata = model
    window = (pairs["m_gg"] >= 120.) & (pairs["m_gg"] <= 130.)
    matrix = np.column_stack([pairs[key] for key in FEATURE_NAMES]).astype(np.float32)
    if not np.isfinite(matrix).all():
        raise ValueError("nonfinite pair features")
    score = classifier.predict_proba(matrix)[:,1] if len(pairs) else np.array([])
    if not np.isfinite(score).all():
        raise ValueError("nonfinite frozen classifier score")
    scored = []
    for row in rows:
        params = (int(row["events_read"]), float(row["fake_scale"]), row["scenario"])
        check = moments(pairs, window, *params)
        for value, key in zip(check, ("sum_probability","sum_probability_squared","contributing_events")):
            if not math.isclose(value, float(row[key]), rel_tol=1.e-10, abs_tol=1.e-15):
                raise ValueError("pair features do not reproduce the mass-window calculation")
        sw, sw2, count = moments(pairs, score >= metadata["threshold"], *params)
        scored.append(dict(row, scenario="frozen_"+row["scenario"], sum_probability=sw,
                           sum_probability_squared=sw2, contributing_events=count))
    return rows + scored
