"""Optional XGBoost analysis for gamma-gamma LO ``_var.root`` files."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from read_root_varfiles import FEATURE_NAMES, read_ROOT_varfile


def _require_xgboost_dependencies():
    try:
        import xgboost as xgb
        from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score, roc_curve
        from sklearn.model_selection import train_test_split
        import tqdm  # noqa: F401 - imported to validate the documented optional dependency.

        os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on optional environment.
        raise ImportError("xgboost mode requires xgboost, scikit-learn, and tqdm") from exc

    return {
        "plt": plt,
        "xgb": xgb,
        "accuracy_score": accuracy_score,
        "confusion_matrix": confusion_matrix,
        "roc_auc_score": roc_auc_score,
        "roc_curve": roc_curve,
        "train_test_split": train_test_split,
    }


def _as_path_list(paths) -> list[Path]:
    if paths is None:
        return []
    if isinstance(paths, (str, Path)):
        return [Path(paths)]
    return [Path(path) for path in paths]


def _expand_per_file(values, files: Sequence[Path], default):
    if values is None:
        return [default for _ in files]
    if isinstance(values, (int, float)):
        return [float(values) for _ in files]
    expanded = [default if value is None else float(value) for value in values]
    if len(expanded) == 1 and len(files) > 1:
        return expanded * len(files)
    if len(expanded) != len(files):
        raise ValueError(f"expected {len(files)} values, got {len(expanded)}")
    return expanded


def _expand_metadata(metadata, files: Sequence[Path]) -> list[dict[str, Any]]:
    if metadata is None:
        return [{} for _ in files]
    rows = [dict(item or {}) for item in metadata]
    if len(rows) == 1 and len(files) > 1:
        return rows * len(files)
    if len(rows) != len(files):
        raise ValueError(f"expected {len(files)} metadata rows, got {len(rows)}")
    return rows


def _normalisation_denominator(normalisation_weight, raw_weights) -> tuple[float, str]:
    if normalisation_weight is not None and normalisation_weight > 0:
        return float(normalisation_weight), "input_weight_sum"
    fallback = float(np.sum(raw_weights))
    if fallback > 0.0:
        return fallback, "loaded_weight_sum"
    return float(len(raw_weights)), "entry_count"


def _balanced_training_weights(labels, raw_weights):
    labels = np.asarray(labels)
    raw_weights = np.asarray(raw_weights, dtype=float)
    base = np.abs(raw_weights)
    base[base == 0.0] = 1.0
    balanced = np.zeros_like(base, dtype=float)
    for label in np.unique(labels):
        mask = labels == label
        class_sum = np.sum(base[mask])
        balanced[mask] = base[mask] * np.sum(mask) / class_sum if class_sum > 0.0 else 1.0
    return balanced


def _load_group(
    files: Sequence[Path],
    label: int,
    xsecs_fb: Sequence[float],
    rate_factors: Sequence[float],
    normalisation_weights: Sequence[float | None],
    luminosity: float,
    max_events=None,
    metadata: Sequence[dict[str, Any]] | None = None,
):
    rows = []
    labels = []
    raw_weights = []
    physical_weights = []
    sources = []
    sample_summaries = []
    metadata = metadata or [{} for _ in files]

    for path, xsec_fb, rate_factor, normalisation_weight, file_metadata in zip(
        files,
        xsecs_fb,
        rate_factors,
        normalisation_weights,
        metadata,
    ):
        features, sample_labels, weights = read_ROOT_varfile(path, label, 1.0, max_events=max_events)
        weights = np.asarray(weights, dtype=float)
        normalisation, normalisation_source = _normalisation_denominator(normalisation_weight, weights)
        effective_xsec_fb = float(xsec_fb) * float(rate_factor)
        physical = float(luminosity) * effective_xsec_fb * weights / normalisation

        rows.extend(features)
        labels.extend(sample_labels)
        raw_weights.extend(weights.tolist())
        physical_weights.extend(physical.tolist())
        sources.extend([str(path)] * len(features))
        sample_summaries.append(
            {
                **dict(file_metadata or {}),
                "input_file": str(path),
                "xsec_fb": float(xsec_fb),
                "rate_factor": float(rate_factor),
                "effective_xsec_fb": effective_xsec_fb,
                "normalisation_weight": normalisation,
                "normalisation_source": normalisation_source,
                "entries_read": len(features),
                "sum_weight": float(np.sum(weights)),
            }
        )

    return rows, labels, raw_weights, physical_weights, sources, sample_summaries


def _best_threshold(scores, labels, physical_weights, systematics=0.0):
    thresholds = np.linspace(0.0, 1.0, 501)
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    physical_weights = np.asarray(physical_weights, dtype=float)
    total_signal = float(np.sum(physical_weights[labels == 1]))
    total_background = float(np.sum(physical_weights[labels == 0]))
    best = {
        "threshold": 0.5,
        "signal_events": 0.0,
        "background_events": 0.0,
        "significance": 0.0,
        "signal_efficiency": 0.0,
        "background_efficiency": 0.0,
    }

    for threshold in thresholds:
        selected = scores >= threshold
        signal = float(np.sum(physical_weights[(labels == 1) & selected]))
        background = float(np.sum(physical_weights[(labels == 0) & selected]))
        if background <= 0.0:
            continue
        denominator = math.sqrt(background + (float(systematics) * background) ** 2)
        significance = signal / denominator if denominator > 0.0 else 0.0
        if significance > best["significance"]:
            best = {
                "threshold": float(threshold),
                "signal_events": signal,
                "background_events": background,
                "significance": float(significance),
                "signal_efficiency": signal / total_signal if total_signal > 0.0 else 0.0,
                "background_efficiency": background / total_background if total_background > 0.0 else 0.0,
            }
    return best


def _write_scores_csv(path, labels, scores, physical_weights, sources):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["label", "score", "physical_weight", "source"])
        for label, score, weight, source in zip(labels, scores, physical_weights, sources):
            writer.writerow([int(label), float(score), float(weight), source])


def _write_roc_plot(path, plt, roc_curve, labels, scores, physical_weights):
    fpr, tpr, _ = roc_curve(labels, scores, sample_weight=physical_weights)
    plt.figure(figsize=(6, 5))
    plt.plot(fpr, tpr, label="XGBoost")
    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.xlabel("Background efficiency")
    plt.ylabel("Signal efficiency")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def _write_feature_importance_plot(path, plt, model):
    importances = np.asarray(model.feature_importances_, dtype=float)
    order = np.argsort(importances)
    plt.figure(figsize=(7, 4.8))
    plt.barh(np.asarray(FEATURE_NAMES)[order], importances[order])
    plt.xlabel("Feature importance")
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def _summarize_full_sample(
    model,
    sample: dict[str, Any],
    label: int,
    threshold: float,
    luminosity: float,
    max_events=None,
) -> dict[str, Any]:
    features, _, weights = read_ROOT_varfile(sample["input_file"], label, 1.0, max_events=max_events)
    weights = np.asarray(weights, dtype=float)
    scores = model.predict_proba(np.asarray(features, dtype=float))[:, 1]
    selected = scores >= threshold
    normalisation = float(sample["normalisation_weight"])
    physical_weights = float(luminosity) * float(sample["effective_xsec_fb"]) * weights / normalisation
    sum_weight = float(np.sum(weights))
    sum_selected_weight = float(np.sum(weights[selected]))
    efficiency = sum_selected_weight / sum_weight if sum_weight > 0.0 else 0.0
    selected_xsec_pb = (float(sample["effective_xsec_fb"]) * efficiency) / 1000.0

    return {
        "sample": sample.get("sample", Path(sample["input_file"]).parent.name),
        "category": sample.get("category", "Signal" if label == 1 else "Backgrounds"),
        "input_file": sample["input_file"],
        "cross_section_pb": float(sample["xsec_fb"]) / 1000.0,
        "weight_scale": float(sample["rate_factor"]),
        "entries_read": len(features),
        "selected_entries": int(np.sum(selected)),
        "sum_weight": sum_weight,
        "sum_selected_weight": sum_selected_weight,
        "analysis_efficiency": efficiency,
        "selected_cross_section_pb": selected_xsec_pb,
        "expected_events": float(np.sum(physical_weights[selected])),
        "mc_events_after_analysis": int(np.sum(selected)),
    }


def run_signal_background_analysis(
    signal_files,
    background_files,
    output_dir="xgboost_results",
    signal_xsecs_fb=None,
    background_xsecs_fb=None,
    signal_rate_factors=None,
    background_rate_factors=None,
    signal_normalisation_weights=None,
    background_normalisation_weights=None,
    signal_metadata=None,
    background_metadata=None,
    luminosity=100.0,
    test_size=0.35,
    seed=12345,
    systematics=0.0,
    max_events=None,
    model_params=None,
):
    """Train and evaluate a binary gamma-gamma signal-vs-background classifier."""

    deps = _require_xgboost_dependencies()
    plt = deps["plt"]
    xgb = deps["xgb"]
    accuracy_score = deps["accuracy_score"]
    confusion_matrix = deps["confusion_matrix"]
    roc_auc_score = deps["roc_auc_score"]
    roc_curve = deps["roc_curve"]
    train_test_split = deps["train_test_split"]

    signal_files = _as_path_list(signal_files)
    background_files = _as_path_list(background_files)
    if not signal_files:
        raise ValueError("at least one signal ROOT variable file is required")
    if not background_files:
        raise ValueError("at least one background ROOT variable file is required")

    signal_xsecs_fb = _expand_per_file(signal_xsecs_fb, signal_files, 1.0)
    background_xsecs_fb = _expand_per_file(background_xsecs_fb, background_files, 1.0)
    signal_rate_factors = _expand_per_file(signal_rate_factors, signal_files, 1.0)
    background_rate_factors = _expand_per_file(background_rate_factors, background_files, 1.0)
    signal_normalisation_weights = _expand_per_file(signal_normalisation_weights, signal_files, None)
    background_normalisation_weights = _expand_per_file(background_normalisation_weights, background_files, None)
    signal_metadata = _expand_metadata(signal_metadata, signal_files)
    background_metadata = _expand_metadata(background_metadata, background_files)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    s_rows, s_labels, s_raw, s_phys, s_sources, signal_summary = _load_group(
        signal_files,
        1,
        signal_xsecs_fb,
        signal_rate_factors,
        signal_normalisation_weights,
        luminosity,
        max_events,
        signal_metadata,
    )
    b_rows, b_labels, b_raw, b_phys, b_sources, background_summary = _load_group(
        background_files,
        0,
        background_xsecs_fb,
        background_rate_factors,
        background_normalisation_weights,
        luminosity,
        max_events,
        background_metadata,
    )

    X = np.asarray(s_rows + b_rows, dtype=float)
    y = np.asarray(s_labels + b_labels, dtype=int)
    raw_weights = np.asarray(s_raw + b_raw, dtype=float)
    physical_weights = np.asarray(s_phys + b_phys, dtype=float)
    sources = np.asarray(s_sources + b_sources)
    if len(np.unique(y)) != 2:
        raise ValueError("training sample must contain both signal and background events")

    training_weights = _balanced_training_weights(y, raw_weights)
    split = train_test_split(
        X,
        y,
        physical_weights,
        training_weights,
        sources,
        test_size=float(test_size),
        random_state=int(seed),
        stratify=y,
    )
    X_train, X_test, y_train, y_test, phys_train, phys_test, train_w, test_w, src_train, src_test = split

    params = {
        "objective": "binary:logistic",
        "eval_metric": "logloss",
        "n_estimators": 300,
        "max_depth": 3,
        "learning_rate": 0.05,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "random_state": int(seed),
        "n_jobs": 1,
    }
    if model_params:
        params.update(model_params)

    model = xgb.XGBClassifier(**params)
    model.fit(X_train, y_train, sample_weight=train_w)

    scores = model.predict_proba(X_test)[:, 1]
    predictions = (scores >= 0.5).astype(int)
    best = _best_threshold(scores, y_test, phys_test, systematics)
    best_predictions = (scores >= best["threshold"]).astype(int)

    model_file = output_dir / "signal_background_xgboost.json"
    metrics_file = output_dir / "metrics.json"
    scores_file = output_dir / "scores.csv"
    roc_file = output_dir / "roc.png"
    feature_file = output_dir / "feature_importance.png"

    model.save_model(str(model_file))
    _write_scores_csv(scores_file, y_test, scores, phys_test, src_test)
    _write_roc_plot(roc_file, plt, roc_curve, y_test, scores, phys_test)
    _write_feature_importance_plot(feature_file, plt, model)

    summary_rows = [
        *[_summarize_full_sample(model, sample, 1, best["threshold"], luminosity, max_events) for sample in signal_summary],
        *[_summarize_full_sample(model, sample, 0, best["threshold"], luminosity, max_events) for sample in background_summary],
    ]
    total_signal = float(sum(row["expected_events"] for row in summary_rows if row["category"] == "Signal"))
    total_background = float(sum(row["expected_events"] for row in summary_rows if row["category"] != "Signal"))
    denominator = math.sqrt(total_background + (float(systematics) * total_background) ** 2) if total_background > 0.0 else 0.0

    metadata = {
        "n_events": int(len(y)),
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "test_size": float(test_size),
        "seed": int(seed),
        "luminosity_fb_inverse": float(luminosity),
        "systematics": float(systematics),
        "accuracy_threshold_0p5": float(accuracy_score(y_test, predictions)),
        "auc_unweighted": float(roc_auc_score(y_test, scores)),
        "auc_weighted": float(roc_auc_score(y_test, scores, sample_weight=phys_test)),
        "best_threshold": best["threshold"],
        "best_threshold_test_metrics": best,
        "confusion_matrix_at_best_threshold": confusion_matrix(y_test, best_predictions, labels=[0, 1]).tolist(),
        "expected_selected_signal_events": total_signal,
        "expected_selected_background_events": total_background,
        "significance_full_samples": total_signal / denominator if denominator > 0.0 else 0.0,
        "feature_names": FEATURE_NAMES,
        "outputs": {
            "model": str(model_file),
            "metrics": str(metrics_file),
            "scores": str(scores_file),
            "roc": str(roc_file),
            "feature_importance": str(feature_file),
        },
    }

    metrics_file.write_text(json.dumps(metadata, indent=2))
    print("XGBoost gamma-gamma analysis complete")
    print("AUC (weighted) =", metadata["auc_weighted"])
    print("Best threshold =", metadata["best_threshold"])
    print("Expected selected S, B =", total_signal, total_background)
    print("Wrote outputs to", output_dir)
    return {"model": model, "metadata": metadata, "summary_rows": summary_rows}
