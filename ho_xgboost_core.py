"""NumPy-only splitting and signed-weight statistics for the HO classifier.

Detector-response rows from the same source event are correlated.  Splits and
Monte Carlo variances therefore use ``(sample_id, source_event)`` groups.  Pass
the complete response trees to :func:`grouped_sample_split`, before applying
the two-photon eligibility selection.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Any

import numpy as np


PARTITIONS = ("train", "validation", "test")


@dataclass(frozen=True)
class GroupedSplit:
    """Row assignments and actual source-event sampling fractions per sample."""

    partitions: np.ndarray
    sample_fractions: dict[str, dict[str, float]]
    sample_group_counts: dict[str, dict[str, int]]


def _vector(values, name: str, *, dtype=None) -> np.ndarray:
    result = np.asarray(values, dtype=dtype)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return result


def _weights(values) -> np.ndarray:
    result = _vector(values, "weights", dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError("weights must be finite")
    return result


def _identities(sample_ids, source_events) -> tuple[np.ndarray, np.ndarray]:
    samples = _vector(sample_ids, "sample_ids")
    sources = _vector(source_events, "source_events")
    if len(samples) != len(sources):
        raise ValueError("sample_ids and source_events must have equal lengths")
    if any(not isinstance(value, (str, np.str_)) or not value for value in samples):
        raise ValueError("sample_ids must be nonempty strings")
    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
           or value < 0 for value in sources):
        raise ValueError("source_events must be non-negative integers")
    return samples.astype(str), sources


def _group_indices(samples: np.ndarray, sources: np.ndarray) -> np.ndarray:
    groups = {}
    indices = np.empty(len(samples), dtype=np.int64)
    for index, (sample, source) in enumerate(zip(samples, sources)):
        key = (str(sample), int(source))
        if key not in groups:
            groups[key] = len(groups)
        indices[index] = groups[key]
    return indices


def grouped_sample_split(
    sample_ids,
    source_events,
    *,
    fractions=(0.6, 0.2, 0.2),
    seed=12345,
) -> GroupedSplit:
    """Split complete source-event groups separately within every sample.

    Largest-remainder allocation approximates the requested fractions while
    keeping at least one source event in each partition.  Every sample must
    therefore contain at least three source events.  Row reordering or adding
    another sample does not change an existing source event's assignment.
    """

    samples, sources = _identities(sample_ids, source_events)
    requested = _vector(fractions, "fractions", dtype=float)
    if len(requested) != 3 or not np.all(np.isfinite(requested)) or np.any(requested <= 0):
        raise ValueError("fractions must contain three finite, positive values")
    if not math.isclose(float(requested.sum()), 1.0, rel_tol=0.0, abs_tol=1.e-12):
        raise ValueError("split fractions must sum to one")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if not len(samples):
        raise ValueError("cannot split an empty sample")

    assignments = np.empty(len(samples), dtype="<U10")
    sample_fractions = {}
    sample_counts = {}
    for sample in sorted(set(samples)):
        mask = samples == sample
        unique_sources = np.unique(sources[mask])
        size = len(unique_sources)
        if size < 3:
            raise ValueError(f"sample {sample!r} needs at least three source events for train/validation/test")
        expected = size * requested
        counts = np.floor(expected).astype(int)
        remainder_order = np.argsort(-(expected - counts), kind="stable")
        for index in remainder_order[:size - int(counts.sum())]:
            counts[index] += 1
        for index in np.flatnonzero(counts == 0):
            donor = int(np.argmax(counts))
            counts[donor] -= 1
            counts[index] += 1

        digest = hashlib.sha256(f"{int(seed)}\0{sample}".encode("utf-8")).digest()
        rng = np.random.default_rng(int.from_bytes(digest[:16], "little"))
        ordered_sources = rng.permutation(unique_sources)
        group_partition = {}
        start = 0
        for partition, count in zip(PARTITIONS, counts):
            for source in ordered_sources[start:start + count]:
                group_partition[int(source)] = partition
            start += int(count)
        assignments[mask] = [group_partition[int(source)] for source in sources[mask]]
        sample_counts[str(sample)] = {
            partition: int(count) for partition, count in zip(PARTITIONS, counts)
        }
        sample_fractions[str(sample)] = {
            partition: int(count) / size for partition, count in zip(PARTITIONS, counts)
        }
    return GroupedSplit(assignments, sample_fractions, sample_counts)


def project_partition_weights(physical_weights, sample_ids, split: GroupedSplit, partition: str) -> np.ndarray:
    """Return full-length weights, zero off-partition and divided by its fraction.

    The divisor is the actual fraction of *all source events* assigned to that
    sample's partition.  It is never a fraction of selected detector outcomes
    or a ratio of signed weight sums.
    """

    weights = _weights(physical_weights)
    samples = _vector(sample_ids, "sample_ids", dtype=str)
    assignments = _vector(split.partitions, "partitions")
    if partition not in PARTITIONS:
        raise ValueError(f"unknown partition {partition!r}")
    if len(weights) != len(samples) or len(weights) != len(assignments):
        raise ValueError("weights, sample_ids and split assignments must have equal lengths")
    if np.any(~np.isin(assignments, PARTITIONS)):
        raise ValueError("split contains an unknown partition")
    result = np.zeros_like(weights)
    for sample in set(samples):
        try:
            fraction = float(split.sample_fractions[str(sample)][partition])
        except KeyError as exc:
            raise ValueError(f"missing split fraction for sample {sample!r}") from exc
        if not math.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError(f"invalid split fraction for sample {sample!r}")
        selected = (samples == sample) & (assignments == partition)
        result[selected] = weights[selected] / fraction
    return result


def class_balanced_abs_weights(labels, physical_weights) -> np.ndarray:
    """Mean-one training weights with equal class totals and zeroes preserved.

    Call on the training partition only.  Within each class, relative absolute
    physical event weights retain the expected mixture of background samples.
    """

    labels = _vector(labels, "labels")
    weights = _weights(physical_weights)
    if len(labels) != len(weights):
        raise ValueError("labels and weights must have equal lengths")
    if not len(labels):
        raise ValueError("cannot balance an empty training sample")
    if np.any(~np.isin(labels, (0, 1))) or set(labels) != {0, 1}:
        raise ValueError("training labels must contain both binary classes 0 and 1")
    absolute = np.abs(weights)
    result = np.zeros_like(absolute)
    for label in (0, 1):
        mask = labels == label
        total = float(absolute[mask].sum())
        if not math.isfinite(total) or total <= 0:
            raise ValueError(f"training class {label} has no positive absolute weight")
        result[mask] = absolute[mask] / total * (len(labels) / 2.0)
    return result


def _summary(weights, groups, selected) -> dict[str, Any]:
    selected_weights = weights[selected]
    selected_groups = groups[selected]
    group_weights = np.bincount(selected_groups, weights=selected_weights)
    total = float(math.fsum(selected_weights))
    variance = float(np.dot(group_weights, group_weights))
    return {
        "sum_weight": total,
        "variance": variance,
        "mc_error": math.sqrt(variance),
        "sum_abs_weight": float(np.abs(selected_weights).sum()),
        "positive_weight": float(selected_weights[selected_weights > 0].sum()),
        "negative_weight": float(selected_weights[selected_weights < 0].sum()),
        "effective_count": total * total / variance if variance > 0 else 0.0,
        "selected_entries": int(np.count_nonzero(selected)),
        "selected_source_events": int(len(np.unique(selected_groups))),
    }


def signed_yield_summary(weights, sample_ids, source_events, *, selected=None) -> dict[str, Any]:
    """Sum signed yields and square each selected source event's total weight.

    ``variance`` is the finite-Monte-Carlo variance, not the expected Poisson
    counting variance.  Detector hypotheses belonging to one source event are
    summed before squaring; same numeric source IDs in different samples are
    independent.  An empty selection has zero yield and zero variance.
    """

    weights = _weights(weights)
    samples, sources = _identities(sample_ids, source_events)
    if len(weights) != len(samples):
        raise ValueError("weights and event identities must have equal lengths")
    selected = np.ones(len(weights), dtype=bool) if selected is None else _vector(selected, "selected", dtype=bool)
    if len(selected) != len(weights):
        raise ValueError("selected and weights must have equal lengths")
    return _summary(weights, _group_indices(samples, sources), selected)


def scan_thresholds(
    scores,
    labels,
    weights,
    sample_ids,
    source_events,
    *,
    systematics=0.0,
    min_background_effective_events=25.0,
) -> tuple[dict[str, Any] | None, list[dict[str, Any]]]:
    """Choose a threshold using eligible validation rows and projected weights.

    Scan 501 inclusive thresholds from 0 to 1 with acceptance ``score >= t``.
    Require positive signed signal/background and sufficient background Monte
    Carlo effective events.  The highest ``S/sqrt(B + (systematics*B)**2)``
    wins; exact ties retain the lowest threshold.  An invalid point has null
    significance and an explicit reason.  No valid point returns ``None``.
    """

    scores = _vector(scores, "scores", dtype=float)
    labels = _vector(labels, "labels")
    weights = _weights(weights)
    samples, sources = _identities(sample_ids, source_events)
    if not (len(scores) == len(labels) == len(weights) == len(samples)):
        raise ValueError("scores, labels, weights and identities must have equal lengths")
    if not np.all(np.isfinite(scores)) or np.any((scores < 0) | (scores > 1)):
        raise ValueError("scores must be finite and between zero and one")
    if np.any(~np.isin(labels, (0, 1))):
        raise ValueError("labels must be binary 0 and 1")
    systematics = float(systematics)
    minimum = float(min_background_effective_events)
    if not math.isfinite(systematics) or systematics < 0:
        raise ValueError("systematics must be finite and non-negative")
    if not math.isfinite(minimum) or minimum <= 0:
        raise ValueError("min_background_effective_events must be finite and positive")
    groups = _group_indices(samples, sources)
    signal_mask, background_mask = labels == 1, labels == 0
    if len(groups):
        signal_groups = set(groups[signal_mask])
        if signal_groups.intersection(groups[background_mask]):
            raise ValueError("a source event cannot have both signal and background labels")
    total_signal = float(math.fsum(weights[signal_mask]))
    total_background = float(math.fsum(weights[background_mask]))
    scan = []
    best = None
    for threshold in np.linspace(0.0, 1.0, 501):
        selected = scores >= threshold
        signal = _summary(weights, groups, selected & signal_mask)
        background = _summary(weights, groups, selected & background_mask)
        signal_yield, background_yield = signal["sum_weight"], background["sum_weight"]
        reason = None
        if signal_yield <= 0:
            reason = "nonpositive_signal"
        elif background_yield <= 0:
            reason = "nonpositive_background"
        elif background["effective_count"] < minimum and not math.isclose(
            background["effective_count"], minimum, rel_tol=1.e-12, abs_tol=0.0,
        ):
            reason = "insufficient_background_effective_events"
        significance = (
            signal_yield / math.sqrt(background_yield + (systematics * background_yield) ** 2)
            if reason is None else None
        )
        row = {
            "threshold": float(threshold),
            "signal_events": signal_yield,
            "background_events": background_yield,
            "signal_variance": signal["variance"],
            "background_variance": background["variance"],
            "background_effective_count": background["effective_count"],
            "signal_efficiency": signal_yield / total_signal if total_signal > 0 else None,
            "background_efficiency": background_yield / total_background if total_background > 0 else None,
            "significance": significance,
            "valid": reason is None,
            "invalid_reason": reason,
        }
        scan.append(row)
        if reason is None and (best is None or significance > best["significance"]):
            best = row.copy()
    return best, scan
