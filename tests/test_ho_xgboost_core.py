"""Source-event split and signed-weight checks without ROOT or ML packages."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ho_xgboost_core import (
    PARTITIONS,
    class_balanced_abs_weights,
    grouped_sample_split,
    project_partition_weights,
    scan_thresholds,
    signed_yield_summary,
)


class GroupedSplitTests(unittest.TestCase):
    def test_deterministic_disjoint_groups_and_actual_fractions(self):
        samples = np.array(["signal"] * 33 + ["background"] * 21)
        sources = np.concatenate((np.repeat(np.arange(11), 3), np.repeat(np.arange(7), 3)))
        split = grouped_sample_split(samples, sources)
        again = grouped_sample_split(samples, sources)
        np.testing.assert_array_equal(split.partitions, again.partitions)
        self.assertEqual(split.sample_group_counts["signal"], {"train": 7, "validation": 2, "test": 2})
        for sample in ("signal", "background"):
            total = len(set(sources[samples == sample]))
            for source in set(sources[samples == sample]):
                self.assertEqual(len(set(split.partitions[(samples == sample) & (sources == source)])), 1)
            for partition in PARTITIONS:
                mask = (samples == sample) & (split.partitions == partition)
                self.assertEqual(split.sample_fractions[sample][partition], len(set(sources[mask])) / total)
        np.testing.assert_array_equal(
            grouped_sample_split(samples[::-1], sources[::-1]).partitions[::-1], split.partitions,
        )
        signal_only = grouped_sample_split(samples[:33], sources[:33])
        np.testing.assert_array_equal(signal_only.partitions, split.partitions[:33])

    def test_minimum_three_groups_and_invalid_fractions(self):
        split = grouped_sample_split(["a"] * 3, [0, 1, 2])
        self.assertEqual(split.sample_group_counts["a"], dict.fromkeys(PARTITIONS, 1))
        with self.assertRaisesRegex(ValueError, "at least three"):
            grouped_sample_split(["a", "a", "a"], [0, 0, 1])
        for fractions in ((.6, .2, .1), (.8, .2, 0), (.5, .5)):
            with self.assertRaises(ValueError):
                grouped_sample_split(["a"] * 3, [0, 1, 2], fractions=fractions)

    def test_projection_uses_full_group_fraction_and_closes_constant_yield(self):
        samples = ["signal"] * 11 + ["background"] * 7
        sources = list(range(11)) + list(range(7))
        weights = np.array([2.] * 11 + [-3.] * 7)
        split = grouped_sample_split(samples, sources)
        for partition in PARTITIONS:
            projected = project_partition_weights(weights, samples, split, partition)
            self.assertAlmostEqual(projected[:11].sum(), 22.)
            self.assertAlmostEqual(projected[11:].sum(), -21.)
            self.assertTrue(np.all(projected[split.partitions != partition] == 0))
            # Eligibility removes rows only after the full-tree split/projection.
            eligible = np.arange(len(samples)) % 2 == 0
            retained = eligible & (split.partitions == partition)
            for index in np.flatnonzero(retained):
                self.assertEqual(projected[index], weights[index] / split.sample_fractions[samples[index]][partition])

    def test_bad_source_identity_rejected(self):
        for sources in ([0., 1., 2.], [0, -1, 2], [True, False, True]):
            with self.assertRaisesRegex(ValueError, "source_events"):
                grouped_sample_split(["a"] * 3, sources)


class BalancedWeightTests(unittest.TestCase):
    def test_equal_class_totals_zero_preserved_and_negative_abs_used(self):
        labels = [1, 1, 0, 0, 0]
        balanced = class_balanced_abs_weights(labels, [2., -2., 0., 6., -2.])
        self.assertAlmostEqual(balanced[:2].sum(), 2.5)
        self.assertAlmostEqual(balanced[2:].sum(), 2.5)
        self.assertEqual(balanced[2], 0.)
        self.assertEqual(balanced[3] / balanced[4], 3.)
        self.assertAlmostEqual(balanced.mean(), 1.)

    def test_relative_physical_background_mixture_is_preserved(self):
        labels = np.array([1, 1, 0, 0, 0, 0])
        base = class_balanced_abs_weights(labels, [1, 1, 10, 10, 1, 1])
        changed = class_balanced_abs_weights(labels, [1, 1, 20, 20, 1, 1])
        self.assertEqual(base[2:4].sum() / base[4:].sum(), 10.)
        self.assertEqual(changed[2:4].sum() / changed[4:].sum(), 20.)
        self.assertAlmostEqual(changed[labels == 0].sum(), changed[labels == 1].sum())

    def test_zero_class_nonfinite_and_missing_class_rejected(self):
        for labels, weights in (([0, 1], [1, 0]), ([0, 1], [1, np.nan]), ([1, 1], [1, 1])):
            with self.assertRaises(ValueError):
                class_balanced_abs_weights(labels, weights)


class SignedYieldTests(unittest.TestCase):
    def test_correlated_hypotheses_sum_before_squaring(self):
        summary = signed_yield_summary([2., 3., -1., 4.], ["a", "a", "a", "b"], [0, 0, 1, 0])
        self.assertEqual(summary["sum_weight"], 8.)
        self.assertEqual(summary["variance"], 25. + 1. + 16.)
        self.assertEqual(summary["positive_weight"], 9.)
        self.assertEqual(summary["negative_weight"], -1.)
        self.assertEqual(summary["sum_abs_weight"], 10.)
        self.assertEqual(summary["selected_source_events"], 3)
        self.assertEqual(summary["selected_entries"], 4)

    def test_negative_and_empty_selected_yields_are_retained(self):
        summary = signed_yield_summary([2., -3.], ["a", "a"], [0, 1])
        self.assertEqual(summary["sum_weight"], -1.)
        self.assertEqual(summary["variance"], 13.)
        empty = signed_yield_summary([2., -3.], ["a", "a"], [0, 1], selected=[False, False])
        self.assertEqual(empty["sum_weight"], 0.)
        self.assertEqual(empty["variance"], 0.)
        self.assertEqual(empty["effective_count"], 0.)


class ThresholdTests(unittest.TestCase):
    def test_best_threshold_uses_signed_yield_and_systematics(self):
        # At score > .2 the twenty low-score background events disappear.
        scores = [.9] * 5 + [.8] * 30 + [.2] * 20
        labels = [1] * 5 + [0] * 50
        weights = [2.] * 5 + [1.] * 50
        samples = ["signal"] * 5 + ["background"] * 50
        sources = list(range(5)) + list(range(50))
        best, scan = scan_thresholds(scores, labels, weights, samples, sources, systematics=.1)
        self.assertEqual(len(scan), 501)
        self.assertIsNotNone(best)
        self.assertAlmostEqual(best["threshold"], .202)
        self.assertEqual(best["signal_events"], 10.)
        self.assertEqual(best["background_events"], 30.)
        self.assertAlmostEqual(best["significance"], 10. / np.sqrt(30. + 9.))
        self.assertEqual(best["background_effective_count"], 30.)
        self.assertFalse(scan[-1]["valid"])

    def test_no_valid_threshold_for_low_statistics(self):
        best, scan = scan_thresholds([.9] + [.8] * 24, [1] + [0] * 24,
                                     [1.] * 25, ["s"] + ["b"] * 24, [0] + list(range(24)))
        self.assertIsNone(best)
        self.assertEqual(scan[0]["invalid_reason"], "insufficient_background_effective_events")
        self.assertTrue(all(row["significance"] is None for row in scan))

    def test_exact_effective_count_boundary_survives_roundoff(self):
        best, _ = scan_thresholds([.9] + [.8] * 25, [1] + [0] * 25,
                                  [1.] + [.1] * 25, ["s"] + ["b"] * 25,
                                  [0] + list(range(25)))
        self.assertIsNotNone(best)
        self.assertAlmostEqual(best["background_effective_count"], 25.)

    def test_valid_negative_weights_use_signed_yields_and_squared_variance(self):
        best, _ = scan_thresholds([.9] + [.8] * 125, [1] + [0] * 125,
                                  [1.] + [.1] * 100 + [-.1] * 25,
                                  ["s"] + ["b"] * 125, [0] + list(range(125)))
        self.assertIsNotNone(best)
        self.assertEqual(best["threshold"], 0.)
        self.assertAlmostEqual(best["background_events"], 7.5)
        self.assertAlmostEqual(best["background_variance"], 1.25)
        self.assertAlmostEqual(best["background_effective_count"], 45.)

    def test_grouped_variance_prevents_false_statistics(self):
        # Fifty detector rows represent only ten independent source events.
        best, scan = scan_thresholds([.9] + [.8] * 50, [1] + [0] * 50,
                                     [1.] * 51, ["s"] + ["b"] * 50,
                                     [0] + list(np.repeat(np.arange(10), 5)))
        self.assertIsNone(best)
        self.assertEqual(scan[0]["background_effective_count"], 10.)

    def test_negative_signal_and_background_fail_without_clipping(self):
        for weights, reason in (([-2., 1.], "nonpositive_signal"), ([2., -1.], "nonpositive_background")):
            best, scan = scan_thresholds([.9, .8], [1, 0], weights, ["s", "b"], [0, 0])
            self.assertIsNone(best)
            self.assertEqual(scan[0]["invalid_reason"], reason)
            self.assertEqual(scan[0]["signal_events"], weights[0])
            self.assertEqual(scan[0]["background_events"], weights[1])

    def test_empty_input_returns_complete_invalid_scan(self):
        best, scan = scan_thresholds([], [], [], [], [])
        self.assertIsNone(best)
        self.assertEqual(len(scan), 501)

    def test_mixed_labels_for_source_and_nonfinite_scores_rejected(self):
        with self.assertRaisesRegex(ValueError, "both signal and background"):
            scan_thresholds([.9, .8], [1, 0], [1., 1.], ["a", "a"], [0, 0])
        with self.assertRaisesRegex(ValueError, "scores"):
            scan_thresholds([np.nan], [1], [1.], ["a"], [0])


if __name__ == "__main__":
    unittest.main()
