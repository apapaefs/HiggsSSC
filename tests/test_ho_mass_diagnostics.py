"""Mass selection, signed MC errors and portable score exports; no ML fit needed."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from hgammagamma import plot_ho_mass_diagnostics as subject


class MassDiagnosticsTests(unittest.TestCase):
    def test_only_mass_cuts_are_removed_and_score_selection_matches_training_scope(self):
        data = {"m_gg": np.array([119., 120., 125., 130., 131., 125.]),
                "pt_gamma1": np.array([60., 60., 30., 60., 60., 60.]),
                "n_selected_photons": np.array([2., 2., 2., 2., 2., 1.]),
                "score": np.array([.6, .5, .8, .4, .7, .9])}
        cuts = [{"variable": "m_gg", "min": 120}, {"variable": "m_gg", "max": 130},
                {"variable": "pt_gamma1", "min": 50}]
        masks = subject.selection_masks(data, cuts, .5)
        np.testing.assert_array_equal(masks["other_cuts"], [1, 1, 0, 1, 1, 0])
        np.testing.assert_array_equal(masks["all_cuts"], [0, 1, 0, 1, 0, 0])
        np.testing.assert_array_equal(masks["xgboost"], [1, 1, 1, 0, 1, 0])
        self.assertNotIn("xgboost", subject.selection_masks(data, cuts, None))

    def test_histogram_groups_response_hypotheses_per_source_and_sample(self):
        # Two rows of the same source cancel in bin 0. The same integer ID in
        # another sample is independent; the upper endpoint belongs to bin 1.
        values, errors = subject.histogram(
            np.array([.2, .3, .4, 1.2, 2., -1., 3.]),
            np.array([3., -1., 4., -2., 1., 100., 100.]),
            np.array(["a", "a", "b", "a", "a", "a", "a"]),
            np.array([1, 1, 1, 2, 3, 4, 5]), np.array([0., 1., 2.]))
        np.testing.assert_allclose(values, [6., -1.])
        np.testing.assert_allclose(errors, [np.sqrt(20), np.sqrt(5)])
        values, errors = subject.histogram(np.array([]), np.array([]), np.array([], dtype=str),
                                           np.array([], dtype=int), np.array([0., 1., 2.]))
        np.testing.assert_array_equal(values, [0, 0])
        np.testing.assert_array_equal(errors, [0, 0])

    def test_unbinned_width_distinguishes_signed_moments_and_absolute_quantiles(self):
        result = subject.mass_summary(np.array([124., 125., 126.]), np.array([2., -1., 2.]))
        self.assertAlmostEqual(result["mean_gev"], 125.)
        self.assertAlmostEqual(result["rms_gev"], np.sqrt(4/3))
        self.assertAlmostEqual(result["negative_absolute_weight_fraction"], .2)
        self.assertEqual(result["quantile_weights"], "absolute physical weights")
        self.assertEqual(result["signed_events"], 3.)
        self.assertIsNone(subject.mass_summary(np.array([]), np.array([]))["sigma68_gev"])

    def fixture(self, directory):
        records = []
        for label, sample in enumerate(("background", "signal")):
            for entry, partition in enumerate(("train", "validation", "test")):
                record = dict(sample=sample, partition=partition, sourceevent=entry, tree_entry=entry,
                              label=label, score=.7, raw_signed_weight=2., physical_weight=2.,
                              partition_fraction=1/3, projected_physical_weight=6.)
                record.update({name: 1. for name in subject.FEATURE_NAMES})
                record.update(m_gg=125., n_selected_photons=2.)
                records.append(record)
        metadata = dict(luminosity_fb=1., eligible_rows_per_partition={p: 2 for p in ("train", "validation", "test")},
                        sample_split_fractions={s: {p: 1/3 for p in ("train", "validation", "test")}
                                                for s in ("background", "signal")},
                        normalization_provenance=[dict(sample=s, cross_section_pb=.001, weight_scale=1.,
                                                       normalization_sum_weight=1.) for s in ("background", "signal")])
        (directory / "model_metadata.json").write_text(json.dumps(metadata))
        self.save_scores(directory, records)
        return records

    def save_scores(self, directory, records):
        with (directory / "scores.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    def test_embedded_variables_work_without_original_root_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.fixture(directory)
            _, data, provenance = subject.load_events(directory)
            self.assertEqual(provenance["variables_source"], "scores.csv")
            self.assertEqual(provenance["root_files"], {})
            np.testing.assert_array_equal(data["m_gg"], [125.] * 6)
            np.testing.assert_array_equal(data["projected_physical_weight"], [6.] * 6)

    def test_incomplete_duplicate_and_misnormalized_exports_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            records = self.fixture(directory)
            for corrupted in (records[:-1], records + [records[0]],
                              [{**r, "projected_physical_weight": 99.} for r in records]):
                self.save_scores(directory, corrupted)
                with self.assertRaises(ValueError):
                    subject.load_events(directory)


if __name__ == "__main__":
    unittest.main()
