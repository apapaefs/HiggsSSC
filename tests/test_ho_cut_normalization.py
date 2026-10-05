"""HO cut yields with mocked tree reads; no ROOT or event generation needed.

The reader is replaced only while importing an isolated copy of the subject.
These tests never open the synthetic ``_var.root`` placeholders.
"""

from contextlib import ExitStack
import importlib.util
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import ModuleType
import unittest
from unittest.mock import Mock, patch

from hgammagamma import make_gammagamma_report as report


FEATURES = ["m_gg", "pt_gamma1", "eta_gamma1", "pt_gamma2", "eta_gamma2",
            "deltaR_gg", "deltaPhi_gg", "pt_gg", "y_gg", "n_selected_photons"]
BR = 0.00227
RUN_TAG = "test"


def load_cut_subject():
    reader = ModuleType("read_root_varfiles")
    reader.FEATURE_NAMES = FEATURES
    reader.read_named_ROOT_varfile = Mock()
    name = "_ho_cut_normalization_subject"
    path = Path(__file__).resolve().parents[1] / "analyze_lo_varfiles.py"
    spec = importlib.util.spec_from_file_location(name, path)
    subject = importlib.util.module_from_spec(spec)
    with patch.dict(sys.modules, {"read_root_varfiles": reader, name: subject}):
        spec.loader.exec_module(subject)
    return subject


class HOCutSummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subject = load_cut_subject()

    def sample(self, *, signal=True, **overrides):
        values = dict(
            name="signal_gg_h_aa" if signal else "bkg_prompt_aa",
            category="Signal" if signal else "Backgrounds",
            sample_dir=Path("sample"), var_file=Path("sample_var.root"),
            dat_file=Path("sample.dat"), cross_section_pb=200. if signal else 50.,
            cross_section_error_pb=.01, weight_scale=BR if signal else 1.,
            events_read=3., sum_weight=16. * BR if signal else 16.,
            analysis_name="SSC_GEM_weighted_response", detector_response="ssc",
            response_mode="genuine", normalization_sum_weight=20. * BR if signal else 20.,
            normalization_kind="ihixs_n3lo" if signal else "generator",
            ihixs_record_sha256="a" * 64 if signal else "",
            shower_quality={"termination": "source_exhausted", "attempted_events": 4,
                            "saved_events": 3, "discarded_events": 1},
            requires_full_sample=True,
            tree_entries=3, sum_tree_weight=16. * BR if signal else 16.,
            sum_abs_weight=20. * BR if signal else 20.,
        )
        values.update(overrides)
        return self.subject.SampleInfo(**values)

    def rows(self):
        return [{"m_gg": 125.}, {"m_gg": 124.}, {"m_gg": 145.}]

    def summarize(self, sample, weights, rows=None, **kwargs):
        with patch.object(self.subject, "read_named_ROOT_varfile",
                          return_value=(self.rows() if rows is None else rows, weights)):
            return self.subject.summarize_sample(
                sample, [self.subject.Cut("m_gg", 123., 127.)], 10., **kwargs)

    def test_ihixs_rate_uses_signed_source_denominator_and_br_once(self):
        result = self.summarize(self.sample(), [10. * BR, -2. * BR, 8. * BR])
        self.assertEqual(result["selected_entries"], 2)
        self.assertAlmostEqual(result["sum_weight"], 16. * BR)
        self.assertAlmostEqual(result["sum_selected_weight"], 8. * BR)
        self.assertAlmostEqual(result["analysis_efficiency"], .4)
        self.assertAlmostEqual(result["selected_cross_section_pb"], 200. * BR * .4)
        self.assertAlmostEqual(result["expected_events"], 200. * BR * .4 * 10. * 1000.)
        self.assertEqual(result["normalization_kind"], "ihixs_n3lo")
        self.assertEqual(result["ihixs_record_sha256"], "a" * 64)
        self.assertAlmostEqual(result["normalization_sum_weight"], 20. * BR)

    def test_negative_selected_weight_yield_remains_signed(self):
        rows = [{"m_gg": 145.}, {"m_gg": 124.}, {"m_gg": 145.}]
        result = self.summarize(self.sample(), [10. * BR, -2. * BR, 8. * BR], rows)
        self.assertAlmostEqual(result["analysis_efficiency"], -.1)
        self.assertAlmostEqual(result["selected_cross_section_pb"], -20. * BR)
        self.assertAlmostEqual(result["expected_events"], -20. * BR * 10. * 1000.)

    def test_nlo_background_keeps_its_signed_generator_rate(self):
        result = self.summarize(self.sample(signal=False), [10., -2., 8.])
        self.assertAlmostEqual(result["selected_cross_section_pb"], 50. * .4)
        self.assertAlmostEqual(result["expected_events"], 50. * .4 * 10. * 1000.)
        self.assertEqual(result["normalization_kind"], "generator")

    def test_ho_partial_tree_limit_is_rejected_before_reading(self):
        for signal in (False, True):
            with self.subTest(signal=signal):
                with patch.object(self.subject, "read_named_ROOT_varfile") as reader:
                    with self.assertRaises(ValueError):
                        self.subject.summarize_sample(self.sample(signal=signal), [], 10., max_events=1)
                    reader.assert_not_called()

    def test_tree_weight_sum_must_match_saved_ho_analysis(self):
        for signal in (False, True):
            scale = BR if signal else 1.
            with self.subTest(signal=signal):
                with self.assertRaises(ValueError):
                    self.summarize(self.sample(signal=signal), [10. * scale, -2. * scale, 7. * scale])

    def test_response_hypotheses_need_not_have_one_row_per_saved_event(self):
        rows = self.rows() + [{"m_gg": -999.}]
        result = self.summarize(self.sample(tree_entries=4), [9. * BR, -2. * BR, 8. * BR, BR], rows)
        self.assertEqual(result["entries_read"], 4)
        self.assertAlmostEqual(result["sum_weight"], 16. * BR)
        self.assertAlmostEqual(result["selected_cross_section_pb"], 200. * BR * 7. / 20.)

    def test_ho_reader_uses_strict_full_tree_validation(self):
        with patch.object(self.subject, "read_named_ROOT_varfile",
                          return_value=(self.rows(), [10. * BR, -2. * BR, 8. * BR])) as reader:
            sample = self.sample()
            self.subject.summarize_sample(sample, [], 10.)
            reader.assert_called_once_with(sample.var_file, max_events=None, strict=True)

    def test_missing_tree_rows_cannot_be_normalized_as_a_complete_population(self):
        with self.assertRaises(ValueError):
            self.summarize(self.sample(tree_entries=4), [10. * BR, -2. * BR, 8. * BR])

    def test_ho_source_denominator_must_be_positive_and_finite(self):
        for value in (0., -1., float("nan"), float("inf")):
            with self.subTest(denominator=value), self.assertRaises(ValueError):
                self.summarize(self.sample(normalization_sum_weight=value), [10. * BR, -2. * BR, 8. * BR])

    def test_legacy_lo_partial_reads_keep_their_existing_normalization(self):
        sample = self.sample(requires_full_sample=False, normalization_kind="generator",
                             normalization_sum_weight=None, ihixs_record_sha256="")
        rows = self.rows()[:2]
        result = self.summarize(sample, [10. * BR, -2. * BR], rows, max_events=2)
        self.assertAlmostEqual(result["analysis_efficiency"], 1.)
        self.assertAlmostEqual(result["selected_cross_section_pb"], 200. * BR)

    def test_ho_xgboost_is_rejected_before_loading_the_lo_classifier(self):
        inputs = ({}, Path("campaign"), "classifier", RUN_TAG, 10., [self.sample()], Path("out"))
        with patch.object(self.subject, "load_analysis_inputs", return_value=inputs):
            with self.assertRaisesRegex(ValueError, "HO campaigns currently support the cuts"):
                self.subject.run_xgboost(Path("card.yaml"))


class HOCutDiscoveryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.subject = load_cut_subject()

    def fixture(self, root, *, signal=True):
        name = "signal_gg_h_aa" if signal else "bkg_prompt_aa"
        category = "Signal" if signal else "Backgrounds"
        directory = Path(root) / category / "events" / name
        directory.mkdir(parents=True)
        scale = BR if signal else 1.
        dat = {
            "analysis": "SSC_GEM_weighted_response", "detector_response": "ssc",
            "response_mode": "genuine", "events_read": "3", "weight_scale": str(scale),
            "sum_weight": str(16. * scale), "sum_diphoton_weight": str(8. * scale),
            "tree_entries": "3", "sum_tree_weight": str(16. * scale),
            "sum_abs_weight": str(20. * scale),
        }
        dat_file = directory / f"sample-{RUN_TAG}.dat"
        dat_file.write_text("".join(f"{key} {value}\n" for key, value in dat.items()))
        (directory / f"sample-{RUN_TAG}_var.root").write_bytes(b"mocked tree read")
        completion = {
            "fingerprint": "b" * 64, "termination": "source_exhausted",
            "requested_events": 4, "lhe_events": 4, "attempted_events": 4,
            "generated_events": 3, "saved_events": 3, "discarded_events": 1,
            "exception_counts": {}, "max_momentum_violation_mev": 1.,
        }
        manifest = {
            "sample": name, "run_tag": RUN_TAG, "nevents_requested": 4,
            "completed": ["generate", "shower", "analyze"],
            "hard_accuracy": "NNLO+PS (HJMiNNLO)" if signal else "NLO QCD + PS",
            "matching": "POWHEG" if signal else "MC@NLO",
            "process": "HJMiNNLO; h -> gamma gamma" if signal else "p p > a a [QCD]",
            "weight_scale": scale, "analysis": dat, "shower_completion": completion,
            "lhe": {"cross_section_pb": 190. if signal else 50.},
        }
        (directory / "campaign.json").write_text(json.dumps(manifest))
        sidecar = {
            "run_tag": RUN_TAG, "sample": name,
            "cross_section_pb": 200. if signal else 50., "cross_section_error_pb": .01,
            "weight_scale": scale, "shower_completion_sha256": completion["fingerprint"],
        }
        (directory / f"normalization-{RUN_TAG}.json").write_text(json.dumps(sidecar))
        return directory, manifest, sidecar

    def mocked_proofs(self, *, signal=True):
        """Keep discovery fixtures small; proof validation has its own test suites."""
        stack = ExitStack()
        sidecar = {"normalization_kind": "ihixs_n3lo", "cross_section_pb": 200.,
                   "cross_section_error_pb": .01, "native_cross_section_pb": 190.,
                   "weight_scale": BR, "ihixs_record_sha256": "a" * 64,
                   "shower_completion_sha256": "b" * 64,
                   "ihixs": {"uncertainties": {}}} if signal else None
        stack.enter_context(patch.object(report, "load_ho_signal_sidecar", return_value=sidecar))
        stack.enter_context(patch.object(report, "parse_cross_section",
                                        return_value=(200. if signal else 50., .01)))
        stack.enter_context(patch.object(report, "validate_ho_analysis_summary"))
        stack.enter_context(patch.object(report.ho_shower, "validate_completion"))
        stack.enter_context(patch.object(report.ho_shower, "expected_analysis_events", return_value=3))
        stack.enter_context(patch.object(report.ho_shower, "normalization_denominator",
                                        return_value=20. * (BR if signal else 1.)))
        return stack

    def test_signal_discovers_ihixs_rate_without_the_lo_k_factor(self):
        with TemporaryDirectory() as tmp:
            self.fixture(tmp)
            with self.mocked_proofs():
                sample = self.subject.discover_samples(Path(tmp), RUN_TAG)[0]
            self.assertEqual(sample.cross_section_pb, 200.)
            self.assertEqual(sample.weight_scale, BR)
            self.assertEqual(sample.normalization_kind, "ihixs_n3lo")
            self.assertEqual(sample.ihixs_record_sha256, "a" * 64)
            self.assertAlmostEqual(sample.normalization_sum_weight, 20. * BR)
            self.assertTrue(sample.requires_full_sample)
            self.assertEqual(sample.shower_quality["discarded_events"], 1)

    def test_nlo_background_discovery_keeps_generator_rate_and_completion(self):
        with TemporaryDirectory() as tmp:
            self.fixture(tmp, signal=False)
            with self.mocked_proofs(signal=False):
                sample = self.subject.discover_samples(Path(tmp), RUN_TAG)[0]
            self.assertEqual(sample.cross_section_pb, 50.)
            self.assertEqual(sample.weight_scale, 1.)
            self.assertEqual(sample.normalization_kind, "generator")
            self.assertEqual(sample.ihixs_record_sha256, "")
            self.assertEqual(sample.normalization_sum_weight, 20.)
            self.assertTrue(sample.requires_full_sample)

    def test_yaml_sample_and_category_factors_cannot_change_ho_rates(self):
        for signal in (False, True):
            for category_key in (False, True):
                with self.subTest(signal=signal, category_key=category_key), TemporaryDirectory() as tmp:
                    directory, _, _ = self.fixture(tmp, signal=signal)
                    key = ("Signal" if signal else "Backgrounds") if category_key else directory.name
                    factor = 2. * (BR if signal else 1.)
                    with self.mocked_proofs(signal=signal):
                        with self.assertRaises(ValueError):
                            self.subject.discover_samples(Path(tmp), RUN_TAG, rate_factors={key: factor})

    def test_matching_sample_factor_does_not_hide_a_conflicting_category_factor(self):
        with TemporaryDirectory() as tmp:
            self.fixture(tmp)
            with self.mocked_proofs():
                with self.assertRaises(ValueError):
                    self.subject.discover_samples(
                        Path(tmp), RUN_TAG, rate_factors={"signal_gg_h_aa": BR, "Signal": 2. * BR})

    def test_saved_ho_scale_cannot_be_changed_in_dat_metadata(self):
        with TemporaryDirectory() as tmp:
            directory, _, _ = self.fixture(tmp)
            dat_file = directory / f"sample-{RUN_TAG}.dat"
            dat_file.write_text(dat_file.read_text().replace(f"weight_scale {BR}\n", "weight_scale 1\n"))
            with self.mocked_proofs():
                with self.assertRaisesRegex(ValueError, "weight scale"):
                    self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_missing_ihixs_signal_sidecar_cannot_fall_back_to_a_generator_rate(self):
        with TemporaryDirectory() as tmp:
            directory, _, _ = self.fixture(tmp)
            (directory / f"normalization-{RUN_TAG}.json").unlink()
            with self.assertRaises(ValueError):
                self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_altered_analysis_summary_is_rejected(self):
        for signal in (False, True):
            with self.subTest(signal=signal), TemporaryDirectory() as tmp:
                self.fixture(tmp, signal=signal)
                with self.mocked_proofs(signal=signal):
                    with patch.object(report, "validate_ho_analysis_summary", side_effect=ValueError("changed summary")):
                        with self.assertRaisesRegex(ValueError, "changed summary"):
                            self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_completion_failure_is_not_silently_ignored(self):
        with TemporaryDirectory() as tmp:
            self.fixture(tmp, signal=False)
            with self.mocked_proofs(signal=False):
                with patch.object(report.ho_shower, "validate_completion", side_effect=ValueError("changed shower")):
                    with self.assertRaisesRegex(ValueError, "changed shower"):
                        self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_rate_sidecar_must_match_its_saved_shower_fingerprint(self):
        with TemporaryDirectory() as tmp:
            directory, _, sidecar = self.fixture(tmp, signal=False)
            sidecar["shower_completion_sha256"] = "c" * 64
            (directory / f"normalization-{RUN_TAG}.json").write_text(json.dumps(sidecar))
            with self.mocked_proofs(signal=False):
                with self.assertRaisesRegex(ValueError, "saved shower population"):
                    self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_missing_explicitly_requested_sample_is_rejected(self):
        with TemporaryDirectory() as tmp:
            self.fixture(tmp)
            with self.mocked_proofs():
                with self.assertRaisesRegex(FileNotFoundError, "bkg_gamma_j"):
                    self.subject.discover_samples(
                        Path(tmp), RUN_TAG, requested={"signal_gg_h_aa", "bkg_gamma_j"})

    def test_ho_cut_sample_requires_a_completed_analysis(self):
        with TemporaryDirectory() as tmp:
            directory, manifest, _ = self.fixture(tmp)
            manifest["completed"] = ["generate", "shower"]
            (directory / "campaign.json").write_text(json.dumps(manifest))
            with self.mocked_proofs():
                with self.assertRaisesRegex(ValueError, "completed analysis"):
                    self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_nonfinite_or_inconsistent_weight_summary_is_rejected(self):
        for field, value in (("sum_abs_weight", "inf"), ("sum_abs_weight", "nan"),
                             ("sum_abs_weight", "-1"), ("sum_abs_weight", "0"),
                             ("sum_tree_weight", "inf"), ("sum_tree_weight", "nan")):
            with self.subTest(field=field, value=value), TemporaryDirectory() as tmp:
                directory, _, _ = self.fixture(tmp, signal=False)
                path = directory / f"sample-{RUN_TAG}.dat"
                lines = path.read_text().splitlines()
                path.write_text("\n".join(f"{field} {value}" if line.startswith(field + " ") else line
                                           for line in lines) + "\n")
                with self.mocked_proofs(signal=False):
                    with self.assertRaisesRegex(ValueError, "signed-weight summary"):
                        self.subject.discover_samples(Path(tmp), RUN_TAG)

    def test_same_signal_name_without_ho_manifest_keeps_legacy_lo_default(self):
        with TemporaryDirectory() as tmp:
            directory, _, _ = self.fixture(tmp)
            (directory / "campaign.json").unlink()
            (directory / f"normalization-{RUN_TAG}.json").unlink()
            dat_file = directory / f"sample-{RUN_TAG}.dat"
            dat_file.write_text(dat_file.read_text().replace(f"weight_scale {BR}\n", "weight_scale 1\n"))
            banner_dir = directory / "mg5_process/Events" / RUN_TAG
            banner_dir.mkdir(parents=True)
            (banner_dir / f"{RUN_TAG}_tag_1_banner.txt").write_text("Integrated weight (pb) : 50.0\n")
            sample = self.subject.discover_samples(Path(tmp), RUN_TAG)[0]
            self.assertEqual(sample.cross_section_pb, 50.)
            self.assertEqual(sample.weight_scale, 2. * BR)
            self.assertEqual(sample.normalization_kind, "generator")
            self.assertIsNone(sample.normalization_sum_weight)
            self.assertFalse(sample.requires_full_sample)


if __name__ == "__main__":
    unittest.main()
