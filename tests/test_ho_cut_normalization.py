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

    def test_ho_xgboost_dispatches_to_the_separate_ho_workflow(self):
        analysis = {"_resolved_detector_response": "ssc"}
        samples = [self.sample(), self.sample(signal=False)]
        inputs = (analysis, Path("campaign"), "classifier", RUN_TAG, 10., samples, Path("out"))
        result = object()
        with patch.object(self.subject, "load_analysis_inputs", return_value=inputs) as loader:
            with patch.object(self.subject, "run_ho_xgboost", return_value=result) as ho_runner:
                with patch.dict(sys.modules, {"xgboost_root_varfiles_module": None}):
                    actual = self.subject.run_xgboost(
                        Path("card.yaml"), run_tag_override="override", progress_enabled=False)
        self.assertIs(actual, result)
        loader.assert_called_once_with(Path("card.yaml"), "override")
        ho_runner.assert_called_once_with(
            analysis, "classifier", RUN_TAG, 10., samples, Path("out"), progress_enabled=False)


class HOXGBoostCLITests(unittest.TestCase):
    """Exercise HO orchestration with real population checks and a mocked trainer."""

    @classmethod
    def setUpClass(cls):
        cls.subject = load_cut_subject()

    sample = HOCutSummaryTests.sample

    def samples(self):
        return [self.sample(tree_entries=4, var_file=Path("signal_var.root")),
                self.sample(signal=False, tree_entries=4, var_file=Path("background_var.root"))]

    def population(self, sample):
        selected = dict(zip(FEATURES, [125., 60., -.5, 50., .5, 2.1, 2.9, 25., -.4, 2.]))
        unselected = dict(zip(FEATURES, [-999.] * 9 + [0.]))
        rows = [unselected, selected.copy(), selected.copy(), selected.copy()]
        scale = sample.weight_scale
        return rows, [3. * scale, 7. * scale, -2. * scale, 8. * scale], [0, 0, 1, 2], [0, 1, 2, 3]

    def modules(self, samples):
        populations = {sample.var_file: self.population(sample) for sample in samples}
        reader = ModuleType("read_root_varfiles")
        reader.read_ho_ROOT_varfile = Mock(side_effect=lambda path: populations[path])
        trainer = ModuleType("ho_xgboost_analysis")
        trainer.validate_config = Mock()
        trainer.run_ho_signal_background_analysis = Mock(return_value={
            "metadata": {"name": "classifier"}, "summary_rows": [], "assets": [],
        })
        return reader, trainer, populations

    def run_workflow(self, analysis=None, samples=None):
        return self.subject.run_ho_xgboost(
            {"_resolved_detector_response": "ssc"} if analysis is None else analysis,
            "classifier", RUN_TAG, 10., self.samples() if samples is None else samples,
            Path("out"), progress_enabled=False)

    def test_partial_population_caps_are_rejected_before_importing_ml(self):
        for analysis in ({"max_events": 1}, {"xgboost": {"max_events": 1}}):
            with self.subTest(analysis=analysis):
                with patch.object(self.subject, "load_config") as baseline_loader:
                    with patch.dict(sys.modules, {"ho_xgboost_analysis": None,
                                                  "read_root_varfiles": None}):
                        with self.assertRaisesRegex(ValueError, "requires the full sample"):
                            self.run_workflow(analysis=analysis)
                    baseline_loader.assert_not_called()

    def test_mixed_ho_and_lo_samples_are_rejected_before_importing_ml(self):
        samples = self.samples()
        samples[1].requires_full_sample = False
        with patch.object(self.subject, "load_config") as baseline_loader:
            with patch.dict(sys.modules, {"ho_xgboost_analysis": None,
                                          "read_root_varfiles": None}):
                with self.assertRaisesRegex(ValueError, "cannot mix HO and legacy LO"):
                    self.run_workflow(samples=samples)
            baseline_loader.assert_not_called()

    def test_missing_signal_or_background_is_rejected_before_importing_ml(self):
        for sample in self.samples():
            with self.subTest(category=sample.category):
                with patch.object(self.subject, "load_config") as baseline_loader:
                    with patch.dict(sys.modules, {"ho_xgboost_analysis": None,
                                                  "read_root_varfiles": None}):
                        with self.assertRaisesRegex(ValueError, "requires signal and background"):
                            self.run_workflow(samples=[sample])
                    baseline_loader.assert_not_called()

    def test_all_full_populations_are_validated_before_training(self):
        samples = self.samples()
        reader, trainer, populations = self.modules(samples)
        original_validator = self.subject.validate_ho_sample_population
        validated = []

        def validate(sample, rows, weights):
            original_validator(sample, rows, weights)
            validated.append(sample.name)

        def train(loaded, **kwargs):
            self.assertEqual(validated, [sample.name for sample in samples])
            self.assertEqual(loaded, [(sample, *populations[sample.var_file]) for sample in samples])
            self.assertTrue(all(item[1][0]["n_selected_photons"] == 0. for item in loaded))
            return {"metadata": kwargs["metadata"], "summary_rows": [], "assets": []}

        trainer.run_ho_signal_background_analysis.side_effect = train
        baseline = {"analysis": {"cuts": [{"variable": "m_gg", "min": 123., "max": 127.}]}}
        with patch.dict(sys.modules, {"read_root_varfiles": reader, "ho_xgboost_analysis": trainer}):
            with patch.object(self.subject, "load_config", return_value=baseline):
                with patch.object(self.subject, "validate_ho_sample_population", side_effect=validate):
                    result = self.run_workflow(samples=samples)
        self.assertEqual(reader.read_ho_ROOT_varfile.call_count, len(samples))
        self.assertEqual(result.index_html, Path("out/index.html"))
        trainer.run_ho_signal_background_analysis.assert_called_once()

    def test_invalid_later_population_stops_before_training(self):
        samples = self.samples()
        for failure in ("count", "signed_sum", "source_count", "absolute_sum", "denominator"):
            with self.subTest(failure=failure):
                reader, trainer, populations = self.modules(samples)
                rows, weights, sources, entries = populations[samples[1].var_file]
                if failure == "count":
                    rows.pop()
                elif failure == "signed_sum":
                    weights[-1] -= 1.
                elif failure == "source_count":
                    sources[-1] = sources[-2]
                elif failure == "absolute_sum":
                    weights[0] += 1.
                    weights[2] -= 1.
                else:
                    samples[1].normalization_sum_weight = 0.
                baseline = {"analysis": {"cuts": [{"variable": "m_gg", "min": 123.}]}}
                with patch.dict(sys.modules, {"read_root_varfiles": reader,
                                              "ho_xgboost_analysis": trainer}):
                    with patch.object(self.subject, "load_config", return_value=baseline):
                        with self.assertRaises(ValueError):
                            self.run_workflow(samples=samples)
                self.assertEqual(reader.read_ho_ROOT_varfile.call_count, 2)
                trainer.run_ho_signal_background_analysis.assert_not_called()

    def test_baseline_card_contributes_only_its_cuts(self):
        samples = self.samples()
        reader, trainer, _ = self.modules(samples)
        baseline_path = Path("comparison.yaml")
        analysis = {"name": "actual", "run_tag": RUN_TAG, "luminosity_fb": 10.,
                    "samples": [sample.name for sample in samples],
                    "_resolved_detector_response": "ssc",
                    "xgboost": {"baseline_cuts_config": str(baseline_path), "seed": 42}}
        baseline = {"analysis": {
            "name": "other", "run_tag": "other-run", "luminosity_fb": 999.,
            "analysis_root": "/unused", "output_dir": "/unused", "samples": ["other"],
            "rate_factors": {"Signal": 200.}, "max_events": 1,
            "detector_response": "none", "xgboost": {"seed": 0},
            "cuts": [{"variable": "m_gg", "min": 123., "max": 127.}],
        }}
        with patch.dict(sys.modules, {"read_root_varfiles": reader, "ho_xgboost_analysis": trainer}):
            with patch.object(self.subject, "load_config", return_value=baseline) as loader:
                with patch.object(self.subject, "discover_samples") as discover:
                    self.run_workflow(analysis=analysis, samples=samples)
        loader.assert_called_once_with(baseline_path)
        discover.assert_not_called()
        args, kwargs = trainer.run_ho_signal_background_analysis.call_args
        self.assertEqual([item[0] for item in args[0]], samples)
        self.assertEqual(kwargs["baseline_cuts"], [self.subject.Cut("m_gg", 123., 127.)])
        self.assertEqual(kwargs["config"], analysis["xgboost"])
        self.assertIs(kwargs["metadata"]["resolved_analysis"], analysis)
        self.assertEqual(kwargs["metadata"]["run_tag"], RUN_TAG)
        self.assertEqual(kwargs["metadata"]["luminosity_fb"], 10.)
        self.assertEqual(kwargs["metadata"]["detector_response"], "ssc")
        self.assertEqual(kwargs["output_dir"], Path("out"))


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
