"""Real ROOT/XGBoost end-to-end checks; optional analysis packages may be absent."""

import array
import builtins
import csv
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import ho_xgboost_analysis as subject

try:
    import ROOT
    import xgboost as xgb
    import sklearn.metrics
    import matplotlib.pyplot
    import analyze_lo_varfiles as cli
    import read_root_varfiles as reader
except (ImportError, RuntimeError) as exc:
    OPTIONAL_REASON = str(exc)
else:
    OPTIONAL_REASON = ""


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path):
    def reject_constant(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(Path(path).read_text(), parse_constant=reject_constant)


def write_sample(directory, name, sample_index, source_count):
    """Make signed sources with three exclusive detector-response outcomes."""
    directory.mkdir(parents=True)
    path = directory / "sample-test_var.root"
    output = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("Data2", "Synthetic HO detector outcomes")
    variables = array.array("d", [0.] * 10)
    eventweight = array.array("d", [0.])
    sourceevent = array.array("q", [0])
    tree.Branch("variables", variables, "variables[10]/D")
    tree.Branch("eventweight", eventweight, "eventweight[1]/D")
    tree.Branch("sourceevent", sourceevent, "sourceevent[1]/L")
    signal = sample_index == 0
    scale = .00227 if signal else 1.
    rows, weights, sources, entries = [], [], [], []
    for source in range(source_count):
        signed = -1. if source % 17 == 0 else 1.
        mass = 125. + .15 * (source % 11 - 5) if signal else 110. + source % 41
        features = [mass, 60. + source * .04, .2, 45. + source * .02,
                    -.3, 1.7, 3.1, 20. + source * .01, .1]
        for photons, probability in ((0, .4), (1, .2), (2, .4)):
            values = features + [2.] if photons == 2 else [-999.] * 9 + [float(photons)]
            if photons == 1:
                values[1], values[2] = features[1], features[2]
            for index, value in enumerate(values):
                variables[index] = value
            eventweight[0] = signed * scale * probability
            sourceevent[0] = source
            tree.Fill()
            rows.append(dict(zip(reader.FEATURE_NAMES, values)))
            weights.append(eventweight[0])
            sources.append(source)
            entries.append(len(entries))
    tree.Write()
    output.Close()
    total = math.fsum(weights)
    # Eleven consumed but discarded positive source events remain in the denominator.
    denominator = total + 11. * scale
    dat = directory / "sample-test.dat"
    dat.write_text(f"sum_weight = {total}\nsum_tree_weight = {total}\n")
    (directory / "campaign.json").write_text(json.dumps({"sample": name, "source_events": source_count + 11}))
    sample = cli.SampleInfo(
        name=name, category="Signal" if signal else "Backgrounds", sample_dir=directory,
        var_file=path, dat_file=dat, cross_section_pb=30. if signal else 1. + sample_index,
        cross_section_error_pb=.01, weight_scale=scale, events_read=float(source_count),
        sum_weight=total, analysis_name="SSC_GEM_weighted_response", detector_response="ssc",
        response_mode=("genuine", "genuine", "jet_fake", "electron_fake")[sample_index],
        normalization_sum_weight=denominator,
        normalization_kind="ihixs_n3lo" if signal else "generator",
        ihixs_record_sha256="a" * 64 if signal else "", requires_full_sample=True,
        shower_quality={"termination": "source_exhausted", "attempted_events": source_count + 11,
                        "saved_events": source_count, "discarded_events": 11},
        tree_entries=len(rows), sum_tree_weight=total, sum_abs_weight=math.fsum(abs(w) for w in weights),
    )
    return sample, rows, weights, sources, entries


@unittest.skipIf(bool(OPTIONAL_REASON), f"optional ROOT/ML dependencies unavailable: {OPTIONAL_REASON}")
class HOXGBoostEndToEndTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="ho-ml-test-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.root = Path(cls.temporary.name)
        names = ("signal_gg_h_aa", "bkg_prompt_aa", "bkg_gamma_j", "bkg_dy_ee")
        cls.source_counts = dict(zip(names, (253, 257, 263, 269)))
        cls.expected_split_counts = dict(zip(names, (
            {"train": 152, "validation": 51, "test": 50},
            {"train": 154, "validation": 52, "test": 51},
            {"train": 158, "validation": 53, "test": 52},
            {"train": 161, "validation": 54, "test": 54},
        )))
        cls.loaded = [write_sample(cls.root / name, name, index, cls.source_counts[name])
                      for index, name in enumerate(names)]
        cls.samples = [loaded[0] for loaded in cls.loaded]
        cls.analysis = {
            "_resolved_detector_response": "ssc",
            "xgboost": {"seed": 12345, "validation_size": .2, "test_size": .2,
                        "model_params": {"n_estimators": 5, "max_depth": 2, "n_jobs": 1}},
        }
        cls.luminosity = 10.
        cls.output = cls.root / "complete"
        cls.fit_calls, cls.scan_calls = [], []
        real_fit, real_scan = xgb.XGBClassifier.fit, subject.scan_thresholds

        def record_fit(model, features, labels, **kwargs):
            cls.fit_calls.append((features.copy(), labels.copy(), kwargs["sample_weight"].copy()))
            return real_fit(model, features, labels, **kwargs)

        def record_scan(*args, **kwargs):
            cls.scan_calls.append((tuple(value.copy() for value in args), kwargs.copy()))
            return real_scan(*args, **kwargs)

        with mock.patch.object(xgb.XGBClassifier, "fit", record_fit), \
                mock.patch.object(subject, "scan_thresholds", record_scan), \
                mock.patch.object(reader, "read_ho_ROOT_varfile", wraps=reader.read_ho_ROOT_varfile) as strict_reader:
            cls.result = cli.run_ho_xgboost(cls.analysis, "test_ml", "test", cls.luminosity,
                                          cls.samples, cls.output, progress_enabled=False)
            cls.reader_calls = strict_reader.call_args_list
        cls.partitions = read_csv(cls.output / "partitions.csv")
        cls.scores = read_csv(cls.output / "scores.csv")
        cls.metrics = read_json(cls.output / "metrics.json")
        cls.comparison = read_json(cls.output / "comparison.json")

    def test_strict_reader_and_complete_artifact_contract(self):
        self.assertEqual(len(self.reader_calls), 4)
        self.assertEqual({call.args[0] for call in self.reader_calls}, {sample.var_file for sample in self.samples})
        self.assertEqual(self.metrics["status"], "complete")
        self.assertIsNotNone(self.metrics["best_threshold"])
        expected = {"index.html", "summary.csv", "summary.json", "comparison.csv", "comparison.json",
                    "metrics.json", "signal_background_xgboost.json", "model_metadata.json", "scores.csv",
                    "partitions.csv", "threshold_scan.csv", "roc.png", "feature_importance.png",
                    "score_distributions.png", "mass_distributions.png"}
        self.assertEqual({path.name for path in self.output.iterdir()}, expected)
        for path in self.output.iterdir():
            self.assertGreater(path.stat().st_size, 0)
            if path.suffix == ".json":
                read_json(path)
            elif path.suffix == ".png":
                self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
        self.assertEqual(len(read_csv(self.output / "threshold_scan.csv")), 501)
        self.assertEqual(len(read_csv(self.output / "summary.csv")), 4)

    def test_model_reload_matches_saved_scores_and_nine_kinematic_features(self):
        self.assertEqual(self.metrics["feature_names"], reader.FEATURE_NAMES[:9])
        self.assertNotIn("n_selected_photons", self.metrics["feature_names"])
        self.assertNotIn("sourceevent", self.metrics["feature_names"])
        model = xgb.XGBClassifier()
        model.load_model(self.output / "signal_background_xgboost.json")
        self.assertEqual(model.n_features_in_, 9)
        input_rows = {sample.name: rows for sample, rows, _, _, _ in self.loaded}
        features = np.asarray([[input_rows[row["sample"]][int(row["tree_entry"])][name]
                               for name in self.metrics["feature_names"]] for row in self.scores])
        np.testing.assert_allclose(model.predict_proba(features)[:, 1],
                                   [float(row["score"]) for row in self.scores], rtol=0., atol=1.e-7)

    def test_all_outcomes_share_group_and_fractions_use_full_population(self):
        self.assertEqual(len(self.partitions), sum(self.source_counts.values()) * 3)
        self.assertEqual(len(self.scores), sum(self.source_counts.values()))
        groups = {}
        for row in self.partitions:
            key = row["sample"], row["sourceevent"]
            groups.setdefault(key, []).append(row)
        for (sample, _), rows in groups.items():
            self.assertEqual(len(rows), 3)
            self.assertEqual(len({row["partition"] for row in rows}), 1)
            self.assertEqual(sum(row["eligible"] == "True" for row in rows), 1)
        for sample in self.samples:
            self.assertEqual(self.metrics["sample_group_counts"][sample.name],
                             self.expected_split_counts[sample.name])
            for partition, count in self.metrics["sample_group_counts"][sample.name].items():
                self.assertAlmostEqual(self.metrics["sample_split_fractions"][sample.name][partition],
                                       count / self.source_counts[sample.name])
        self.assertTrue(any(float(row["physical_weight"]) < 0 for row in self.scores))

    def test_only_training_rows_are_fit_and_only_validation_rows_choose_threshold(self):
        self.assertEqual(len(self.fit_calls), 1)
        self.assertEqual(len(self.scan_calls), 1)
        features, labels, train_weights = self.fit_calls[0]
        records = {(row["sample"], int(row["tree_entry"])): row for row in self.partitions}
        expected_features, expected_labels, expected_sources, expected_names, expected_weights = [], [], [], [], []
        expected_training_physical = []
        for sample, rows, weights, sources, entries in self.loaded:
            factor = self.luminosity * 1000 * sample.cross_section_pb * sample.weight_scale / sample.normalization_sum_weight
            fractions = {partition: count / self.source_counts[sample.name]
                         for partition, count in self.expected_split_counts[sample.name].items()}
            for row, weight, source, entry in zip(rows, weights, sources, entries):
                partition = records[sample.name, entry]["partition"]
                if row["n_selected_photons"] < 2:
                    continue
                if partition == "train":
                    expected_features.append([row[name] for name in subject.FEATURE_NAMES])
                    expected_labels.append(int(sample.category == "Signal"))
                    expected_training_physical.append(weight * factor / fractions["train"])
                elif partition == "validation":
                    expected_sources.append(source)
                    expected_names.append(sample.name)
                    expected_weights.append(weight * factor / fractions["validation"])
        np.testing.assert_array_equal(features, expected_features)
        np.testing.assert_array_equal(labels, expected_labels)
        self.assertTrue(np.all(train_weights >= 0))
        self.assertAlmostEqual(train_weights[labels == 1].sum(), train_weights[labels == 0].sum())
        # Different sample sizes round to different actual training fractions;
        # their physical mixture must be projected before class balancing.
        expected_absolute = np.abs(expected_training_physical)
        expected_balanced = np.empty_like(expected_absolute)
        for label in (0, 1):
            mask = labels == label
            expected_balanced[mask] = expected_absolute[mask] / expected_absolute[mask].sum() * len(labels) / 2
        np.testing.assert_allclose(train_weights, expected_balanced, rtol=1.e-12)
        scan_args, _ = self.scan_calls[0]
        np.testing.assert_array_equal(scan_args[3], expected_names)
        np.testing.assert_array_equal(scan_args[4], expected_sources)
        np.testing.assert_allclose(scan_args[2], expected_weights)
        self.assertEqual(len(scan_args[0]), sum(counts["validation"] for counts in self.expected_split_counts.values()))

    def test_cut_and_ml_yields_use_same_heldout_events_and_consumed_source_denominator(self):
        self.assertEqual(self.comparison["evaluation_partition"], "test")
        threshold = self.metrics["best_threshold"]
        score_by_entry = {(row["sample"], int(row["tree_entry"])): row for row in self.scores}
        for selection in ("baseline", "xgboost"):
            summaries = {row["sample"]: row for row in self.comparison[selection]["samples"]}
            for sample, rows, weights, sources, entries in self.loaded:
                accepted = []
                for row, weight, source, entry in zip(rows, weights, sources, entries):
                    score = score_by_entry.get((sample.name, entry))
                    if score is None or score["partition"] != "test":
                        continue
                    selected = (120. <= row["m_gg"] <= 130. if selection == "baseline"
                                else float(score["score"]) >= threshold)
                    if selected:
                        accepted.append((source, weight))
                summary = summaries[sample.name]
                test_count = self.expected_split_counts[sample.name]["test"]
                fraction = test_count / self.source_counts[sample.name]
                factor = self.luminosity * 1000 * sample.cross_section_pb * sample.weight_scale / sample.normalization_sum_weight
                projected = [weight * factor / fraction for _, weight in accepted]
                self.assertAlmostEqual(summary["expected_events"], math.fsum(projected))
                self.assertAlmostEqual(summary["mc_variance_events"], math.fsum(weight * weight for weight in projected))
                self.assertEqual(summary["selected_entries"], len(accepted))
                self.assertEqual(summary["evaluation_partition"], "test")
                self.assertEqual(summary["partition_entries"], test_count * 3)
                self.assertAlmostEqual(summary["partition_fraction"], fraction)
                self.assertGreater(summary["normalization_sum_weight"], summary["sum_weight"])
                self.assertAlmostEqual(summary["analysis_efficiency"],
                                       math.fsum(weight for _, weight in accepted) / fraction / sample.normalization_sum_weight)

    def test_insufficient_statistics_has_null_threshold_and_no_ml_rows(self):
        output = self.root / "insufficient"
        analysis = {**self.analysis, "xgboost": {**self.analysis["xgboost"],
                                                 "min_background_effective_count": 10000.}}
        result = cli.run_ho_xgboost(analysis, "insufficient", "test", self.luminosity,
                                   self.samples, output, progress_enabled=False)
        metrics = read_json(output / "metrics.json")
        self.assertEqual(metrics["status"], "insufficient_statistics")
        self.assertIsNone(metrics["best_threshold"])
        self.assertIsNone(metrics["totals"])
        self.assertEqual(result.rows, [])
        self.assertEqual(read_csv(output / "summary.csv"), [])
        comparison = read_json(output / "comparison.json")
        self.assertIsNone(comparison["xgboost"])
        self.assertEqual(len(comparison["baseline"]["samples"]), 4)
        self.assertEqual({row["selection"] for row in read_csv(output / "comparison.csv")}, {"baseline"})
        self.assertIn("insufficient_statistics", (output / "index.html").read_text())
        self.assertEqual(len(read_csv(output / "threshold_scan.csv")), 501)
        for path in output.glob("*.json"):
            read_json(path)

    def test_selected_sentinel_features_are_rejected(self):
        sample, rows, weights, sources, entries = self.loaded[0]
        invalid_rows = [dict(row) for row in rows]
        invalid_rows[0]["n_selected_photons"] = 2.
        with self.assertRaisesRegex(ValueError, "invalid selected diphoton features"):
            subject._dataset([(sample, invalid_rows, weights, sources, entries)], self.luminosity)


class HOXGBoostDependencyTests(unittest.TestCase):
    def test_invalid_training_configuration_fails_before_fitting(self):
        for config in ({"validation_size": 0.}, {"test_size": .9, "validation_size": .2},
                       {"min_background_effective_count": 0.}, {"systematics": -1.},
                       {"systematics": float("nan")}, {"seed": True},
                       {"model_params": {"objective": "multi:softprob"}},
                       {"model_params": {"early_stopping_rounds": 5}}):
            with self.subTest(config=config), self.assertRaises(ValueError):
                subject.validate_config(config)

    def test_missing_optional_ml_dependency_has_actionable_error(self):
        original_import = builtins.__import__

        def without_xgboost(name, *args, **kwargs):
            if name == "xgboost":
                raise ImportError("test dependency unavailable")
            return original_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=without_xgboost):
            with self.assertRaisesRegex(RuntimeError, "requires xgboost, scikit-learn, and matplotlib"):
                subject._dependencies()

    def test_publish_failure_restores_old_report_and_preserves_unrelated_files(self):
        with tempfile.TemporaryDirectory(prefix="ho-ml-publish-") as temporary:
            root = Path(temporary)
            stage, output = root / "stage", root / "output"
            stage.mkdir()
            output.mkdir()
            previous = {"index.html": "old report", "metrics.json": '{"old": true}',
                        "notes.txt": "user notes"}
            for name, content in previous.items():
                (output / name).write_text(content)
            for name, content in {"index.html": "new report", "metrics.json": '{"new": true}',
                                  "scores.csv": "new generated asset"}.items():
                (stage / name).write_text(content)
            real_replace = os.replace

            def fail_last_publication(source, destination):
                if Path(source) == stage / "index.html":
                    raise OSError("simulated mid-publication failure")
                return real_replace(source, destination)

            with mock.patch.object(subject.os, "replace", side_effect=fail_last_publication):
                with self.assertRaisesRegex(OSError, "mid-publication"):
                    subject._publish(stage, output)
            self.assertEqual({path.name: path.read_text() for path in output.iterdir()}, previous)


if __name__ == "__main__":
    unittest.main()
