"""HO rate/report checks; these do not need ROOT, Herwig, or an ihixs build."""

import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hgammagamma import make_gammagamma_report as report


class HONormalizationReportTests(unittest.TestCase):
    def sample(self, selected_pb=0.227, normalization_kind="ihixs_n3lo"):
        return report.SampleResult(
            name="signal_gg_h_aa", label="Higgs", category="Signal",
            analysis_name="SSC_GEM_weighted_response", detector_response="ssc",
            response_mode="genuine", weighted_hypotheses=True, metadata_inferred=False,
            sample_dir=Path("sample"), top_file=Path("sample.top"), dat_file=Path("sample.dat"),
            cross_section_pb=200.0, cross_section_error_pb=0.02,
            events_read=10., selected_events=5., sum_weight=0.0227,
            sum_diphoton_weight=0.01135, weight_scale=0.00227, efficiency=0.5,
            selected_cross_section_pb=selected_pb, histograms={},
            normalization_kind=normalization_kind, native_cross_section_pb=190.,
            ihixs_record_sha256="a" * 64,
            rate_uncertainties={"scale_up_pb": 5., "scale_down_pb": 7., "pdf_stddev_pb": 3.},
        )

    def test_yield_has_pb_to_fb_conversion_and_br_once(self):
        sample = self.sample()
        self.assertEqual(report.expected_events(sample, 1.), 227.)
        self.assertEqual(report.expected_events(sample, 10.), 2270.)

    def test_signed_selected_yields_remain_signed(self):
        self.assertEqual(report.expected_events(self.sample(selected_pb=-0.1), 2.), -200.)

    def test_loaded_signed_efficiency_applies_production_rate_and_br_once(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            directory = root / "Signal/events/signal_gg_h_aa"
            directory.mkdir(parents=True)
            dat_file = directory / "signal-test.dat"
            dat_file.write_text("events_read 3\nevents_with_two_selected_photons 2\n"
                                "sum_weight 0.03632\nsum_diphoton_weight -0.00454\n"
                                "weight_scale 0.00227\n")
            dat_file.with_suffix(".top").write_text("")
            sidecar = {"native_cross_section_pb": 190., "ihixs_record_sha256": "a" * 64,
                       "ihixs": {"uncertainties": {"scale_up_pb": 5.}}}
            args = SimpleNamespace(analysis_root=root, run_tag="test", samples="all",
                                   allow_missing_xsec=False)
            with patch.object(report, "load_ho_signal_sidecar", return_value=sidecar), \
                    patch.object(report, "validate_ho_analysis_summary"), \
                    patch.object(report, "parse_cross_section", return_value=(200., .01)):
                sample = report.load_samples(args)[0]
            self.assertAlmostEqual(sample.efficiency, -.125)
            self.assertAlmostEqual(sample.selected_cross_section_pb, -.05675)
            self.assertAlmostEqual(report.expected_events(sample, 10.), -567.5)

    def test_luminosity_rejects_nonfinite_and_nonpositive_values(self):
        for value in (0., -1., float("nan"), float("inf")):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    report.expected_events(self.sample(), value)

    def test_summary_separates_physical_yield_and_inclusive_rate_uncertainties(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = report.write_summary_csv([self.sample()], path, luminosity_fb=10.)
            with (path / summary).open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(float(rows[0]["expected_events"]), 2270.)
            self.assertEqual(float(rows[0]["cross_section_pb"]), 200.)
            self.assertEqual(float(rows[0]["native_cross_section_pb"]), 190.)
            self.assertEqual(rows[0]["normalization_kind"], "ihixs_n3lo")
            self.assertEqual(json.loads(rows[0]["rate_uncertainties_pb"])["pdf_stddev_pb"], 3.)

    def test_default_lo_summary_does_not_add_ho_or_yield_columns(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp)
            summary = report.write_summary_csv([self.sample(normalization_kind="generator")], path)
            with (path / summary).open() as stream:
                columns = csv.DictReader(stream).fieldnames
            self.assertNotIn("expected_events", columns)
            self.assertNotIn("normalization_kind", columns)

    def test_ho_signal_cannot_fall_back_to_native_sidecar(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp) / "signal_gg_h_aa"
            directory.mkdir()
            (directory / "campaign.json").write_text(json.dumps({
                "sample": directory.name, "run_tag": "test", "matching": "POWHEG",
                "hard_accuracy": "NNLO+PS (HJMiNNLO)", "process": "HJMiNNLO; h -> gamma gamma",
                "ebeam_gev": 20000., "pdf": "NNPDF40_nnlo_as_01180_qed",
                "lhaid": 336100, "weight_scale": 0.00227,
            }))
            with self.assertRaisesRegex(ValueError, "requires ihixs"):
                report.parse_cross_section(directory, "test")
            (directory / "normalization-test.json").write_text(json.dumps({
                "run_tag": "test", "sample": directory.name, "cross_section_pb": 190.,
            }))
            with self.assertRaises(ValueError):
                report.parse_cross_section(directory, "test")

    def test_same_sample_name_without_ho_manifest_keeps_lo_normalization(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp) / "signal_gg_h_aa"
            directory.mkdir()
            (directory / "normalization-test.json").write_text(json.dumps({
                "run_tag": "test", "sample": directory.name,
                "cross_section_pb": 50., "cross_section_error_pb": 1.,
            }))
            self.assertEqual(report.parse_cross_section(directory, "test"), (50., 1.))

    def test_ihixs_sidecar_cannot_be_reused_for_lo_or_backgrounds(self):
        with TemporaryDirectory() as tmp:
            for name in ("signal_gg_h_aa", "bkg_prompt_aa"):
                directory = Path(tmp) / name
                directory.mkdir()
                (directory / "normalization-test.json").write_text(json.dumps({
                    "normalization_kind": "ihixs_n3lo", "run_tag": "test",
                    "sample": name, "cross_section_pb": 200.,
                }))
                with self.assertRaisesRegex(ValueError, "compatible HO campaign"):
                    report.parse_cross_section(directory, "test")

    def test_replaced_analysis_summary_is_rejected(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            dat_file = directory / "signal-test.dat"
            dat_file.write_text("# HwSimPostAnalysis_gammagamma_SSC summary\n"
                                "events_read 10\nsum_weight 1\nsum_diphoton_weight 0.5\n")
            (directory / "campaign.json").write_text(json.dumps({"analysis": {
                "#": "HwSimPostAnalysis_gammagamma_SSC summary",
                "events_read": "10", "sum_weight": "1", "sum_diphoton_weight": "0.5",
            }}))
            report.validate_ho_analysis_summary(directory, dat_file)
            dat_file.write_text("events_read 10\nsum_weight 1\nsum_diphoton_weight 0.9\n")
            with self.assertRaisesRegex(ValueError, "differs from its completed campaign"):
                report.validate_ho_analysis_summary(directory, dat_file)


if __name__ == "__main__":
    unittest.main()
