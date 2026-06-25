import array
import io
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import ROOT

import read_root_varfiles


class FakeReportAxes:
    def __init__(self) -> None:
        self.transAxes = object()
        self.fill_between_calls = []
        self.texts = []

    def fill_between(self, edges, lower, upper, **kwargs) -> None:
        self.fill_between_calls.append(
            {
                "edges": list(edges),
                "lower": list(lower),
                "upper": list(upper),
                "kwargs": kwargs,
            }
        )

    def step(self, *args, **kwargs) -> None:
        pass

    def text(self, *args, **kwargs) -> None:
        self.texts.append(args[2])

    def set_xlim(self, *args, **kwargs) -> None:
        pass

    def set_xlabel(self, *args, **kwargs) -> None:
        pass

    def set_ylabel(self, *args, **kwargs) -> None:
        pass

    def minorticks_on(self) -> None:
        pass

    def tick_params(self, *args, **kwargs) -> None:
        pass

    def set_ylim(self, *args, **kwargs) -> None:
        pass

    def legend(self, *args, **kwargs) -> None:
        pass


class FakeReportFigure:
    def tight_layout(self) -> None:
        pass

    def savefig(self, *args, **kwargs) -> None:
        pass


class FakeReportPlot:
    def __init__(self, axes: FakeReportAxes) -> None:
        self.axes = axes

    def subplots(self, *args, **kwargs):
        return FakeReportFigure(), self.axes

    def close(self, *args, **kwargs) -> None:
        pass


def write_gamma_gamma_varfile(path: Path, rows: list[tuple[list[float], float]]) -> None:
    root_file = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("Data2", "Gamma Gamma Data Tree")
    variables = array.array("d", [0.0] * 10)
    eventweight = array.array("d", [0.0])
    tree.Branch("variables", variables, "variables[10]/D")
    tree.Branch("eventweight", eventweight, "eventweight[1]/D")
    for values, weight in rows:
        for index, value in enumerate(values):
            variables[index] = value
        eventweight[0] = weight
        tree.Fill()
    tree.Write()
    root_file.Close()


class GammaGammaReaderTests(unittest.TestCase):
    def test_reads_all_gamma_gamma_features_and_eventweight(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample_var.root"
            write_gamma_gamma_varfile(
                path,
                [
                    ([125.0, 60.0, 0.2, 50.0, -0.3, 1.7, 3.1, 20.0, 0.1, 2.0], 0.5),
                    ([124.0, 45.0, -1.2, 30.0, 0.4, 2.4, 2.9, 12.0, -0.5, 2.0], 1.5),
                ],
            )

            features, labels, weights = read_root_varfiles.read_ROOT_varfile(path, sample_id=1, xsec=2.0)

        self.assertEqual(
            read_root_varfiles.FEATURE_NAMES,
            [
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
            ],
        )
        self.assertEqual(labels, [1, 1])
        self.assertEqual(weights, [1.0, 3.0])
        self.assertEqual(features[0], [125.0, 60.0, 0.2, 50.0, -0.3, 1.7, 3.1, 20.0, 0.1, 2.0])
        self.assertEqual(features[1], [124.0, 45.0, -1.2, 30.0, 0.4, 2.4, 2.9, 12.0, -0.5, 2.0])


class CutFlowTests(unittest.TestCase):
    def test_progress_bar_writes_current_sample_and_count(self) -> None:
        from analyze_lo_varfiles import ProgressBar

        stream = io.StringIO()
        progress = ProgressBar(total=2, label="Analyzing samples", stream=stream, enabled=True, width=10)
        progress.update(1, "signal_gg_h_aa")
        progress.update(2, "bkg_prompt_aa")
        progress.finish()

        output = stream.getvalue()
        self.assertIn("Analyzing samples", output)
        self.assertIn("1/2", output)
        self.assertIn("2/2", output)
        self.assertIn("signal_gg_h_aa", output)
        self.assertIn("bkg_prompt_aa", output)
        self.assertTrue(output.endswith("\n"))

    def test_terminal_summary_reports_outputs_and_totals(self) -> None:
        from analyze_lo_varfiles import build_terminal_summary

        summary = build_terminal_summary(
            metadata={"name": "baseline_cuts", "mode": "cuts", "run_tag": "run_02", "luminosity_fb": 100.0},
            rows=[
                {
                    "sample": "signal_gg_h_aa",
                    "category": "Signal",
                    "selected_cross_section_pb": 0.4,
                    "expected_events": 40000.0,
                    "mc_events_after_analysis": 10,
                },
                {
                    "sample": "bkg_prompt_aa",
                    "category": "Backgrounds",
                    "selected_cross_section_pb": 90.0,
                    "expected_events": 9000000.0,
                    "mc_events_after_analysis": 1,
                },
            ],
            output_dir=Path("/tmp/out"),
        )

        self.assertIn("Analysis summary", summary)
        self.assertIn("run tag: run_02", summary)
        self.assertIn("Signal", summary)
        self.assertIn("Backgrounds", summary)
        self.assertIn("expected events=40000", summary)
        self.assertIn("approx_significance_s_over_sqrt_b=13.3333", summary)
        self.assertIn("/tmp/out/summary.csv", summary)

    def test_analysis_totals_include_s_over_sqrt_b(self) -> None:
        from analyze_lo_varfiles import compute_analysis_totals

        totals = compute_analysis_totals(
            [
                {
                    "category": "Signal",
                    "selected_cross_section_pb": 0.25,
                    "expected_events": 25.0,
                    "mc_events_after_analysis": 3,
                },
                {
                    "category": "Backgrounds",
                    "selected_cross_section_pb": 1.0,
                    "expected_events": 100.0,
                    "mc_events_after_analysis": 7,
                },
            ]
        )

        self.assertEqual(totals["signal_expected_events"], 25.0)
        self.assertEqual(totals["background_expected_events"], 100.0)
        self.assertEqual(totals["approx_significance_s_over_sqrt_b"], 2.5)

    def test_zero_background_significance_is_none(self) -> None:
        from analyze_lo_varfiles import compute_analysis_totals

        totals = compute_analysis_totals(
            [
                {"category": "Signal", "expected_events": 25.0},
                {"category": "Backgrounds", "expected_events": 0.0},
            ]
        )

        self.assertIsNone(totals["approx_significance_s_over_sqrt_b"])

    def test_cli_run_tag_override_is_parsed(self) -> None:
        from analyze_lo_varfiles import build_parser

        args = build_parser().parse_args(["cuts", "--config", "card.yaml", "--run-tag", "run_02"])

        self.assertEqual(args.run_tag, "run_02")

    def test_report_signal_scale_cli_defaults_to_one_and_accepts_override(self) -> None:
        from hgammagamma.make_gammagamma_report import build_parser

        default_args = build_parser().parse_args([])
        scaled_args = build_parser().parse_args(["--signal-scale", "25"])

        self.assertEqual(default_args.signal_scale, 1.0)
        self.assertEqual(scaled_args.signal_scale, 25.0)

    def test_resolve_run_tag_prefers_cli_then_yaml_then_environment(self) -> None:
        from analyze_lo_varfiles import resolve_run_tag

        self.assertEqual(resolve_run_tag({"run_tag": "run_01"}, "run_02"), "run_02")
        self.assertEqual(resolve_run_tag({"run_tag": "run_01"}, None), "run_01")

    def test_rectangular_cuts_are_inclusive_and_anded(self) -> None:
        from analyze_lo_varfiles import Cut, apply_cuts

        rows = [
            {"m_gg": 120.0, "pt_gamma1": 40.0, "n_selected_photons": 2.0},
            {"m_gg": 125.0, "pt_gamma1": 39.9, "n_selected_photons": 2.0},
            {"m_gg": 130.0, "pt_gamma1": 50.0, "n_selected_photons": 1.0},
            {"m_gg": 131.0, "pt_gamma1": 60.0, "n_selected_photons": 2.0},
        ]
        cuts = [
            Cut(variable="m_gg", minimum=120.0, maximum=130.0),
            Cut(variable="pt_gamma1", minimum=40.0, maximum=None),
            Cut(variable="n_selected_photons", minimum=2.0, maximum=None),
        ]

        self.assertEqual(apply_cuts(rows, cuts), [True, False, False, False])

    def test_invalid_cut_variable_is_rejected(self) -> None:
        from analyze_lo_varfiles import Cut

        with self.assertRaises(ValueError):
            Cut(variable="not_a_variable", minimum=0.0, maximum=1.0)

    def test_nan_or_sentinel_values_do_not_pass_physics_cuts(self) -> None:
        from analyze_lo_varfiles import Cut, apply_cuts

        rows = [
            {"m_gg": -999.0},
            {"m_gg": math.nan},
            {"m_gg": 125.0},
        ]

        self.assertEqual(apply_cuts(rows, [Cut(variable="m_gg", minimum=120.0, maximum=130.0)]), [False, False, True])

    def test_summary_cross_section_includes_weight_scale(self) -> None:
        from analyze_lo_varfiles import Cut, SampleInfo, summarize_sample

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample_var.root"
            write_gamma_gamma_varfile(
                path,
                [
                    ([125.0, 60.0, 0.2, 50.0, -0.3, 1.7, 3.1, 20.0, 0.1, 2.0], 1.0),
                    ([110.0, 45.0, -1.2, 30.0, 0.4, 2.4, 2.9, 12.0, -0.5, 2.0], 1.0),
                ],
            )
            dat_file = Path(tmpdir) / "sample.dat"
            dat_file.write_text("")
            sample = SampleInfo(
                name="signal",
                category="Signal",
                sample_dir=Path(tmpdir),
                var_file=path,
                dat_file=dat_file,
                cross_section_pb=100.0,
                cross_section_error_pb=None,
                weight_scale=0.01,
                events_read=2.0,
                sum_weight=2.0,
            )

            row = summarize_sample(sample, [Cut(variable="m_gg", minimum=120.0, maximum=130.0)], luminosity_fb=10.0)

        self.assertEqual(row["raw_cross_section_pb"], 100.0)
        self.assertEqual(row["weight_scale"], 0.01)
        self.assertEqual(row["cross_section_pb"], 1.0)
        self.assertEqual(row["selected_cross_section_pb"], 0.5)
        self.assertEqual(row["expected_events"], 5000.0)

    def test_legacy_signal_uses_default_branching_and_k_factor(self) -> None:
        from analyze_lo_varfiles import SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT, resolve_weight_scale

        weight_scale = resolve_weight_scale("signal_gg_h_aa", "Signal", dat_weight_scale=1.0, configured_rate_factors=None)

        self.assertEqual(weight_scale, SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT)

    def test_configured_rate_factors_override_defaults_and_dat_values(self) -> None:
        from analyze_lo_varfiles import resolve_weight_scale

        configured = {
            "Signal": 0.01,
            "Backgrounds": 1.4,
            "signal_gg_h_aa": 0.00454,
            "bkg_prompt_aa": 1.2,
        }

        self.assertEqual(resolve_weight_scale("signal_gg_h_aa", "Signal", 1.0, configured), 0.00454)
        self.assertEqual(resolve_weight_scale("new_signal", "Signal", 1.0, configured), 0.01)
        self.assertEqual(resolve_weight_scale("bkg_prompt_aa", "Backgrounds", 1.0, configured), 1.2)
        self.assertEqual(resolve_weight_scale("bkg_gamma_j", "Backgrounds", 1.0, configured), 1.4)
        self.assertEqual(resolve_weight_scale("unconfigured_bkg", "Backgrounds", 2.0, {}), 2.0)

    def test_xgboost_summary_cross_section_includes_rate_factor(self) -> None:
        from xgboost_root_varfiles_module import _summarize_full_sample

        class FakeModel:
            def predict_proba(self, rows):
                import numpy as np

                return np.asarray([[0.2, 0.8], [0.7, 0.3]])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "sample_var.root"
            write_gamma_gamma_varfile(
                path,
                [
                    ([125.0, 60.0, 0.2, 50.0, -0.3, 1.7, 3.1, 20.0, 0.1, 2.0], 1.0),
                    ([124.0, 45.0, -1.2, 30.0, 0.4, 2.4, 2.9, 12.0, -0.5, 2.0], 1.0),
                ],
            )
            sample = {
                "sample": "signal",
                "category": "Signal",
                "input_file": str(path),
                "xsec_fb": 1000.0,
                "rate_factor": 0.5,
                "effective_xsec_fb": 500.0,
                "normalisation_weight": 2.0,
            }

            row = _summarize_full_sample(FakeModel(), sample, label=1, threshold=0.5, luminosity=10.0)

        self.assertEqual(row["raw_cross_section_pb"], 1.0)
        self.assertEqual(row["weight_scale"], 0.5)
        self.assertEqual(row["cross_section_pb"], 0.5)
        self.assertEqual(row["selected_cross_section_pb"], 0.25)
        self.assertEqual(row["expected_events"], 2500.0)

    def make_report_sample(
        self,
        name: str,
        category: str,
        selected_xsec: float,
        y_values: list[float],
    ):
        from hgammagamma import make_gammagamma_report as report

        histogram = report.Histogram(
            title="diphoton invariant mass",
            x_min=0.0,
            x_max=2.0,
            x=[0.5, 1.5],
            y=y_values,
        )
        return report.SampleResult(
            name=name,
            label=name,
            category=category,
            sample_dir=Path("/tmp"),
            top_file=Path("/tmp/sample.top"),
            dat_file=Path("/tmp/sample.dat"),
            cross_section_pb=selected_xsec,
            cross_section_error_pb=None,
            events_read=10.0,
            selected_events=10.0,
            sum_weight=10.0,
            sum_diphoton_weight=10.0,
            weight_scale=1.0,
            efficiency=1.0,
            selected_cross_section_pb=selected_xsec,
            histograms={"diphoton invariant mass": histogram},
        )

    def capture_report_plot(self, *, density: bool, signal_scale: float = 1.0) -> FakeReportAxes:
        from hgammagamma import make_gammagamma_report as report

        axes = FakeReportAxes()
        samples = [
            self.make_report_sample("background", "Backgrounds", 10.0, [3.0, 7.0]),
            self.make_report_sample("signal", "Signal", 20.0, [5.0, 5.0]),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.object(report, "ensure_matplotlib", return_value=FakeReportPlot(axes)):
                report.plot_histogram(
                    "diphoton invariant mass",
                    samples,
                    Path(tmpdir),
                    "selected_xsec",
                    density=density,
                    signal_scale=signal_scale,
                )
        return axes

    def test_density_report_plot_overlays_transparent_histograms(self) -> None:
        axes = self.capture_report_plot(density=True)

        self.assertEqual(len(axes.fill_between_calls), 2)
        signal_call = axes.fill_between_calls[1]
        self.assertEqual(signal_call["lower"], [0.0, 0.0, 0.0])
        self.assertEqual(signal_call["upper"], [0.5, 0.5, 0.5])
        self.assertLess(signal_call["kwargs"]["alpha"], 0.6)
        self.assertIn("Norm: unit area", axes.texts)

    def test_non_density_report_plot_stacks_and_labels_scaled_signal(self) -> None:
        axes = self.capture_report_plot(density=False, signal_scale=5.0)

        self.assertEqual(len(axes.fill_between_calls), 2)
        signal_call = axes.fill_between_calls[1]
        self.assertEqual(signal_call["lower"], [3.0, 7.0, 7.0])
        self.assertEqual(signal_call["upper"], [53.0, 57.0, 57.0])
        self.assertEqual(signal_call["kwargs"]["label"], "signal x5")


if __name__ == "__main__":
    unittest.main()
