import array
import math
import tempfile
import unittest
from pathlib import Path

import ROOT

import read_root_varfiles


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


if __name__ == "__main__":
    unittest.main()
