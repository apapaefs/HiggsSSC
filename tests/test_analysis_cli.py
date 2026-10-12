"""Campaign-independent downstream commands; fixtures do not require PyROOT."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import analyze_lo_varfiles as analysis


REPO = Path(__file__).resolve().parents[1]
SAMPLES = ("signal_gg_h_aa", "bkg_prompt_aa", "bkg_gamma_j", "bkg_dy_ee")


class AnalysisCLITests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.campaign = self.root / "campaign"
        self.card = self.root / "cuts.yaml"
        self.card.write_text("""analysis:
  name: test
  run_tag: obsolete
  analysis_root: missing_old_campaign
  output_dir: missing_old_output
  detector_response: ssc
  luminosity_fb: 100
  samples: signal_gg_h_aa,bkg_prompt_aa,bkg_gamma_j,bkg_dy_ee
  cuts:
    - variable: m_gg
      min: 120
      max: 130
""")

    def add_sample(self, name):
        category = "Signal" if name.startswith("signal_") else "Backgrounds"
        folder = self.campaign / category / "events" / name
        folder.mkdir(parents=True)
        (folder / "events-current_var.root").write_bytes(b"mock tree")
        (folder / "events-current.dat").write_text(
            "analysis SSC_GEM_weighted_response\ndetector_response ssc\n"
            "response_mode genuine\nevents_read 2\nsum_weight 2\nweight_scale 1\n")
        banner = folder / "mg5_process/Events/current/current_tag_1_banner.txt"
        banner.parent.mkdir(parents=True)
        banner.write_text("Integrated weight (pb) : 10\n")

    def test_override_paths_tag_and_luminosity_use_selected_campaign_without_editing_card(self):
        for name in SAMPLES:
            self.add_sample(name)
        before = self.card.read_bytes()
        output = self.root / "selected_output"
        with patch.object(analysis, "read_named_ROOT_varfile", return_value=(
                [{"m_gg": 125.}, {"m_gg": 150.}], [1., 1.])):
            result = analysis.run_cuts(
                self.card, run_tag_override="current", progress_enabled=False,
                analysis_root_override=self.campaign, output_dir_override=output,
                luminosity_fb_override=25.)
        self.assertEqual(result.output_dir, output)
        self.assertEqual(result.metadata["run_tag"], "current")
        self.assertEqual(result.metadata["luminosity_fb"], 25.)
        self.assertEqual(len(result.rows), 4)
        self.assertEqual(self.card.read_bytes(), before)
        saved = json.loads((output / "summary.json").read_text())
        self.assertIn("metadata", saved)
        self.assertTrue(result.index_html.is_file())
        for row in result.rows:
            if row["category"] == "Backgrounds":
                self.assertAlmostEqual(row["expected_events"], 125000.)

    def test_incomplete_campaign_fails_before_tree_reads_or_outputs(self):
        self.add_sample("signal_gg_h_aa")
        output = self.root / "output"
        with patch.object(analysis, "read_named_ROOT_varfile") as reader:
            with self.assertRaisesRegex(FileNotFoundError, "missing requested analyzed samples"):
                analysis.run_cuts(self.card, "current", False,
                                  analysis_root_override=self.campaign, output_dir_override=output)
        reader.assert_not_called()
        self.assertFalse(output.exists())

    def test_invalid_luminosity_fails_before_sample_discovery(self):
        for value in (0., -1., float("nan"), float("inf")):
            with self.subTest(value=value), patch.object(analysis, "discover_samples") as discover:
                with self.assertRaisesRegex(ValueError, "finite and positive"):
                    analysis.load_analysis_inputs(self.card, luminosity_fb_override=value)
                discover.assert_not_called()

    def test_cli_options_forward_to_cuts_and_classifier(self):
        for action in ("cuts", "xgboost"):
            with self.subTest(action=action), patch.object(analysis, "run_" + action) as runner:
                with patch.object(analysis, "build_terminal_summary", return_value="done"):
                    analysis.main([action, "--config", str(self.card), "--run-tag", "current",
                                   "--analysis-root", str(self.campaign), "--output-dir", "out",
                                   "--luminosity-fb", "25", "--no-progress"])
                runner.assert_called_once_with(
                    self.card, run_tag_override="current", progress_enabled=False,
                    analysis_root_override=self.campaign, output_dir_override=Path("out"),
                    luminosity_fb_override=25.)

    def test_all_student_cards_require_the_complete_campaign_and_have_no_account_paths(self):
        for card in sorted((REPO / "hgammagamma/analysis_cards").glob("*.yaml")):
            with self.subTest(card=card.name):
                config = analysis.load_config(card)["analysis"]
                self.assertEqual(analysis.requested_samples(config["samples"]), set(SAMPLES))
                self.assertNotIn("analysis_root", config)
                self.assertNotIn("output_dir", config)
                self.assertNotIn("run_tag", config)

    def test_help_without_site_packages(self):
        result = subprocess.run([sys.executable, "-S", str(REPO / "analyze_lo_varfiles.py"),
                                 "cuts", "--help"], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--analysis-root", result.stdout)
        self.assertIn("--luminosity-fb", result.stdout)


if __name__ == "__main__":
    unittest.main()
