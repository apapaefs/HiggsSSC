import io
import unittest
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
from unittest.mock import patch

from hgammagamma import run_gammagamma_campaign as campaign


class GammaGammaCampaignRunCardTests(unittest.TestCase):
    def test_default_samples_include_ssc_fake_background_modes(self) -> None:
        samples = {sample.name: sample for sample in campaign.SAMPLES}

        self.assertEqual(samples["signal_gg_h_aa"].response_mode, "genuine")
        self.assertEqual(samples["bkg_prompt_aa"].response_mode, "genuine")
        self.assertEqual(samples["bkg_gamma_j"].process, "p p > a j")
        self.assertEqual(samples["bkg_gamma_j"].response_mode, "gammajet")
        self.assertEqual(samples["bkg_dy_ee"].process, "p p > e+ e-")
        self.assertEqual(samples["bkg_dy_ee"].response_mode, "dielectron")

    def test_config_uses_ssc_response_executable_and_four_sample_default(self) -> None:
        with patch.dict("os.environ", {}, clear=True):
            cfg = campaign.parse_config(["--dry-run", "--no-herwig-module"])

        self.assertEqual(cfg.detector_response, "ssc")
        self.assertEqual(cfg.analysis_target, "HwSimPostAnalysis_gammagamma_SSC")
        self.assertEqual(cfg.analysis_exe.name, "HwSimPostAnalysis_gammagamma_SSC")
        self.assertEqual(
            cfg.run_samples,
            "signal_gg_h_aa,bkg_prompt_aa,bkg_gamma_j,bkg_dy_ee",
        )

    def test_none_response_uses_legacy_executable_and_environment_default(self) -> None:
        with patch.dict("os.environ", {"DETECTOR_RESPONSE": "none"}, clear=True):
            cfg = campaign.parse_config(["--dry-run", "--no-herwig-module"])

        self.assertEqual(cfg.detector_response, "none")
        self.assertEqual(cfg.analysis_target, "HwSimPostAnalysis_gammagamma")
        self.assertEqual(cfg.analysis_exe.name, "HwSimPostAnalysis_gammagamma")

    def test_invalid_environment_detector_response_is_rejected(self) -> None:
        with patch.dict("os.environ", {"DETECTOR_RESPONSE": "typo"}, clear=True):
            with patch("sys.stderr", new=io.StringIO()):
                with self.assertRaises(SystemExit):
                    campaign.parse_config(["--dry-run", "--no-herwig-module"])

    def test_ssc_analysis_command_includes_response_mode_and_seed(self) -> None:
        cfg = campaign.parse_config(["--dry-run", "--no-herwig-module"])
        sample = next(sample for sample in campaign.SAMPLES if sample.name == "bkg_gamma_j")

        command = campaign.analysis_command(sample, Path("roots.input"), 1234, cfg)

        self.assertIn("--response-mode", command)
        self.assertIn("gammajet", command)
        self.assertIn("--seed", command)
        self.assertIn(1234, command)

    def test_none_analysis_command_omits_ssc_only_options(self) -> None:
        cfg = campaign.parse_config(
            ["--dry-run", "--no-herwig-module", "--detector-response", "none"]
        )
        sample = next(sample for sample in campaign.SAMPLES if sample.name == "signal_gg_h_aa")

        command = campaign.analysis_command(sample, Path("roots.input"), 1234, cfg)

        self.assertNotIn("--response-mode", command)
        self.assertNotIn("--seed", command)
        self.assertEqual(command[0], cfg.analysis_exe)

    def test_none_response_accepts_genuine_samples_and_rejects_fake_samples(self) -> None:
        cfg = campaign.parse_config(
            [
                "--dry-run",
                "--no-herwig-module",
                "--detector-response",
                "none",
                "--run-samples",
                "signal_gg_h_aa,bkg_prompt_aa",
            ]
        )
        campaign.validate_detector_response(cfg, campaign.selected_samples(cfg))

        fake_cfg = campaign.parse_config(
            [
                "--dry-run",
                "--no-herwig-module",
                "--detector-response",
                "none",
                "--run-samples",
                "bkg_gamma_j,bkg_dy_ee",
            ]
        )
        with self.assertRaisesRegex(SystemExit, "cannot model the photon fakes.*bkg_gamma_j,bkg_dy_ee"):
            campaign.validate_detector_response(fake_cfg, campaign.selected_samples(fake_cfg))

    def test_module_takes_precedence_over_exported_herwig_prefix(self) -> None:
        with patch.dict("os.environ", {"HERWIG_ENV": "/home/shared/Herwig"}, clear=False):
            cfg = campaign.parse_config(["--herwig-module", "herwig/stable"])

        self.assertEqual(cfg.herwig_module, "herwig/stable")
        self.assertIsNone(cfg.herwig_env)

    def test_exported_herwig_directory_uses_linux_module_default(self) -> None:
        with TemporaryDirectory() as tmpdir:
            with patch.dict("os.environ", {"HERWIG_ENV": tmpdir}, clear=True):
                with patch.object(campaign.platform, "system", return_value="Linux"):
                    cfg = campaign.parse_config([])

        self.assertEqual(cfg.herwig_module, "herwig/stable")
        self.assertIsNone(cfg.herwig_env)

    def test_rewrite_run_card_include_updates_pdf_assignments(self) -> None:
        stale_include = """      GRIDPACK = .FALSE.

      PDLABEL = 'nn23lo1'

      PDSUBLABEL(1) = 'nn23lo1'

      PDSUBLABEL(2) = 'nn23lo1'

      LHAID = 230000

      PTJ = 2.000000000000000D+01
"""

        updated = campaign.rewrite_run_card_include(stale_include)

        self.assertIn("      PDLABEL = 'lhapdf'", updated)
        self.assertIn("      PDSUBLABEL(1) = 'lhapdf'", updated)
        self.assertIn("      PDSUBLABEL(2) = 'lhapdf'", updated)
        self.assertIn("      LHAID = 331900", updated)
        self.assertIn("      PTJ = 2.000000000000000D+01", updated)
        self.assertNotIn("nn23lo1", updated)

    def test_patch_run_card_include_removes_stale_setrun_object(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            source_dir = process_dir / "Source"
            source_dir.mkdir()
            include_path = source_dir / "run_card.inc"
            include_path.write_text("      PDLABEL = 'nn23lo1'\n      LHAID = 230000\n")
            setrun_object = source_dir / "setrun.o"
            setrun_object.write_bytes(b"old object")

            campaign.patch_run_card_include(process_dir, SimpleNamespace(dry_run=False))

            updated = include_path.read_text()
            self.assertIn("      PDLABEL = 'lhapdf'", updated)
            self.assertIn("      LHAID = 331900", updated)
            self.assertFalse(setrun_object.exists())

    def test_patch_run_card_include_allows_missing_generated_include(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            source_dir = process_dir / "Source"
            source_dir.mkdir()
            setrun_object = source_dir / "setrun.o"
            setrun_object.write_bytes(b"old object")

            campaign.patch_run_card_include(process_dir, SimpleNamespace(dry_run=False))

            self.assertFalse((source_dir / "run_card.inc").exists())
            self.assertFalse(setrun_object.exists())

    def test_generate_run_card_include_invokes_mg5_make_target(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            source_dir = process_dir / "Source"
            source_dir.mkdir()
            include_path = source_dir / "run_card.inc"

            def fake_run(args, cfg):
                include_path.write_text("      PDLABEL = 'nn23lo1'\n")

            with patch.object(campaign, "run_runtime_command", side_effect=fake_run) as run_mock:
                campaign.generate_run_card_include(process_dir, SimpleNamespace(dry_run=False))

            run_mock.assert_called_once_with(["make", "-C", source_dir, "run_card.inc"], SimpleNamespace(dry_run=False))
            self.assertTrue(include_path.exists())

    def test_patch_mg5_pdf_label_sync_updates_generated_driver_once(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            interface = process_dir / "bin" / "internal" / "madevent_interface.py"
            interface.parent.mkdir(parents=True)
            interface.write_text(
                "        # set  lhapdf.\n"
                "        if self.run_card['pdlabel'] == \"lhapdf\":\n"
                "            self.make_opts_var['lhapdf'] = 'True'\n"
                "        # create param_card.inc and run_card.inc\n"
                "        self.do_treatcards('')\n"
                "        \n"
                "        logger.info(\"compile Source Directory\")\n"
            )

            campaign.patch_mg5_pdf_label_sync(process_dir, SimpleNamespace(dry_run=False))
            campaign.patch_mg5_pdf_label_sync(process_dir, SimpleNamespace(dry_run=False))

            updated = interface.read_text()
            self.assertEqual(updated.count(campaign.MG5_PDF_LABEL_SYNC_MARKER), 1)
            self.assertEqual(updated.count(campaign.MG5_POST_TREATCARDS_PDF_PATCH_MARKER), 1)
            self.assertIn("self.run_card['pdlabel1'] == \"lhapdf\"", updated)
            self.assertIn("self.run_card['pdlabel2'] == \"lhapdf\"", updated)
            self.assertLess(
                updated.index(campaign.MG5_PDF_LABEL_SYNC_MARKER),
                updated.index("        if self.run_card['pdlabel'] == \"lhapdf\":"),
            )
            self.assertLess(
                updated.index("        self.do_treatcards('')"),
                updated.index(campaign.MG5_POST_TREATCARDS_PDF_PATCH_MARKER),
            )
            self.assertIn("run_card_inc = pjoin(self.me_dir, 'Source', 'run_card.inc')", updated)
            self.assertIn("setrun_object = pjoin(self.me_dir, 'Source', 'setrun.o')", updated)

    def test_patch_mg5_pdf_label_sync_upgrades_driver_with_existing_marker(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            interface = process_dir / "bin" / "internal" / "madevent_interface.py"
            interface.parent.mkdir(parents=True)
            interface.write_text(
                f"        {campaign.MG5_PDF_LABEL_SYNC_MARKER}\n"
                "        # create param_card.inc and run_card.inc\n"
                "        self.do_treatcards('')\n"
            )

            campaign.patch_mg5_pdf_label_sync(process_dir, SimpleNamespace(dry_run=False))

            updated = interface.read_text()
            self.assertEqual(updated.count(campaign.MG5_PDF_LABEL_SYNC_MARKER), 1)
            self.assertEqual(updated.count(campaign.MG5_POST_TREATCARDS_PDF_PATCH_MARKER), 1)

    def test_patch_mg5_pdf_defaults_updates_hidden_banner_defaults(self) -> None:
        with TemporaryDirectory() as tmpdir:
            process_dir = Path(tmpdir)
            banner = process_dir / "bin" / "internal" / "banner.py"
            banner.parent.mkdir(parents=True)
            banner.write_text(
                '        self.add_param("pdlabel", "nn23lo1", hidden=True, allowed=valid_pdf)\n'
                '        self.add_param("pdlabel1", "nn23lo1", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(1)")\n'
                '        self.add_param("pdlabel2", "nn23lo1", hidden=True, allowed=valid_pdf, fortran_name="pdsublabel(2)")\n'
                '        self.add_param("lhaid", 230000, hidden=True)\n'
            )

            campaign.patch_mg5_pdf_defaults(process_dir, SimpleNamespace(dry_run=False))
            campaign.patch_mg5_pdf_defaults(process_dir, SimpleNamespace(dry_run=False))

            updated = banner.read_text()
            self.assertIn('self.add_param("pdlabel", "lhapdf"', updated)
            self.assertIn('self.add_param("pdlabel1", "lhapdf"', updated)
            self.assertIn('self.add_param("pdlabel2", "lhapdf"', updated)
            self.assertIn('self.add_param("lhaid", 331900, hidden=True)', updated)
            self.assertNotIn("nn23lo1", updated)
            self.assertNotIn("230000", updated)


if __name__ == "__main__":
    unittest.main()
