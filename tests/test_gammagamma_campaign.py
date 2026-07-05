import unittest
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
from unittest.mock import patch

from hgammagamma import run_gammagamma_campaign as campaign


class GammaGammaCampaignRunCardTests(unittest.TestCase):
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

            with patch.object(campaign, "run", side_effect=fake_run) as run_mock:
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
            )

            campaign.patch_mg5_pdf_label_sync(process_dir, SimpleNamespace(dry_run=False))
            campaign.patch_mg5_pdf_label_sync(process_dir, SimpleNamespace(dry_run=False))

            updated = interface.read_text()
            self.assertEqual(updated.count(campaign.MG5_PDF_LABEL_SYNC_MARKER), 1)
            self.assertIn("self.run_card['pdlabel1'] == \"lhapdf\"", updated)
            self.assertIn("self.run_card['pdlabel2'] == \"lhapdf\"", updated)
            self.assertLess(
                updated.index(campaign.MG5_PDF_LABEL_SYNC_MARKER),
                updated.index("        if self.run_card['pdlabel'] == \"lhapdf\":"),
            )


if __name__ == "__main__":
    unittest.main()
