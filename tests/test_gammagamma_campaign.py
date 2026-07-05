import unittest
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory

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


if __name__ == "__main__":
    unittest.main()
