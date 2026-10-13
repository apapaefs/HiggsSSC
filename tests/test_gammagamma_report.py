"""LO and HO report metadata share manifests without sharing accuracy labels."""
import unittest
from hgammagamma import make_gammagamma_report as report


class CampaignContextTests(unittest.TestCase):
    def test_lo_manifest_keeps_lo_label_and_recorded_energy(self):
        self.assertEqual(report.campaign_context([
            {"matching": "LO", "hard_accuracy": "LO", "ebeam_gev": 20000.}
        ]), (False, 40.))

    def test_early_lo_manifest_reads_frozen_configuration(self):
        self.assertEqual(report.campaign_context([
            {"configuration": {"ebeam": 14000.}}
        ]), (False, 28.))

    def test_ho_signal_and_background_identify_higher_order(self):
        for metadata in ({"matching": "POWHEG"}, {"matching": "MC@NLO"},
                         {"hard_accuracy": "NLO QCD + PS"}):
            with self.subTest(metadata=metadata):
                self.assertEqual(report.campaign_context([
                    {**metadata, "ebeam_gev": 20000.}
                ]), (True, 40.))

    def test_legacy_unmanifested_outputs_keep_historical_defaults(self):
        self.assertEqual(report.campaign_context([]), (False, 40.))

    def test_incompatible_energies_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "different beam energies"):
            report.campaign_context([{"ebeam_gev": 20000.}, {"configuration": {"ebeam": 6500.}}])
        for invalid in (float("nan"), float("inf"), -1.):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "invalid beam"):
                report.campaign_context([{"ebeam_gev": invalid}])
        with self.assertRaisesRegex(ValueError, "does not record"):
            report.campaign_context([{}])


if __name__ == "__main__":
    unittest.main()
