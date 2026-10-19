"""Report metadata, paired plots, and their physical normalization."""
import contextlib
import csv
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import zipfile
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


class ReportViewTests(unittest.TestCase):
    def campaign(self, root, background=(3., -1., 8.)):
        for name, category, values, cross_section, weight_scale in (
            ("bkg_prompt_aa", "Backgrounds", background, 10., 1.),
            ("signal_gg_h_aa", "Signal", (1., 3., 0.), 2., .25),
        ):
            directory = root / category / "events" / name
            directory.mkdir(parents=True)
            stem = directory / f"{name}-test"
            stem.with_suffix(".top").write_text(
                'NEW FRAME\nTITLE TOP "diphoton invariant mass"\nSET LIMITS X 100 160\n'
                + "".join(f"{x} {value}\n" for x, value in zip((110, 130, 150), values))
                + "HIST\n"
            )
            stem.with_suffix(".dat").write_text(
                f"events_read 10\nevents_with_two_selected_photons 4\nsum_weight 10\n"
                f"sum_diphoton_weight {sum(values)}\nweight_scale {weight_scale}\n"
            )
            (directory / "normalization-test.json").write_text(json.dumps({
                "run_tag": "test", "sample": name, "cross_section_pb": cross_section,
            }))

    def render(self, campaign, output, *options):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(report.main([
                "--analysis-root", str(campaign), "--run-tag", "test",
                "--output-dir", str(output), *options,
            ]), 0)
        return (output / "index.html").read_text()

    def rows(self, path, sample=None):
        with path.open() as stream:
            return [row for row in csv.DictReader(stream)
                    if sample is None or row["sample"] == sample]

    def test_default_pairs_shapes_and_rates_with_independent_downloads(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.campaign(root / "campaign")
            output = root / "report"
            page = self.render(root / "campaign", output,
                               "--signal-scale", "5", "--luminosity-fb", "2")
            self.assertIn('<div class="plot-views paired">', page)
            self.assertLess(page.index("Unit-area shapes"), page.index("Physical cross sections"))
            self.assertLess(page.index('src="plots/m_gg_density.png"'),
                            page.index('src="plots/m_gg.png"'))
            for name in ("m_gg", "m_gg_density"):
                for folder, suffix in (("plots", "png"), ("plots", "svg"), ("data", "csv")):
                    asset = f"{folder}/{name}.{suffix}"
                    self.assertIn(f'href="{asset}"', page)
                    self.assertGreater((output / asset).stat().st_size, 0)
            with zipfile.ZipFile(output / "gammagamma_report_assets.zip") as archive:
                for name in ("m_gg", "m_gg_density"):
                    for suffix in ("png", "svg"):
                        self.assertIn(f"plots/{name}.{suffix}", archive.namelist())

            for sample in ("bkg_prompt_aa", "signal_gg_h_aa"):
                shapes = self.rows(output / "data/m_gg_density.csv", sample)
                self.assertAlmostEqual(sum(float(row["scaled_bin_value"]) * float(row["bin_width"])
                                           for row in shapes), 1.)
                self.assertTrue(all(float(row["display_scale"]) == 1. for row in shapes))
            signal = self.rows(output / "data/m_gg.csv", "signal_gg_h_aa")
            self.assertAlmostEqual(sum(float(row["scaled_bin_value"]) for row in signal), 1.)
            self.assertTrue(all(float(row["display_scale"]) == 5. for row in signal))
            background = self.rows(output / "data/m_gg.csv", "bkg_prompt_aa")
            self.assertEqual([float(row["scaled_bin_value"]) for row in background], [3., -1., 8.])
            summary = self.rows(output / "data/sample_summary.csv", "signal_gg_h_aa")[0]
            self.assertAlmostEqual(float(summary["selected_cross_section_pb"]), .2)
            self.assertAlmostEqual(float(summary["expected_events"]), 400.)

    def test_single_view_options_keep_original_filenames_and_scaling(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.campaign(root / "campaign")
            for name, options in (("density", ["--plot-mode", "density"]),
                                  ("stacked", ["--plot-mode", "stacked"]),
                                  ("legacy", ["--no-density"])):
                with self.subTest(mode=name):
                    output = root / name
                    page = self.render(root / "campaign", output, *options)
                    self.assertNotIn('<div class="plot-views paired">', page)
                    self.assertEqual(page.count('<img src="'), 1)
                    self.assertFalse((output / "data/m_gg_density.csv").exists())
                    signal = self.rows(output / "data/m_gg.csv", "signal_gg_h_aa")
                    total = sum(float(row["scaled_bin_value"]) for row in signal)
                    self.assertAlmostEqual(total, 1. / 20. if name == "density" else .2)

    def test_nonpositive_signed_shape_does_not_block_physical_report(self):
        for values in ((3., -3., 0.), (-3., 1., 0.)):
            with self.subTest(weights=values), TemporaryDirectory() as tmp:
                root = Path(tmp)
                self.campaign(root / "campaign", background=values)
                output = root / "report"
                page = self.render(root / "campaign", output)
                self.assertIn("Unit-area shapes unavailable", page)
                self.assertIn("Non-positive sums for: Prompt gamma gamma", page)
                self.assertIn('src="plots/m_gg.png"', page)
                self.assertNotIn('src="plots/m_gg_density.png"', page)
                self.assertFalse((output / "data/m_gg_density.csv").exists())
                rows = self.rows(output / "data/m_gg.csv", "bkg_prompt_aa")
                self.assertEqual([float(row["scaled_bin_value"]) for row in rows], list(values))


if __name__ == "__main__":
    unittest.main()
