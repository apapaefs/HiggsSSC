import contextlib
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from hgammagamma import run_gammagamma_ho_campaign as ho
from hgammagamma import make_gammagamma_report as report


def lhe_text(signal=False, weights=(10., -2., 8.), xsec=16. / 3.):
    pdf = 336100 if signal else 335900
    provenance = "minnlo 1\nlhans1 336100\nlhans2 336100" if signal else "HERWIGPP = parton_shower"
    particle = "25 1 1 2 0 0 0 0 0 125 125 0 9" if signal else "11 1 1 2 0 0 0 0 0 25 0 0 9"
    events = "".join(f"<event>\n1 1 {w} 100 .007 .118\n{particle}\n</event>\n" for w in weights)
    return (f'<LesHouchesEvents version="3.0">\n<header>\n{provenance}\n</header>\n'
            f"<init>\n2212 2212 20000 20000 -1 -1 {pdf} {pdf} -4 1\n{xsec} .2 10 1\n</init>\n"
            f"{events}</LesHouchesEvents>\n")


class HigherOrderCampaignTests(unittest.TestCase):
    def args(self, directory, *extra):
        return ho.parse_args(["--output-dir", str(directory), "--no-herwig-module", "--nevents", "3", *extra])

    def test_backgrounds_match_detector_sample_table(self):
        args = self.args("/tmp/test-ho", "--run-samples", "backgrounds")
        self.assertEqual([s.name for s in ho.selected_samples(args)],
                         ["bkg_prompt_aa", "bkg_gamma_j", "bkg_dy_ee"])
        with self.assertRaisesRegex(ValueError, "unknown samples"):
            ho.selected_samples(self.args("/tmp/test-ho", "--run-samples", "typo"))

    def test_nlo_cards_use_herwig_subtraction_and_infrared_safe_photon_cuts(self):
        args = self.args("/tmp/test-ho")
        for sample in ho.SAMPLES[1:]:
            process, launch = ho.mg5_cards(args, sample)
            self.assertIn("loop_sm-no_b_mass", process)
            self.assertIn("b~", process)
            self.assertIn("[QCD]", process)
            self.assertNotIn(" -f", process)  # Never force-overwrite an export.
            self.assertIn("aMC@NLO --parton", launch)
            self.assertIn("set parton_shower HERWIGPP", launch)
            self.assertIn("set lhaid 335900", launch)
            self.assertIn("set gamma_is_j False", launch)
            self.assertIn("set ptgmin 10.0", launch)
            self.assertIn("set R0gamma 0.4", launch)
            self.assertIn("set maxjetflavor 5", launch)
            self.assertIn("set ickkw 0", launch)
            # With an explicit aMC@NLO mode, only the card-edit menu remains.
            self.assertEqual(launch.count("\ndone\n"), 1)
            self.assertLess(launch.index("set nevents"), launch.index("\ndone\n"))
        self.assertEqual(ho.run_settings(args, ho.SAMPLES[1])["ptj"], 0.)
        self.assertEqual(ho.run_settings(args, ho.SAMPLES[2])["ptj"], 10.)
        self.assertEqual(ho.run_settings(args, ho.SAMPLES[3])["ptj"], 0.)

    def test_shower_modes_are_distinct_and_br_is_external_once(self):
        args = self.args("/tmp/test-ho")
        signal = ho.herwig_card(args, ho.SAMPLES[0])
        background = ho.herwig_card(args, ho.SAMPLES[1])
        self.assertIn("SelectDecayModes h0->gamma,gamma;", signal)
        self.assertIn("BranchingRatio 1.0", signal)
        self.assertNotIn("LongTransBoost", signal)
        self.assertIn("LongTransBoost", background)
        self.assertIn("SpinCorrelations No", background)
        self.assertNotIn("SelectDecayModes", background)
        for card in (signal, background):
            self.assertIn("WeightOption VarNegWeight", card)
            self.assertIn("AllowedToReOpen No", card)
            self.assertIn("SavePartons Yes", card)
            self.assertNotIn("MPIExtractor:", card)

    def test_signed_lhe_weights_are_not_replaced_by_absolute_weights(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.lhe.gz"
            with gzip.open(path, "wt") as stream:
                stream.write(lhe_text())
            data = ho.inspect_lhe(path, ho.SAMPLES[3], self.args(tmp))
            self.assertEqual(data["negative_events"], 1)
            self.assertEqual(data["sum_weight"], 16.)
            self.assertEqual(data["sum_abs_weight"], 20.)
            self.assertEqual(data["sum_weight_squared"], 168.)
            self.assertAlmostEqual(data["cross_section_pb"], 16. / 3.)

    def test_powheg_unspecified_init_uses_signed_mean_and_header_pdf(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "powheg.lhe"
            path.write_text(lhe_text(True, xsec=-1).replace("336100 336100 -4", "-1 -1 -4"))
            data = ho.inspect_lhe(path, ho.SAMPLES[0], self.args(tmp))
            self.assertAlmostEqual(data["cross_section_pb"], 16. / 3.)
            self.assertIn("mean signed", data["normalization_source"])

    def test_lhe_rejects_wrong_shower_beams_pdf_decay_and_incomplete_events(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.lhe"
            cases = [
                (lhe_text().replace("HERWIGPP", "PYTHIA8"), ho.SAMPLES[3], "HERWIGPP"),
                (lhe_text().replace("20000", "6500"), ho.SAMPLES[3], "beam energy"),
                (lhe_text().replace("335900", "331900"), ho.SAMPLES[3], "PDF IDs"),
                (lhe_text(True).replace("25 1 1 2", "25 2 1 2"), ho.SAMPLES[0], "undecayed"),
                (lhe_text().replace("</LesHouchesEvents>", ""), ho.SAMPLES[3], "incomplete"),
                (lhe_text().replace("1 1 10.0", "2 1 10.0"), ho.SAMPLES[3], "incomplete"),
            ]
            for text, sample, message in cases:
                with self.subTest(message=message):
                    path.write_text(text)
                    with self.assertRaisesRegex(ValueError, message):
                        ho.inspect_lhe(path, sample, self.args(tmp))

    def test_prepare_is_idempotent_but_rejects_changed_physics(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp)
            first = ho.prepare(args, ho.SAMPLES[0])
            second = ho.prepare(args, ho.SAMPLES[0])
            self.assertEqual(first, second)
            self.assertEqual(first["weight_scale"], 0.00227)
            args.higgs_br = 0.003
            with self.assertRaisesRegex(ValueError, "configuration changed"):
                ho.prepare(args, ho.SAMPLES[0])
            args.higgs_br = 0.00227
            args.ebeam = 6500.
            with self.assertRaisesRegex(ValueError, "configuration changed"):
                ho.prepare(args, ho.SAMPLES[0])

    def test_dry_run_creates_no_products(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            out = Path(tmp) / "absent"
            ho.main(["--output-dir", str(out), "--dry-run"])
            self.assertFalse(out.exists())

    def test_report_reads_ho_normalization_without_mg5_banner(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp) / "signal_gg_h_aa"
            directory.mkdir()
            path = directory / "normalization-test.json"
            path.write_text(json.dumps({"run_tag": "test", "sample": directory.name,
                                        "cross_section_pb": 200., "cross_section_error_pb": 2.}))
            self.assertEqual(report.parse_cross_section(directory, "test"), (200., 2.))
            path.write_text(path.read_text().replace('"test"', '"wrong"'))
            with self.assertRaisesRegex(ValueError, "provenance mismatch"):
                report.parse_cross_section(directory, "test")

    def test_report_preserves_negative_and_cancelling_histogram_bins(self):
        sample = SimpleNamespace(name="nlo", sum_weight=10., cross_section_pb=30., weight_scale=1.)
        for values in ([-2., 1.], [-1., 1.]):
            histogram = SimpleNamespace(y=values)
            self.assertEqual(report.scaled_histogram(histogram, sample, "event_xsec", False),
                             [3. * x for x in values])
            with self.assertRaisesRegex(ValueError, "cannot unit-normalize"):
                report.scaled_histogram(histogram, sample, "event_xsec", True)

    def test_resume_detects_missing_or_changed_lhe(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp)
            sample = ho.SAMPLES[3]
            path = ho.lhe_path(args, sample)
            path.parent.mkdir(parents=True)
            with gzip.open(path, "wt") as stream:
                stream.write(lhe_text())
            manifest = {"lhe": ho.inspect_lhe(path, sample, args)}
            ho.verify_completed(args, sample, "generate", manifest)
            with gzip.open(path, "wt") as stream:
                stream.write(lhe_text(weights=(20., -2., 8.)))
            with self.assertRaisesRegex(ValueError, "LHE changed"):
                ho.verify_completed(args, sample, "generate", manifest)

    def test_fortran_compatibility_patch_is_local_and_idempotent(self):
        with TemporaryDirectory() as tmp:
            process = Path(tmp)
            path = process / "SubProcesses/genps_fks.f"
            path.parent.mkdir()
            original = ("      external derivative,ran2,virtgranny,xinv_redvirtgranny\n"
                        "      der=derivative(virtgranny_red,x)\n")
            path.write_text(original)
            ho.patch_mg5_fortran(process)
            first = path.read_text()
            ho.patch_mg5_fortran(process)
            self.assertEqual(path.read_text(), first)
            self.assertEqual(first.count("      external virtgranny_red\n"), 1)

    def test_runtime_paths_replace_stale_exported_configuration(self):
        with TemporaryDirectory() as tmp:
            process = Path(tmp) / "process"
            config = process / "Cards/amcatnlo_configuration.txt"
            config.parent.mkdir(parents=True)
            config.write_text("# keep this comment\nfastjet = /missing/fastjet-config # old\n"
                              "lhapdf = lhapdf-config\nnb_core = 1\n")
            paths = ["/runtime/bin/fastjet-config", "/runtime/bin/lhapdf-config"]
            ho.configure_mg5_runtime(process, paths)
            ho.configure_mg5_runtime(process, paths)
            text = config.read_text()
            self.assertEqual(text.count("fastjet ="), 1)
            self.assertIn("fastjet = /runtime/bin/fastjet-config", text)
            self.assertIn("lhapdf = /runtime/bin/lhapdf-config", text)
            self.assertIn("nb_core = 1", text)
            self.assertIn("# keep this comment", text)
            config.unlink()
            shared = Path(tmp) / "shared-config.txt"
            shared.write_text("fastjet = shared\n")
            config.symlink_to(shared)
            with self.assertRaisesRegex(ValueError, "shared MG5 configuration"):
                ho.configure_mg5_runtime(process, paths)
            self.assertEqual(shared.read_text(), "fastjet = shared\n")


if __name__ == "__main__":
    unittest.main()
