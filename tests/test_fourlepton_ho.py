"""Synthetic bookkeeping checks; no generators, showers or ML fitting."""
from dataclasses import replace
import gzip
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

import run_pipeline
from hfourlepton import higher_order as ho, workflow
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import analyze_fourlepton as analysis
from hfourlepton import validate_lhc
from hfourlepton import make_fourlepton_report as report
from hgammagamma import ho_signal_normalization as norm
from test_fourlepton_campaign import lhe_text
from test_ihixs_normalization import synthetic_record
import test_ho_shower_completion as completion_fixtures


def config(*args):
    with patch.dict("os.environ", {}, clear=True):
        return four.parse_config(["--dry-run", "--no-herwig-module", *args])


def cfg(*args):
    return config("--order", "ho", "--signal-modes", "all", *args)


def powheg_signal_lhe(masses=(125.0011342, 124.85, 125.2)):
    events = "".join(
        f"<event>\n1 1 {weight} 100 .007 .118\n"
        f"25 1 1 2 0 0 0 0 0 {mass} {mass} 0 9\n</event>\n"
        for mass, weight in zip(masses, (10., -2., 8.))
    )
    return ('<LesHouchesEvents version="3.0">\n<!--\n'
            'hmass 125\nhwidth 0.004152d0\nbwcutoff 1000\nhiggsfixedwidth 1\n'
            'minnlo 1\nlhans1 336100\nlhans2 336100\n-->\n<header>\n</header>\n'
            '<init>\n2212 2212 6500 6500 -1 -1 -1 -1 -4 1\n'
            '-1 .1 10 1\n</init>\n' + events + '</LesHouchesEvents>\n')


def powheg_fixture(root):
    record = synthetic_record()
    record.update(profile="lhc13-ho-pure-heft-v1", sqrt_s_gev=13000.)
    record["native_parameters"].update(ebeam1=6500., ebeam2=6500.)
    record_path = root / "normalization.json"
    record_path.write_text(json.dumps(norm.seal_record(record)))
    source = root / "signal.lhe"
    source.write_text(powheg_signal_lhe())
    settings = replace(cfg("--nevents", "1", "--ebeam", "6500",
        "--signal-lhe", str(source), "--signal-normalization", str(record_path)), dry_run=False)
    return source, settings


class FourLeptonHOTests(unittest.TestCase):
    def test_powheg_width_before_and_after_independent_decay_slices(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            source, settings = powheg_fixture(Path(tmp))
            original = source.read_bytes(), source.stat().st_mtime_ns
            samples = ho.all_samples()[:3]
            info = ho.inspect_production(source, samples[0], settings)
            self.assertEqual(info["events"], 3)
            self.assertEqual(info["negative_events"], 1)
            for sample in samples:
                sliced = ho.obtain_production(sample, samples, settings)
                mass, width = four.validate_stable_higgs_lhe(sliced, settings, 1, sample=sample)
                self.assertEqual(mass, 125.)
                self.assertEqual(width, .004152)
                self.assertEqual(gzip.open(sliced, "rt").read().count("<event>"), 1)
            self.assertEqual((source.read_bytes(), source.stat().st_mtime_ns), original)

    def test_powheg_width_requires_explicit_compatible_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, settings = powheg_fixture(Path(tmp))
            text = source.read_text()
            mutations = [("hmass 125", "hmass 125.09"),
                         ("hwidth 0.004152d0", "hwidth 0.005"),
                         ("hmass 125\n", ""), ("hwidth 0.004152d0\n", ""),
                         ("bwcutoff 1000\n", ""), ("minnlo 1", "minnlo 0"),
                         ("bwcutoff 1000", "bwcutoff nan"),
                         ("bwcutoff 1000", "bwcutoff inf"),
                         ("bwcutoff 1000", "bwcutoff -1"),
                         ("bwcutoff 1000", "bwcutoff 0"),
                         ("bwcutoff 1000", "bwcutoff 1000\nbwcutoff 1")]
            for old, new in mutations:
                with self.subTest(replacement=new):
                    source.write_text(text.replace(old, new))
                    with self.assertRaises(ValueError):
                        four.validate_stable_higgs_lhe(source, settings, 1, sample=ho.all_samples()[0])

    def test_powheg_event_masses_must_stay_inside_recorded_window(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, settings = powheg_fixture(Path(tmp))
            sample = ho.all_samples()[0]
            for mass in (120.848, 129.152):
                source.write_text(powheg_signal_lhe((mass,)))
                self.assertEqual(four.validate_stable_higgs_lhe(source, settings, 1, sample=sample), (125., .004152))
            for mass in (120.847, 129.153, 0., float("nan"), float("inf")):
                with self.subTest(mass=mass):
                    source.write_text(powheg_signal_lhe((mass,)))
                    with self.assertRaisesRegex(ValueError, "Higgs mass|Higgs.*mass"):
                        four.validate_stable_higgs_lhe(source, settings, 1, sample=sample)
            source.write_text(powheg_signal_lhe().replace("25 1 1 2", "25 2 1 2"))
            with self.assertRaisesRegex(ValueError, "exactly one undecayed Higgs"):
                four.validate_stable_higgs_lhe(source, settings, 1, sample=sample)
            source.write_text(powheg_signal_lhe((125.0011342,)))
            with self.assertRaisesRegex(ValueError, "only 1 complete stable-H"):
                four.validate_stable_higgs_lhe(source, settings, 2, sample=sample)

    def test_finite_width_exception_does_not_relax_other_signal_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            source, settings = powheg_fixture(Path(tmp))
            for sample in (None, ho.all_samples()[3], four.SAMPLES[0]):
                with self.subTest(sample=sample.name if sample else None):
                    with self.assertRaisesRegex(ValueError, "expected 125.0"):
                        four.validate_stable_higgs_lhe(source, settings, 1, sample=sample)

    def test_inventory_orders_and_independent_decay_statistics(self):
        settings = cfg()
        samples = four.selected_samples(settings)
        self.assertEqual(len(samples), 27)
        self.assertEqual(len({s.production_group_id for s in samples}), 13)
        self.assertEqual({s.production_mode for s in samples if s.is_signal}, {"ggf", "vbf", "wh", "zh", "tth"})
        self.assertEqual(sum(s.matching == "MC@NLO" for s in samples), 18)
        self.assertEqual(sum(s.matching == "POWHEG" for s in samples), 3)
        for s in samples:
            count = four.group_event_count(s.production_group_id, samples, settings)
            self.assertEqual(count, settings.nevents * (3 if s.is_signal else 1))
        self.assertIn("$$ w+ w- z", next(s.process for s in samples if s.production_mode == "vbf"))
        self.assertEqual(settings.higgs_mass, 125.)
        self.assertEqual(config().higgs_mass, 125.09)

    def test_energy_profiles_cannot_share_pilot_evidence(self):
        lhc = run_pipeline.parse_args(["prepare", "--workflow", "fourlepton-ho", "--run-tag", "test", "--energy-tev", "13"])
        ssc = run_pipeline.parse_args(["prepare", "--workflow", "fourlepton-ho", "--run-tag", "test", "--energy-tev", "40"])
        self.assertNotEqual(run_pipeline.pilot_path(lhc), run_pipeline.pilot_path(ssc))
        self.assertEqual(workflow.profiles(lhc), ("perfect",))
        self.assertEqual(workflow.profiles(ssc), ("perfect", "ssc"))
        settings, card = workflow.rate_inputs(lhc)
        self.assertEqual(settings["sqrt_s_gev"], 13000.)
        self.assertEqual(settings["profile"], "lhc13-ho-pure-heft-v1")
        self.assertIn("6500", card)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "powheg.input"
            path.write_text(card)
            self.assertEqual(workflow.ihixs.native_parameters(path, 13000.)["ebeam1"], 6500.)
            path = Path(tmp) / "ihixs.json"
            path.write_text(json.dumps(settings))
            self.assertEqual(workflow.ihixs.load_settings(path)["sqrt_s_gev"], 13000.)

    def test_managed_run_wires_three_disjoint_decay_inputs_without_launching(self):
        for energy, profile_count in ((13, 1), (40, 2)):
            with self.subTest(energy=energy), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                record = root / "record.json"
                record.write_text("{}")
                args = SimpleNamespace(energy_tev=energy, run_tag="synthetic", nevents=30, jobs=4,
                    signal_normalization=record, signal_lhe=None)
                calls = []
                pipeline = SimpleNamespace(args=args, check=lambda: None, work=root, python=Path("/synthetic/python"),
                    powheg=root / "powheg", mg5=root / "mg5", call=lambda label, argv: calls.append((label, argv)))
                with patch.object(workflow, "REPO", root), patch.object(workflow.shared.norm, "load_record", return_value={"sqrt_s_gev": energy*1000.}):
                    workflow.run(pipeline)
                self.assertEqual(len(calls), 1 + profile_count)
                powheg_args = calls[0][1]
                self.assertEqual(powheg_args[powheg_args.index("--nevents")+1], 90)
                for _, command in calls[1:]:
                    self.assertEqual(command[command.index("--nevents")+1], 30)
                    self.assertEqual(command[command.index("--signal-modes")+1], "all")

    def test_nlo_card_matching_pdf_and_loose_cuts(self):
        sample = next(s for s in ho.all_samples() if s.name == "bkg_qq_2e2mu")
        proc, launch = ho.nlo_cards(sample, cfg("--ebeam", "6500"), Path("/tmp/process"), 100)
        self.assertIn("/ h [QCD]", proc)
        for token in ("aMC@NLO --parton", "HERWIGPP", "set ebeam1 6500", "set lhaid 335900", "set ptl 2.0", "set event_norm average"):
            self.assertIn(token, launch)
        shower = ho.herwig_card(sample, Path("/tmp/events.lhe"), 42, 100, cfg())
        for token in ("AllowedToReOpen No", "WeightOption VarNegWeight", "SaveObjects Yes", "SaveOptionalWeights Yes"):
            self.assertIn(token, shower)
        self.assertIn("InitialStateReconOption Rapidity", shower)

    def test_vbf_orders_and_preview_cards_use_the_same_mixed_loop_policy(self):
        for energy in (13, 40):
            with self.subTest(energy=energy), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                args = run_pipeline.parse_args(["prepare", "--workflow", "fourlepton-ho",
                    "--run-tag", "vbf_orders", "--energy-tev", str(energy), "--nevents", "30"])
                settings = cfg("--run-tag", args.run_tag, "--ebeam", str(500 * energy),
                               "--nevents", "30", "--nb-core", str(args.jobs))
                with patch.object(workflow, "REPO", root), patch.object(workflow, "configuration", return_value={}):
                    plan_dir = workflow.prepare(args)
                inventory = json.loads((plan_dir / "plan.json").read_text())
                samples = ho.selected_samples(settings)
                for sample in samples:
                    if sample.matching != "MC@NLO":
                        continue
                    process_dir = root / "hfourlepton/runs" / args.run_tag / "production" / sample.production_group_id / "mg5_process"
                    process, launch = ho.nlo_cards(sample, settings, process_dir,
                        four.group_event_count(sample.production_group_id, samples, settings))
                    self.assertEqual((plan_dir / sample.production_group_id / "process.mg5").read_text(), process)
                    self.assertEqual((plan_dir / sample.production_group_id / "launch.mg5").read_text(), launch)
                    if sample.production_mode == "vbf":
                        declaration = next(x for x in process.splitlines() if x.startswith("generate "))
                        self.assertIn("QCD^2=0 QED^2=6 [QCD]", declaration)
                        self.assertNotIn("QCD=0", declaration)
                        self.assertIn("$$ w+ w- z", declaration)
                        self.assertIn("set nlo_mixed_expansion True\n", process)
                        self.assertEqual(next(x['process'] for x in inventory['samples'] if x['name'] == sample.name), sample.process)
                    else:
                        self.assertNotIn("set nlo_mixed_expansion", process)
                    self.assertNotIn("IRPoleCheckThreshold", process + launch)

    def test_missing_nlo_output_reports_pole_failures_and_debug_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name, content in (("P0_bad", " FAILED 1e-5\n" * 20),
                                  ("P0_good", " PASSED 1e-5\n" * 20),
                                  ("P0_within_tolerance", " FAILED 1e-5\n" + " PASSED 1e-5\n" * 19)):
                folder = root / "SubProcesses" / name
                folder.mkdir(parents=True)
                (folder / "check_poles.log").write_text(content)
            debug = root / "pilot_tag_1_debug.log"
            debug.write_text("Poles do not cancel, run cannot continue\n")
            message = ho.nlo_failure_message(root, "pilot", root / "Events/pilot/events.lhe.gz")
            self.assertIn("NLO sample is invalid", message)
            self.assertIn("P0_bad/check_poles.log (20/20 points failed)", message)
            self.assertIn(str(debug), message)
            self.assertNotIn("P0_good", message)
            self.assertNotIn("P0_within_tolerance", message)
            empty = root / "missing"
            self.assertEqual(ho.nlo_failure_message(empty, "pilot", empty / "events.lhe.gz"),
                             f"MG5 did not produce {empty / 'events.lhe.gz'}")

    def test_vbf_polynomial_build_setting_is_scoped_and_idempotent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            virtual = root / "SubProcesses/P0_uu_huu/V0_uu_huu"
            virtual.mkdir(parents=True)
            makefile = virtual / "makefile"
            original = ("FFLAGS = -O\nPOLYNOMIAL_OPTIMIZATION = -O3\n"
                        "POLYNOMIAL_BOUNDS_CHECK =\n"
                        "%olynomial.o : %olynomial.f\n"
                        "\t$(FC) $(filter-out -O%,$(FFLAGS)) $(POLYNOMIAL_OPTIMIZATION) -c $< -o $@\n")
            makefile.write_text(original)
            other_file = root / "Source/make_opts"
            other_file.parent.mkdir()
            other_file.write_text("FFLAGS = -O\n")
            vbf = next(s for s in ho.all_samples() if s.production_mode == "vbf")
            qq = next(s for s in ho.all_samples() if s.name == "bkg_qq_4e")
            ho.configure_nlo_polynomial_build(qq, root)
            self.assertEqual(makefile.read_text(), original)
            ho.configure_nlo_polynomial_build(vbf, root)
            self.assertEqual(makefile.read_text(), original.replace("= -O3", "= -O1"))
            self.assertEqual(other_file.read_text(), "FFLAGS = -O\n")
            changed_time = makefile.stat().st_mtime_ns
            ho.configure_nlo_polynomial_build(vbf, root)
            self.assertEqual(makefile.stat().st_mtime_ns, changed_time)

    def test_vbf_build_rejects_unknown_settings_and_external_makefiles_before_writing(self):
        vbf = next(s for s in ho.all_samples() if s.production_mode == "vbf")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "export"
            good = root / "SubProcesses/P0_a/V0_a/makefile"
            bad = root / "SubProcesses/P0_b/V0_b/makefile"
            good.parent.mkdir(parents=True)
            bad.parent.mkdir(parents=True)
            good.write_text("POLYNOMIAL_OPTIMIZATION = -O3\n")
            bad.write_text("POLYNOMIAL_OPTIMIZATION = -Ofast\n")
            with self.assertRaisesRegex(ValueError, "unrecognized"):
                ho.configure_nlo_polynomial_build(vbf, root)
            self.assertEqual(good.read_text(), "POLYNOMIAL_OPTIMIZATION = -O3\n")
            external = Path(tmp) / "shared_makefile"
            external.write_text("POLYNOMIAL_OPTIMIZATION = -O3\n")
            bad.unlink()
            bad.symlink_to(external)
            with self.assertRaisesRegex(ValueError, "outside the VBF export"):
                ho.configure_nlo_polynomial_build(vbf, root)
            self.assertEqual(external.read_text(), "POLYNOMIAL_OPTIMIZATION = -O3\n")
            self.assertEqual(good.read_text(), "POLYNOMIAL_OPTIMIZATION = -O3\n")
            with self.assertRaisesRegex(ValueError, "no virtual-subprocess makefiles"):
                ho.configure_nlo_polynomial_build(vbf, Path(tmp) / "empty")

    def test_real_nlo_inspector_contract_with_synthetic_events(self):
        sample = next(s for s in ho.all_samples() if s.name == "bkg_qq_4e")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "events.lhe"
            text = lhe_text([10.], [2., -1.]).replace("<header>", "<header>\nHERWIGPP = parton_shower")
            text = text.replace("331900", "335900").replace(" 2 1 ", " 0 1 ")
            path.write_text(text)
            info = ho.inspect_production(path, sample, cfg("--nevents", "2"))
            self.assertEqual(info["events"], 2)
            self.assertEqual(info["negative_events"], 1)
            self.assertEqual(info["cross_section_pb"], 10.)

    def test_disjoint_slices_keep_event_weights_and_rate(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            settings = replace(cfg("--nevents", "2"), dry_run=False)
            samples = ho.all_samples()[:3]
            source = Path(tmp) / "source.lhe"
            source.write_text(lhe_text([10.], [1., -1., 3., 4., -5., 6.]))
            outputs = [ho.independent_decay_input(source, s, samples, settings) for s in samples]
            payloads = [gzip.open(p, "rt").read() for p in outputs]
            self.assertEqual([text.count("<event>") for text in payloads], [2, 2, 2])
            self.assertTrue(all(four.parse_lhe_init_cross_section_pb(p) == 10. for p in outputs))
            self.assertIn("-1.000000000000e+00", payloads[0])
            self.assertNotIn("-1.000000000000e+00", payloads[1])
            self.assertIn("-5.000000000000e+00", payloads[2])
            self.assertEqual(ho.independent_decay_input(source, samples[0], samples, settings), outputs[0])
            source.write_text(source.read_text().replace("6.000000000000e+00", "7.000000000000e+00"))
            with self.assertRaisesRegex(ValueError, "changed"):
                ho.independent_decay_input(source, samples[0], samples, settings)

    def test_short_stable_pool_cannot_recycle_events(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            settings = replace(cfg("--nevents", "2"), dry_run=False)
            samples = ho.all_samples()[:3]
            source = Path(tmp) / "short.lhe"
            source.write_text(lhe_text([10.], [1., 1.]))
            with self.assertRaisesRegex(ValueError, "at least 4 hard events"):
                ho.independent_decay_input(source, samples[1], samples, settings)

    def test_both_energy_records_are_validated_and_not_interchangeable(self):
        record = synthetic_record()
        record.update(profile="lhc13-ho-pure-heft-v1", sqrt_s_gev=13000.)
        record["native_parameters"].update(ebeam1=6500., ebeam2=6500.)
        record = norm.seal_record(record)
        self.assertIs(norm.validate_record(record), record)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "record.json"
            path.write_text(json.dumps(record))
            self.assertEqual(ho.load_signal_record(cfg("--ebeam", "6500", "--signal-normalization", str(path))), record)
            with self.assertRaisesRegex(ValueError, "energy/mass"):
                ho.load_signal_record(cfg("--signal-normalization", str(path)))
        record["profile"] = norm.PROFILE
        with self.assertRaises(ValueError):
            norm.validate_record(norm.seal_record(record))

    def test_signal_branching_fraction_is_applied_once(self):
        settings = cfg()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "decayed.lhe"
            path.write_text(lhe_text([.01]))
            for sample in (ho.all_samples()[0], ho.all_samples()[3]):
                with patch.object(ho, "inspect_production", return_value={"cross_section_pb": 10.}), patch.object(ho, "load_signal_record", return_value=synthetic_record()):
                    result = ho.normalization(sample, path, settings, path)
                production = 100. if sample.production_mode == "ggf" else 10.
                expected = production * four.signal_target_branching_fraction(sample, settings)
                self.assertAlmostEqual(.01 * result["weight_scale"], expected)
                self.assertEqual("ihixs" in result, sample.production_mode == "ggf")


class FourLeptonAnalysisTests(unittest.TestCase):
    def test_matched_report_keeps_failed_source_weights_in_denominator(self):
        with tempfile.TemporaryDirectory() as tmp:
            fixture_helper = completion_fixtures.ShowerCompletionTests()
            fixture = fixture_helper.fixture(tmp, saved=2, requested=3, attempts=3, termination="source_exhausted")
            completion, manifest = fixture_helper.build(fixture)
            sample = SimpleNamespace(name=manifest["sample"], raw={"order": "ho", "normalization": {},
                "shower_completion": completion, "nevents_requested": 3, "input_lhe": manifest["lhe_file"]},
                branching_ratio_reweighter=False, denominator_sumw=1.6, denominator_sumw2=1.68,
                generated_sumw=None, generated_sumw2=None, native_cross_section_pb=10.,
                target_cross_section_pb=10., weight_scale=1.)
            generated = report.WeightStats()
            for weight in (1., -.2):
                generated.add(weight)
            audit = report.Audit()
            result = report._normalization(sample, generated, generated, audit)
            self.assertFalse(audit.errors, audit.errors)
            self.assertAlmostEqual(result["selected_cross_section_pb"], 5.)
            sample.denominator_sumw = .8
            audit = report.Audit()
            report._normalization(sample, generated, generated, audit)
            self.assertTrue(any("consumed signed source" in e for e in audit.errors))

    def test_cards_cover_all_modes_and_exclude_mass_training(self):
        for energy in ("lhc13", "ssc40"):
            for mode in ("cuts", "xgboost"):
                card = analysis.read_card(analysis.REPO / f"hfourlepton/analysis_cards/{energy}_{mode}.json")
                self.assertEqual(set(card["required_samples"]), {s.name for s in ho.all_samples()})
                self.assertNotIn("m4l", card["features"])
                self.assertEqual(card["cuts"]["m4l"], {"min": 120., "max": 130.})

    def test_cut_boundaries_and_signed_source_variances(self):
        data = {"m4l": np.array([120., 125., 130.]), "label": np.array([1, 1, 0]),
                "group": np.array(["signal"]*2 + ["background"]), "source": np.array([0, 0, 0])}
        np.testing.assert_array_equal(analysis.cut_mask(data, {"m4l": {"min": 120., "max": 130.}}), [False, True, False])
        result = analysis.yields(data, np.array([2., -.5, 1.]), np.ones(3, dtype=bool))
        self.assertEqual(result["signal"]["sum_weight"], 1.5)
        self.assertEqual(result["signal"]["variance"], 2.25)

    def test_overlapping_slice_provenance_fails(self):
        samples = []
        for i, (start, stop) in enumerate(((0, 10), (9, 19))):
            name = f"signal{i}"
            samples.append(SimpleNamespace(name=name, category="signal", production_group_id="common", raw={
                "analysis_source_group": name, "source_lhe": name, "source_lhe_sha256": "a"*64, "nevents_requested": 10,
                "hard_event_slice": {"identity": {"kind": "disjoint_hard_event_slice", "sample": name,
                    "start": start, "stop": stop, "source_sha256": "b"*64}, "artifact": {"path": name, "sha256": "a"*64}}}))
        with self.assertRaisesRegex(ValueError, "overlap"):
            analysis.audit_disjoint_slices(samples)
        samples[1].raw["hard_event_slice"]["identity"].update(start=10, stop=20)
        analysis.audit_disjoint_slices(samples)

    def test_small_pilot_cannot_pass_lhc_validation(self):
        meta = [{"name": s.name, "production_mode": s.production_mode, "target_cross_section_pb": .001}
                for s in ho.all_samples() if s.is_signal]
        names = [s["name"] for s in meta]
        data = {"sample": np.array(names), "group": np.array(names), "source": np.zeros(len(names), dtype=int),
                "eligible": np.ones(len(names), dtype=bool), "weight": np.full(len(names), .4)}
        self.assertTrue(all(r["status"] == "insufficient_statistics" for r in validate_lhc.compare(data, meta)))


if __name__ == "__main__":
    unittest.main()
