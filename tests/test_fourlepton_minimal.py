"""Minimal-preset bookkeeping and executable synthetic checks; no simulation."""
import contextlib
from dataclasses import asdict, replace
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_pipeline as cli
from hfourlepton import presets, minimal_runtime as runtime, workflow, higher_order as ho
from hfourlepton import run_four_lepton_campaign as four
from hgammagamma import run_powheg_hjminnlo as powheg
import test_ho_shower_completion as completion_fixtures


def config(*extra):
    with patch.dict("os.environ", {}, clear=True):
        return four.parse_config(["--order", "ho", "--preset", "minimal", "--run-tag", "minimal_test",
                                  "--no-herwig-module", "--no-herwig-env", "--dry-run", *extra])


def args(*extra):
    return cli.parse_args(["prepare", "--workflow", "fourlepton-ho", "--run-tag", "minimal_test",
                           "--preset", "minimal", "--energy-tev", "13", "--jobs", "40", *extra])


def export_build_files(root, sample):
    """Layouts checked against MG5 3.5.15 and Timur's failed WH export."""
    if sample.matching == "MC@NLO":
        files = {"Source/make_opts": "ifdef madloop\nlibcuttools=-lcts -lninja -lcollier\nelse\nlibcuttools=\nendif\n",
                 "SubProcesses/makefile_loop": "include ../../../Source/make_opts\n"
                     "LINKLIBS = -L$(LIBDIR) -lcts -lninja -lcollier\n"
                     "%.o : %.f\n\t$(FC) $(FFLAGS) -c $< -I /example/include\n"}
    else:
        files = {"SubProcesses/MadLoop_makefile_definitions":
                 "LINK_LOOP_LIBS = -L$(LIBDIR) -lcts -lninja -lcollier\nLOOP_INCLUDE = -I /example/include\n"}
    for relative, text in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    return [root / name for name in files]


class MinimalTests(unittest.TestCase):
    def test_counts_modes_group_sums_and_fractional_multiplier(self):
        samples = ho.all_samples()
        cfg = config("--nb-core", "40")
        self.assertEqual(set(presets.counts()), {s.name for s in samples})
        self.assertEqual(sum(cfg.events_for(s) for s in samples), 231000)
        expected = {"signal_gg_h_stable": 60000, "signal_vbf_h_stable": 30000,
                    "signal_tt_h_stable": 15000, "bkg_qq_4e": 20000, "bkg_gg_4mu": 2000}
        for group, count in expected.items():
            self.assertEqual(four.group_event_count(group, samples, cfg), count)
            twice = config("--event-multiplier", "2")
            self.assertEqual(four.group_event_count(group, samples, twice), 2 * count)
        self.assertEqual(sum(presets.counts("0.001").values()), 231)
        self.assertEqual(presets.counts("0.00051")["bkg_gg_4e"], 2)
        self.assertEqual(presets.max_events(args("--event-multiplier", ".001")), 20)
        self.assertEqual(presets.max_events(args()), 20000)
        self.assertEqual(len(ho.selected_samples(cfg)), 27)
        self.assertEqual(presets.description(cfg)["cores_per_decay"], 5)

    def test_invalid_or_ambiguous_options_and_frozen_identity(self):
        for options in (("--event-multiplier", "nan"), ("--event-multiplier", "inf"),
                        ("--event-multiplier", "-1"), ("--event-multiplier", "0"),
                        ("--event-multiplier", "100000"), ("--nevents", "30")):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                args(*options)
        for options in (("--force",), ("--sample-nevents", "bkg_gg_4e=10"), ("--nb-core", "0")):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                config(*options)
        self.assertEqual(cli.identity(args(), {}), cli.identity(args("--event-multiplier", "1.00"), {}))
        for options in (("--event-multiplier", "2"), ("--wide-continuum",), ("--jobs", "20")):
            self.assertNotEqual(cli.identity(args(), {}), cli.identity(args(*options), {}))
        with tempfile.TemporaryDirectory() as tmp, patch.object(cli, "REPO", Path(tmp)):
            original = args("--work-dir", str(Path(tmp) / "work"))
            cli.prepare_run(original, {})
            changed = args("--work-dir", str(Path(tmp) / "work"), "--event-multiplier", "2")
            with self.assertRaisesRegex(ValueError, "configuration changed"):
                cli.prepare_run(changed, {})

    def test_prepare_is_only_metadata_and_contains_exact_runtime_cards(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(workflow, "REPO", Path(tmp)), \
             patch.object(workflow, "configuration", return_value={}), \
             patch.object(four, "run_runtime_command", side_effect=AssertionError("generator launched")):
            destination = workflow.prepare(args())
            plan = json.loads((destination / "plan.json").read_text())
            self.assertEqual(sum(s["nevents"] for s in plan["samples"]), 231000)
            self.assertEqual(plan["production_events"]["signal_gg_h_stable"], 60000)
            for sample in ho.all_samples():
                card = (destination / sample.name / "herwig.in").read_text()
                self.assertIn(str(presets.counts()[sample.name]), card)
            nlo = (destination / "bkg_qq_4e/launch.mg5").read_text()
            self.assertIn("set req_acc 0.03\n", nlo)
            self.assertNotIn("req_acc_FO", nlo)
            self.assertIn("set output_dependencies external --no_save", (destination / "bkg_qq_4e/process.mg5").read_text())
            self.assertIn("7|6|1", (destination / "bkg_gg_4e/MadLoopParams.override").read_text())
            self.assertIn("istatus", (destination / "bkg_qq_4e/dummy_cuts.f").read_text())
            with self.assertRaises(ValueError):
                workflow.prepare(args("--event-multiplier", "2"))

    def test_managed_commands_use_scaled_production_counts(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as tmp, patch.object(workflow, "REPO", Path(tmp)), \
             patch.object(workflow, "prepare"), patch.object(runtime, "validate_libraries"), \
             patch.object(workflow.shared.norm, "load_record", return_value={"sqrt_s_gev": 13000.}):
            options = args("--event-multiplier", "2")
            calls = []
            pipe = SimpleNamespace(args=options, check=lambda: None, work=Path(tmp), mg5=Path(tmp) / "mg5",
                    python=sys.executable, prefix=Path(tmp), cuba=Path(tmp), boost=Path(tmp),
                    powheg=Path(tmp) / "powheg", config={"cc": "cc", "cxx": "c++"},
                    call=lambda label, command: calls.append((label, command)))
            workflow.run(pipe)
            commands = dict(calls)
            for suffix, count in (("-powheg", 120000), ("-powheg-vbf", 60000)):
                command = commands[options.run_tag + suffix]
                self.assertEqual(command[command.index("--nevents") + 1], count)
            command = commands[options.run_tag + "-fourlepton-perfect"]
            parsed = four.parse_config([str(x) for x in command[2:]])
            self.assertEqual(parsed.sample_nevents, presets.counts(2))
            self.assertEqual(parsed.seed_base, presets.seed_base(options.run_tag))

    def test_managed_external_lhe_uses_normalization_and_checks_all_channels(self):
        from types import SimpleNamespace
        record_data = {"sqrt_s_gev": 13000., "higgs_mass_gev": 125., "native_parameters": {}}
        for explicit_record in (False, True):
            for short_mode in (None, "ggf", "vbf"):
                with self.subTest(explicit_record=explicit_record, short_mode=short_mode), \
                     tempfile.TemporaryDirectory() as tmp, patch.object(workflow, "REPO", Path(tmp)), \
                     patch.object(workflow, "prepare"), patch.object(runtime, "validate_libraries"), \
                     patch.object(workflow.shared.norm, "load_record", return_value=record_data) as load_record:
                    root = Path(tmp).resolve()
                    record = root / ("custom.json" if explicit_record else "normalization/ggf-lhc13-n3lo.json")
                    record.parent.mkdir(parents=True, exist_ok=True)
                    record.write_text("{}")
                    options = args("--signal-lhe", str(root / "ggf.lhe"), "--vbf-lhe", str(root / "vbf.lhe"),
                                   *(["--signal-normalization", str(record)] if explicit_record else []))
                    calls, inspected = [], []
                    pipe = SimpleNamespace(args=options, check=lambda: None, work=root, mg5=root / "mg5",
                            python=sys.executable, call=lambda label, command: calls.append((label, command)))

                    def inspect(source, sample, cfg):
                        self.assertFalse(cfg.dry_run)
                        self.assertEqual(cfg.signal_normalization, record)
                        self.assertEqual(ho.load_signal_record(cfg), record_data)
                        mode = sample.production_mode
                        self.assertEqual(source, root / f"{mode}.lhe")
                        inspected.append(mode)
                        count = {"ggf": 60000, "vbf": 30000}[mode]
                        return {"events": count - (mode == short_mode)}

                    with patch.object(ho, "inspect_production", side_effect=inspect):
                        if short_mode:
                            with self.assertRaisesRegex(ValueError, "across its decay channels"):
                                workflow.run(pipe)
                            self.assertEqual(calls, [])
                        else:
                            workflow.run(pipe)
                            self.assertEqual(inspected, ["ggf", "vbf"])
                            self.assertEqual([label for label, _ in calls], [options.run_tag + "-fourlepton-perfect"])
                            command = calls[0][1]
                            self.assertEqual(command[command.index("--signal-normalization") + 1], record)
                    self.assertTrue(load_record.called)
                    self.assertTrue(all(call.args == (record,) for call in load_record.call_args_list))

    def test_loop_card_and_cut_patch_preserve_stability_and_other_hooks(self):
        for name in ("signal_wp_h_4e", "signal_wm_h_4e", "signal_z_h_4e", "signal_tt_h_4e",
                     "bkg_qq_4mu", "bkg_gg_4e", "signal_gg_z_h_4e"):
            sample = next(s for s in ho.all_samples() if s.name == name)
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / "Cards").mkdir()
                (root / "SubProcesses").mkdir()
                export_build_files(root, sample)
                if sample.matching == "MC@NLO":
                    self.assertFalse((root / "SubProcesses/MadLoop_makefile_definitions").exists())
                ml = root / "Cards/MadLoopParams.dat"
                ml.write_text("#MLReductionLib\n! documentation\n6|1\n#CTModeRun\n-1\n#MLStabThres\n1.0d-3\n")
                hook = root / "SubProcesses/dummy_fct.f"
                original = "      logical function dummy_cuts(p)\n      dummy_cuts=.true.\n      end\n\n      subroutine untouched\n      end\n"
                hook.write_text(original)
                runtime.configure_export(root, sample, replace(config(), dry_run=False))
                self.assertIn("#CTModeRun\n-1\n#MLStabThres\n1.0d-3", ml.read_text())
                self.assertIn("7|6|1" if name.startswith("bkg_gg") else "6|7|1", ml.read_text())
                self.assertIn("      subroutine untouched\n      end", hook.read_text())
                self.assertEqual(hook.read_text() == original, sample.is_signal)

    def test_export_audit_rejects_missing_files_libraries_and_include_paths(self):
        for name in ("signal_wp_h_4e", "bkg_gg_4e"):
            sample = next(s for s in ho.all_samples() if s.name == name)
            with tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                files = export_build_files(root, sample)
                runtime.validate_export_libraries(root, sample)
                for path in files:
                    original = path.read_text()
                    for provider in ("collier", "ninja"):
                        path.write_text(original.replace("-l" + provider, "") + f"\n# -l{provider}\n")
                        with self.assertRaisesRegex(ValueError, "did not enable both"):
                            runtime.validate_export_libraries(root, sample)
                    path.write_text(original + "\nLOOP_INCLUDE = -I <Not_found_define_it_yourself>\n")
                    with self.assertRaisesRegex(ValueError, "include directory"):
                        runtime.validate_export_libraries(root, sample)
                    path.unlink()
                    with self.assertRaisesRegex(ValueError, "incomplete.*MG5 export: missing"):
                        runtime.validate_export_libraries(root, sample)
                    path.write_text(original)

    def test_nlo_export_does_not_accept_a_stale_lo_build_file(self):
        sample = next(s for s in ho.all_samples() if s.name == "signal_wp_h_4e")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_build_files(root, next(s for s in ho.all_samples() if s.name == "bkg_gg_4e"))
            with self.assertRaisesRegex(ValueError, "Source/make_opts"):
                runtime.validate_export_libraries(root, sample)

    def test_nlo_export_requires_flags_on_the_actual_link_assignment(self):
        sample = next(s for s in ho.all_samples() if s.name == "signal_wp_h_4e")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            export_build_files(root, sample)
            (root / "Source/make_opts").write_text("libcuttools=\nUNUSED = -lninja -lcollier\n")
            with self.assertRaisesRegex(ValueError, "did not enable both"):
                runtime.validate_export_libraries(root, sample)

    @unittest.skipUnless(shutil.which("gfortran"), "requires gfortran for the synthetic kinematic check")
    def test_compiled_mass_window_ignores_beams_radiation_and_zero_born_slot(self):
        for nlo in (False, True):
            with self.subTest(nlo=nlo), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                n = 7 if nlo else 6
                (root / "nexternal.inc").write_text(f"      integer nexternal\n      parameter(nexternal={n})\n")
                (root / "cut.f").write_text(runtime.continuum_cut(nlo))
                call = "dummy_cuts(p,istatus,ipdg)" if nlo else "dummy_cuts(p)"
                driver = f"""program test_cut
implicit none
double precision p(0:{4 if nlo else 3},{n}),mass
logical dummy_cuts,answer
integer i,k,istatus({n}),ipdg({n})
double precision, parameter :: masses(5)=[80d0,90d0,125d0,180d0,190d0]
do k=1,5
  p=0d0
  p(0,1:2)=6500d0
  istatus=1
  istatus(1:2)=-1
  ipdg=21
  ipdg(3:6)=[11,-11,13,-13]
  mass=masses(k)
  p(0,3:6)=mass/4d0
  p(1,3)=mass/4d0
  p(1,4)=-mass/4d0
  p(2,5)=mass/4d0
  p(2,6)=-mass/4d0
  answer={call}
  if(answer.neqv.(k.ge.2.and.k.le.4)) stop 1
  {"p(0,7)=1000d0" if nlo else "continue"}
  if({call}.neqv.answer) stop 2
  {"p(:,7)=p(:,4); p(:,4)=0d0; ipdg(7)=-11; ipdg(4)=21" if nlo else "continue"}
  if({call}.neqv.answer) stop 3
enddo
end program
"""
                (root / "test.f90").write_text(driver)
                subprocess.run(["gfortran", "-fcheck=all", "cut.f", "test.f90", "-o", "test"], cwd=root, check=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                subprocess.run([root / "test"], check=True, capture_output=True)

    def test_library_preflight_and_shower_profile_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ninja, collier = root / "ninja/lib", root / "collier"
            for path in (ninja / "libninja.a", ninja / "libavh_olo.a", ninja.parent / "include/mninja.mod",
                         collier / "libcollier.a", collier / "modules/collier.mod"):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic library")
            config_file = root / "input/mg5_configuration.txt"
            config_file.parent.mkdir()
            config_file.write_text(f"ninja = {ninja}\ncollier = {collier}\n")
            self.assertEqual(runtime.validate_libraries(root), dict(ninja=ninja.resolve(), collier=collier.resolve()))
            (collier / "modules/collier.mod").unlink()
            with self.assertRaisesRegex(ValueError, "COLLIER"):
                runtime.validate_libraries(root)
        sample = ho.all_samples()[0]
        cfg = config()
        self.assertEqual(four.shower_work_dir(sample, cfg), four.shower_work_dir(sample, replace(cfg, detector_response="ssc")))
        self.assertNotEqual(four.sample_work_dir(sample, cfg), four.sample_work_dir(sample, replace(cfg, detector_response="ssc")))
        legacy = replace(cfg, preset="standard")
        self.assertNotEqual(four.shower_work_dir(sample, legacy), four.shower_work_dir(sample, replace(legacy, detector_response="ssc")))

    def test_powheg_seeds_independent_between_tags_and_checked_on_resume(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = presets.seed_base("first")
            powheg.copy_or_make_seeds(root, root, 40, False, False, first)
            before = (root / "pwgseeds.dat").read_bytes()
            self.assertEqual(len(set(before.split())), 40)
            powheg.copy_or_make_seeds(root, root, 40, True, False, first)
            with self.assertRaisesRegex(ValueError, "seeds differ"):
                powheg.copy_or_make_seeds(root, root, 40, True, False, presets.seed_base("second"))
            self.assertEqual((root / "pwgseeds.dat").read_bytes(), before)

    def test_shared_shower_reuses_verified_roots_without_rewriting_evidence(self):
        sample = ho.all_samples()[0]
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)), \
             patch.object(completion_fixtures, "SAMPLE", sample.name):
            cfg = replace(config(), dry_run=False, nevents=3, sample_nevents={sample.name: 3})
            work = four.shower_work_dir(sample, cfg)
            helper = completion_fixtures.ShowerCompletionTests()
            fixture = helper.fixture(work)
            directory, source = fixture[:2]
            card = directory / (sample.name + ".in")
            card.write_text(ho.herwig_card(sample, source, four.seed_for(sample.name + ":herwig", cfg), 3, cfg))
            record, _ = helper.build(fixture)
            record["input_lhe_sha256"] = four.sha256_file(source)
            record["fingerprint"] = ho.completion.completion_fingerprint(record)
            marker = directory / "matched_completion.json"
            marker.write_text(json.dumps(record))
            before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (card, source, marker)}
            with patch.object(ho, "logged_command", side_effect=AssertionError("new shower launched")):
                perfect, _ = ho.run_shower(sample, source, cfg)
                ssc, _ = ho.run_shower(sample, source, replace(cfg, detector_response="ssc"))
                self.assertEqual(perfect, ssc)
                self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before})
                source.write_text(source.read_text() + "\n")
                with self.assertRaisesRegex(ValueError, "shared shower LHE changed"):
                    ho.run_shower(sample, source, cfg)

    def test_analysis_uses_saved_preset_statistical_thresholds(self):
        from types import SimpleNamespace
        for action, threshold in (("cuts", 400.), ("xgboost", 100.)):
            with tempfile.TemporaryDirectory() as tmp, patch.object(cli, "REPO", Path(tmp)):
                root = Path(tmp)
                options = cli.parse_args([action, "--workflow", "fourlepton-ho", "--run-tag", "completed",
                    "--energy-tev", "13", "--work-dir", str(root / "work")])
                source = root / f"hfourlepton/analysis_cards/lhc13_{action}.json"
                source.parent.mkdir(parents=True)
                source.write_text('{"min_background_effective_events":25}')
                cli.shared.save(cli.state_path(options), {"identity": {"sampling": presets.description(args())}})
                commands = []
                pipe = SimpleNamespace(args=options, python=sys.executable, call=lambda label, command: commands.append(command))
                with patch.object(cli, "physics_artifacts"):
                    cli.downstream_analysis(pipe)
                command = commands[0]
                card = Path(command[command.index("--config") + 1])
                self.assertEqual(json.loads(card.read_text())["min_background_effective_events"], threshold)
                options.config = source
                with patch.object(cli, "physics_artifacts"):
                    cli.downstream_analysis(pipe)
                self.assertEqual(commands[-1][commands[-1].index("--config") + 1], source)

    def test_postprocessing_workers_are_bounded_and_merge_private_journals(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            root = Path(tmp)
            cfg = replace(config("--nb-core", "2"), dry_run=False)
            samples = ho.all_samples()[-3:]
            manifest = {"configuration_fingerprint": "test", "samples": [{"name": s.name, "status": "planned"} for s in samples]}
            worker = root / "standin.py"
            worker.write_text("""import json,sys,time
from pathlib import Path
job=Path(sys.argv[1]); payload=json.loads(job.read_text()); name=payload['sample']
manifest=payload['manifest']; journal=Path(payload['config']['manifest_override']['__path__'])
timeline=job.parent/'timeline'
def event(phase):
    with timeline.open('a') as stream: stream.write(json.dumps([time.monotonic(),phase,name])+'\\n')
event('start'); time.sleep(.1)
for entry in manifest['samples']:
    if entry['name']==name: entry.update(status='complete',proof=name)
temporary=journal.with_suffix('.tmp'); temporary.write_text(json.dumps(manifest)); temporary.replace(journal)
event('end')
""")
            popen = subprocess.Popen
            def launch(command, **kwargs):
                self.assertTrue(kwargs["start_new_session"])
                self.assertEqual(kwargs["env"]["OMP_NUM_THREADS"], "1")
                return popen([sys.executable, worker, command[-1]], **kwargs)
            with patch.object(four, "obtain_production_lhe") as production, patch.object(runtime.subprocess, "Popen", side_effect=launch):
                runtime.execute_campaign(samples, manifest, cfg)
            self.assertEqual(production.call_count, 3)
            self.assertTrue(all(s["status"] == "complete" and s["proof"] == s["name"] for s in manifest["samples"]))
            timeline = sorted(json.loads(line) for line in (cfg.run_dir / "workers/perfect/timeline").read_text().splitlines())
            active = maximum = 0
            for _, phase, _ in timeline:
                active += 1 if phase == "start" else -1
                maximum = max(maximum, active)
            self.assertEqual(maximum, 2)
            self.assertEqual(active, 0)
            self.assertEqual(json.loads(cfg.manifest_path.read_text())["samples"], manifest["samples"])

    def test_worker_failure_is_journaled_and_production_cannot_restart(self):
        with tempfile.TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            root = Path(tmp)
            cfg = replace(config("--run-samples", "bkg_gg_4e"), dry_run=False, postprocess_only=True,
                          manifest_override=root / "worker.json")
            manifest = {"samples": [{"name": "bkg_gg_4e", "status": "planned"}]}
            job = root / "job.json"
            job.write_text(json.dumps(dict(config=asdict(cfg), sample="bkg_gg_4e", manifest=manifest), default=runtime._encode))
            with patch.object(four, "run_sample", side_effect=ValueError("test failure")), self.assertRaisesRegex(ValueError, "test failure"):
                runtime.run_worker(job)
            saved = json.loads(cfg.manifest_path.read_text())["samples"][0]
            self.assertEqual(saved["status"], "failed")
            self.assertEqual(saved["error"], "test failure")
            sample = ho.selected_samples(cfg)[0]
            with patch.object(four, "run_runtime_command", side_effect=AssertionError("generator launched")), \
                 self.assertRaisesRegex(ValueError, "production is not complete"):
                four.obtain_mg5_lhe(sample, [sample], cfg)


if __name__ == "__main__":
    unittest.main()
