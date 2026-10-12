"""Pipeline orchestration tests; external generators use synthetic boundaries.

These tests do not claim physics validation. The existing generator, ihixs,
shower and normalization tests own validation of the actual file formats.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

from hgammagamma import run_gammagamma_ho_pipeline as pipeline


class PipelineTests(unittest.TestCase):
    def fixture(self, root, *extra):
        args = pipeline.parse_args(["run", "--work-dir", str(root / "work"), *extra])
        config = pipeline.configuration(args)
        obj = pipeline.Pipeline.__new__(pipeline.Pipeline)
        obj.args, obj.config, obj.work = args, config, args.work_dir
        obj.prefix, obj.boost = root / "runtime", root / "runtime/include"
        obj.python, obj.mg5 = root / "python", root / "mg5"
        obj.cuba, obj.powheg = root / "cuba", root / "powheg"
        obj.check = mock.Mock()
        pipeline.save(args.work_dir / "environment.json", config)
        calls = []

        def fake_call(label, command, **kwargs):
            command = [str(x) for x in command]
            calls.append((label, command))
            if label == "ihixs-calculate":
                record = Path(command[command.index("--record") + 1])
                record.parent.mkdir(parents=True, exist_ok=True)
                record.write_text("SYNTHETIC record: validated only by the mocked boundary")
            elif label.endswith("-powheg"):
                folder = Path(command[command.index("--run-dir") + 1])
                folder.mkdir(parents=True, exist_ok=True)
                (folder / "powheg-hjminnlo-merged.lhe").write_text("SYNTHETIC LHE")
            elif label.endswith("-campaign"):
                directory = pipeline.HO / "runs" / args.run_tag
                for sample in pipeline.SAMPLES:
                    category = "Signal" if sample.startswith("signal_") else "Backgrounds"
                    folder = directory / category / "events" / sample
                    pipeline.save(folder / "campaign.json", {"sample": sample,
                                  "completed": ["generate", "shower", "analyze"]})
                    for suffix in (".dat", ".top", "_var.root"):
                        (folder / f"{sample}_hwsim_roots-{args.run_tag}{suffix}").write_text("SYNTHETIC analysis")
                    pipeline.save(folder / f"normalization-{args.run_tag}.json", {"synthetic": True})
            elif label.endswith("-report"):
                folder = Path(command[command.index("--output-dir") + 1])
                folder.mkdir(parents=True, exist_ok=True)
                (folder / "index.html").write_text("SYNTHETIC report")

        obj.call = fake_call
        return obj, calls

    def ho_fixture(self, root):
        ho = root / "HO"
        ho.mkdir()
        (ho / "powheg-hjminnlo-ssc40-nnpdf40nnloqed.input").write_text(
            "ncall2 50000\ncheck_bad_st2 1\n")
        return mock.patch.object(pipeline, "HO", ho)

    def test_fresh_pipeline_calculates_normalization_before_events_and_report(self):
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline.norm, "load_record") as validate:
            root = Path(tmp)
            with self.ho_fixture(root), contextlib.redirect_stdout(io.StringIO()):
                obj, calls = self.fixture(root)
                obj.run()
                self.assertEqual([label for label, _ in calls], [
                    "ihixs-build", "ihixs-benchmark", "ihixs-calculate",
                    "ho_pilot_01-powheg", "ho_pilot_01-campaign", "ho_pilot_01-report"])
                validate.assert_called_once()
                self.assertIn("--refine-failed", calls[2][1])
                card = (pipeline.HO / "runs/ho_pilot_01/powheg-pipeline.input").read_text()
                self.assertIn("ncall2 500000", card)
                self.assertIn("check_bad_st2 1", card)
                pipeline.verify_pilot(obj.work)

    def test_existing_record_and_lhe_skip_expensive_upstream_stages(self):
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline.norm, "load_record"), \
                mock.patch.object(pipeline.ho_campaign, "inspect_lhe") as inspect, \
                mock.patch.object(pipeline.ho_campaign.shower_completion, "source_weight_summary") as weights:
            root = Path(tmp)
            record, lhe = root / "rate.json", root / "events.lhe"
            record.write_text("SYNTHETIC record")
            lhe.write_text("SYNTHETIC LHE")
            with self.ho_fixture(root), contextlib.redirect_stdout(io.StringIO()):
                obj, calls = self.fixture(root, "--signal-lhe", str(lhe), "--signal-normalization", str(record))
                obj.run()
                self.assertEqual([label for label, _ in calls], ["ho_pilot_01-campaign", "ho_pilot_01-report"])
                campaign = calls[0][1]
                self.assertEqual(campaign[campaign.index("--signal-lhe") + 1], str(lhe.resolve()))
                inspect.assert_called_once()
                weights.assert_called_once_with(lhe.resolve(), 30)

    def test_supplied_missing_lhe_prevents_ihixs_and_dependency_work(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, calls = self.fixture(root, "--signal-lhe", str(root / "missing.lhe"))
            with self.assertRaisesRegex(ValueError, "supplied signal LHE is missing"):
                obj.run()
            obj.check.assert_not_called()
            self.assertEqual(calls, [])

    def test_supplied_invalid_lhe_prevents_ihixs_and_dependency_work(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            lhe = root / "incomplete.lhe"
            lhe.write_text("<LesHouchesEvents>")
            obj, calls = self.fixture(root, "--signal-lhe", str(lhe))
            with self.assertRaisesRegex(ValueError, "empty or incomplete LHE"):
                obj.run()
            obj.check.assert_not_called()
            self.assertEqual(calls, [])

    def test_invalid_record_prevents_generation(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            record = root / "bad.json"
            record.write_text("{}")
            with self.ho_fixture(root):
                obj, calls = self.fixture(root, "--signal-normalization", str(record))
                with self.assertRaises(ValueError):
                    obj.run()
                self.assertEqual(calls, [])

    def test_ihixs_failure_prevents_generation_and_preserves_run_identity(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.ho_fixture(root):
                obj, _ = self.fixture(root)
                obj.call = mock.Mock(side_effect=RuntimeError("synthetic build failure"))
                with self.assertRaisesRegex(RuntimeError, "synthetic build failure"):
                    obj.run()
                obj.call.assert_called_once()
                self.assertTrue((pipeline.HO / "runs/ho_pilot_01/pipeline.json").is_file())
                self.assertFalse((obj.work / "pilot-complete.json").exists())

    def test_existing_student_campaign_is_preserved(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.ho_fixture(root):
                obj, calls = self.fixture(root)
                folder = pipeline.HO / "runs/ho_pilot_01"
                folder.mkdir(parents=True)
                saved = folder / "student.root"
                saved.write_bytes(b"saved student events")
                with self.assertRaisesRegex(ValueError, "existing campaign"):
                    obj.run()
                self.assertEqual(saved.read_bytes(), b"saved student events")
                self.assertEqual(calls, [])

    def test_production_requires_completed_unchanged_pilot(self):
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline.norm, "load_record"):
            root = Path(tmp)
            with self.ho_fixture(root), contextlib.redirect_stdout(io.StringIO()):
                obj, calls = self.fixture(root, "--nevents", "10000")
                with self.assertRaisesRegex(ValueError, "first complete a pilot"):
                    obj.run()
                self.assertEqual(calls, [])
                obj.args.nevents = 30
                obj.run()
                pilot = json.loads((obj.work / "pilot-complete.json").read_text())
                Path(next(iter(pilot["artifacts"]))).write_text("changed")
                with self.assertRaisesRegex(ValueError, "pilot artifact changed"):
                    pipeline.verify_pilot(obj.work)

    def test_prefix_activation_can_be_reused_and_changed_compiler_is_rejected(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = root / "runtime"
            prefix.mkdir()
            args = pipeline.parse_args(["setup", "--work-dir", str(root / "work"),
                                        "--herwig-env", str(prefix)])
            config = pipeline.configuration(args)
            pipeline.save(args.work_dir / "environment.json", config)
            self.assertEqual(pipeline.configuration(args), config)
            args.cc = "/different/gcc"
            with self.assertRaisesRegex(ValueError, "cc differs"):
                pipeline.configuration(args)

    def test_child_environment_cannot_reactivate_over_private_pdf_path(self):
        env = {"PATH": "/runtime/bin", "HERWIG_ENV": "/shared", "HERWIG_MODULE": "wrong",
               "MG5_PYTHON": "/other/python", "DEFAULT_HERWIG_ENV": "/other"}
        result = subprocess.CompletedProcess([], 0, "HIGGSSSC_ENV=" + json.dumps(env), "")
        with mock.patch.object(pipeline.subprocess, "run", return_value=result):
            config = pipeline.configuration(pipeline.parse_args(["setup"]))
            actual = pipeline.activated_environment(config)
        self.assertEqual(actual["PATH"], "/runtime/bin")
        self.assertNotIn("HERWIG_ENV", actual)
        self.assertNotIn("HERWIG_MODULE", actual)
        self.assertEqual(actual["FC"], config["fc"])

    def test_dry_run_needs_no_runtime_and_creates_no_files(self):
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline, "Pipeline") as constructor:
            work = Path(tmp) / "not-created"
            with contextlib.redirect_stdout(io.StringIO()):
                pipeline.main(["run", "--work-dir", str(work), "--dry-run"])
            constructor.assert_not_called()
            self.assertFalse(work.exists())

    def test_mg5_startup_fix_keeps_original_and_is_idempotent(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "mg5_aMC"
            original = "import sys\nif sys.version_info < (3, 7):\n    pass\nlogging.warning('test')\n"
            path.write_text(original)
            path.chmod(0o755)
            pipeline.repair_mg5_startup(path)
            once = path.read_text()
            pipeline.repair_mg5_startup(path)
            self.assertEqual(path.read_text(), once)
            self.assertEqual(path.with_name("mg5_aMC.before-pipeline").read_text(), original)
            self.assertEqual(once.count("import logging"), 1)
            self.assertTrue(os.access(path, os.X_OK))

    def test_pdf_install_forces_private_copy_but_preserves_completed_sets(self):
        with TemporaryDirectory() as tmp:
            obj = pipeline.Pipeline.__new__(pipeline.Pipeline)
            obj.pdf_dir = Path(tmp)
            for name in pipeline.PDF_SETS:
                for path in pipeline.pdf_files(obj.pdf_dir, name):
                    path.parent.mkdir(exist_ok=True)
                    path.write_text("synthetic PDF")
            missing = pipeline.pdf_files(obj.pdf_dir, pipeline.PDF_SETS[1])[-1]
            missing.unlink()
            def install(_label, command):
                self.assertEqual(command[-3:], ["install", "--upgrade", pipeline.PDF_SETS[1]])
                missing.write_text("synthetic recovered PDF")
            obj.call = mock.Mock(side_effect=install)
            obj.install_pdfs(["lhapdf", "--pdfdir", obj.pdf_dir])
            obj.call.assert_called_once()
            obj.call.reset_mock()
            obj.install_pdfs(["lhapdf", "--pdfdir", obj.pdf_dir])
            obj.call.assert_not_called()

    def test_lo_setup_then_ho_extension_preserves_installation_and_pilot(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp)
            obj, _ = self.fixture(root)
            obj.args.workflow = "gammagamma-lo"
            obj.pdf_dir = obj.work / "pdfs"
            obj.runtime_check = mock.Mock()
            for name in ("setup_python", "setup_mg5", "setup_cuba", "setup_powheg", "setup_analysis"):
                setattr(obj, name, mock.Mock())
            def fake_pdfs(_command, names):
                for name in names:
                    info = obj.pdf_dir / name / (name + ".info")
                    info.parent.mkdir(parents=True, exist_ok=True)
                    info.write_text("synthetic PDF info")
            obj.install_pdfs = mock.Mock(side_effect=fake_pdfs)
            obj.call = mock.Mock()
            obj.capture = mock.Mock(return_value="synthetic ihixs commit")
            pipeline.save(obj.work / "pilot-complete.json", {"legacy": "untouched"})
            pilot_before = (obj.work / "pilot-complete.json").read_bytes()
            obj.setup()
            obj.install_pdfs.assert_called_with(mock.ANY, pipeline.LO_PDF_SETS)
            obj.setup_cuba.assert_not_called()
            obj.setup_powheg.assert_not_called()
            installed_before = (obj.work / "installed.json").read_bytes()
            obj.args.workflow = "gammagamma-ho"
            obj.setup()
            obj.install_pdfs.assert_called_with(mock.ANY, pipeline.PDF_SETS)
            obj.setup_cuba.assert_called_once()
            obj.setup_powheg.assert_called_once()
            self.assertEqual((obj.work / "installed.json").read_bytes(), installed_before)
            self.assertEqual((obj.work / "pilot-complete.json").read_bytes(), pilot_before)
            self.assertEqual(obj.components()["gammagamma-lo"], {"ready": True})
            self.assertEqual(obj.components()["gammagamma-ho"], {"ready": True})

    def test_optional_ml_is_additive_and_existing_python_packages_are_preserved(self):
        with TemporaryDirectory() as tmp:
            obj, _ = self.fixture(Path(tmp))
            obj.python.write_text("synthetic installed interpreter")
            pipeline.save(obj.work / "installed.json", {"legacy": True})
            obj.call = mock.Mock()
            versions = {"numpy": "1.24.4", "matplotlib": "3.8.4", "PyYAML": "6.0.2", "six": "1.17.0"}
            def capture(command):
                if "-c" in command:
                    compile(command[-1], "dependency-probe", "exec")
                    if command[-1] == pipeline.NUMPY_TYPING_PROBE:
                        return json.dumps({"typing_error": None})
                    return json.dumps(versions)
                return "numpy==1.24.4\n"
            obj.capture = mock.Mock(side_effect=capture)
            obj.setup_python()
            self.assertEqual([call.args[0] for call in obj.call.call_args_list], ["check-python"])
            obj.call.reset_mock()
            obj.args.with_ml = True
            obj.setup_python()
            self.assertEqual([call.args[0] for call in obj.call.call_args_list],
                             ["check-python", "python-ml", "check-python-ml"])
            install = obj.call.call_args_list[1].args[1]
            self.assertEqual(install[-2:], ["--requirement", pipeline.REPO / "requirements-ml.txt"])
            constraints = Path(install[install.index("--constraint") + 1])
            self.assertEqual(constraints.read_text(), "".join(f"{name}=={version}\n" for name, version in versions.items()))
            obj.call.reset_mock()
            obj.setup_python()
            self.assertEqual([call.args[0] for call in obj.call.call_args_list],
                             ["check-python", "check-python-ml"])

    def test_incompatible_optional_ml_keeps_core_and_does_not_record_success(self):
        with TemporaryDirectory() as tmp:
            obj, _ = self.fixture(Path(tmp))
            obj.python.write_text("synthetic installed interpreter")
            obj.args.with_ml = True
            obj.record_component("python_core", {"packages": "saved validated core"})
            obj.ensure_ml_numpy = mock.Mock()
            obj.capture = mock.Mock(return_value=json.dumps({name: "1.0" for name in pipeline.CORE_PACKAGES}))
            def install(label, command):
                if label == "python-ml":
                    self.assertIn("--constraint", command)
                    raise RuntimeError("dependencies cannot be resolved with pinned core versions")
            obj.call = mock.Mock(side_effect=install)
            with self.assertRaisesRegex(RuntimeError, "pinned core"):
                obj.setup_python()
            self.assertNotIn("python_ml", obj.components())
            self.assertEqual(obj.components()["python_core"], {"packages": "saved validated core"})

    def test_ml_repairs_only_incomplete_inherited_numpy_at_the_same_version(self):
        with TemporaryDirectory() as tmp:
            obj, _ = self.fixture(Path(tmp))
            # The real interpreter is inside the venv; only the test fixture
            # normally substitutes a short path for external-command mocks.
            obj.python = obj.work / "venv/bin/python"
            system = {"version": "1.23.5", "distribution_version": "1.23.5",
                      "path": "/system/lib/numpy/__init__.py", "typing_error": "No module named numpy._typing"}
            private = {**system, "path": str(obj.work / "venv/lib/python3.9/site-packages/numpy/__init__.py"),
                       "typing_error": None}
            pipeline.save(obj.work / "installed.json", {"legacy": "unchanged"})
            before = (obj.work / "installed.json").read_bytes()
            obj.capture = mock.Mock(side_effect=[json.dumps(system), json.dumps(private), "numpy==1.23.5\n"])
            obj.call = mock.Mock()
            with contextlib.redirect_stdout(io.StringIO()):
                obj.ensure_ml_numpy()
            command = obj.call.call_args_list[0].args[1]
            self.assertEqual(command, [obj.python, "-m", "pip", "install", "--ignore-installed",
                                       "--no-deps", "numpy==1.23.5"])
            self.assertEqual(obj.components()["numpy_private_repair"],
                             {"version": "1.23.5", "distribution_version": "1.23.5",
                              "from_path": system["path"], "installed_path": private["path"]})
            self.assertEqual(obj.components()["python_core"]["numpy"],
                             {"version": "1.23.5", "path": private["path"]})
            self.assertEqual((obj.work / "installed.json").read_bytes(), before)
            obj.capture = mock.Mock(return_value=json.dumps(private))
            obj.call.reset_mock()
            obj.ensure_ml_numpy()
            obj.call.assert_not_called()

    def test_ml_does_not_overwrite_an_incomplete_private_numpy(self):
        with TemporaryDirectory() as tmp:
            obj, _ = self.fixture(Path(tmp))
            obj.capture = mock.Mock(return_value=json.dumps({"version": "1.23.5", "distribution_version": "1.23.5",
                "path": str(obj.work / "venv/lib/python3.9/site-packages/numpy/__init__.py"),
                "typing_error": "No module named numpy._typing"}))
            obj.call = mock.Mock()
            with self.assertRaisesRegex(ValueError, "preserve this installation"):
                obj.ensure_ml_numpy()
            obj.call.assert_not_called()

    def test_ml_repair_rejects_a_changed_numpy_version(self):
        with TemporaryDirectory() as tmp:
            obj, _ = self.fixture(Path(tmp))
            before = {"version": "1.23.5", "distribution_version": "1.23.5", "path": "/system/numpy/__init__.py",
                      "typing_error": "No module named numpy._typing"}
            after = {**before, "version": "2.0.0", "distribution_version": "2.0.0", "typing_error": None,
                     "path": str(obj.work / "venv/lib/numpy/__init__.py")}
            obj.capture = mock.Mock(side_effect=[json.dumps(before), json.dumps(after)])
            obj.call = mock.Mock()
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "did not pass import checks"):
                obj.ensure_ml_numpy()
            self.assertNotIn("numpy_private_repair", obj.components())

    def test_lo_model_copies_only_tracked_source_and_rejects_changed_installation(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = self.fixture(root)
            source = root / "MG5_aMC_v3_5_15/models/loop_sm_haa"
            source.mkdir(parents=True)
            for name in ("parameters.py", "restrict_default.dat", "py3_model.pkl", "untracked.py"):
                (source / name).write_text("synthetic model source")
            tracked = [source / name for name in ("parameters.py", "restrict_default.dat", "py3_model.pkl")]
            obj.capture = mock.Mock(return_value="\n".join(str(path.relative_to(root)) for path in tracked))
            with mock.patch.object(pipeline, "REPO", root):
                obj.install_lo_model()
                target = obj.mg5 / "models/loop_sm_haa"
                self.assertEqual({path.name for path in target.iterdir()}, {"parameters.py", "restrict_default.dat"})
                obj.install_lo_model()
                (target / "parameters.py").write_text("different saved model")
                with self.assertRaisesRegex(ValueError, "installed LO model differs"):
                    obj.install_lo_model()
                self.assertEqual((target / "parameters.py").read_text(), "different saved model")

    def test_fourlepton_build_selects_its_own_analyzer(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = self.fixture(root)
            obj.args.workflow = "fourlepton-lo"
            target = root / "hfourlepton/LOAnalysis/Code/HwSimPostAnalysis_fourlepton"
            target.parent.mkdir(parents=True)
            target.write_text("synthetic analyzer")
            obj.call = mock.Mock()
            with mock.patch.object(pipeline, "REPO", root):
                obj.setup_analysis()
            self.assertIn("HwSimPostAnalysis_fourlepton", obj.call.call_args.args[1])
            self.assertIn("fourlepton_analysis", obj.components())
            self.assertEqual(obj.pdf_sets, pipeline.LO_PDF_SETS)

    def test_lo_dependency_check_does_not_require_ho_installation(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            root = Path(tmp)
            obj, _ = self.fixture(root)
            del obj.check
            obj.args.workflow = "gammagamma-lo"
            obj.pdf_dir = obj.work / "pdfs"
            paths = [obj.python, obj.mg5 / "bin/mg5_aMC",
                     obj.mg5 / "models/loop_sm_haa/parameters.py", obj.pdf_dir / "pdfsets.index",
                     root / "hgammagamma/LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC"]
            paths.extend(pipeline.pdf_files(obj.pdf_dir, pipeline.LO_PDF_SETS[0]))
            for path in paths:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("synthetic dependency")
            obj.runtime_check = mock.Mock()
            obj.call = mock.Mock()
            with mock.patch.object(pipeline, "REPO", root):
                obj.check()
                obj.args.workflow = "gammagamma-ho"
                with self.assertRaisesRegex(ValueError, "setup is incomplete for gammagamma-ho"):
                    obj.check()

    def test_activated_private_interpreter_is_first_for_child_launchers(self):
        with TemporaryDirectory() as tmp:
            args = pipeline.parse_args(["setup", "--work-dir", str(Path(tmp) / "work")])
            config = pipeline.configuration(args)
            with mock.patch.object(pipeline, "activated_environment", return_value={"PATH": "/shared/bin"}), \
                    mock.patch.object(pipeline.Pipeline, "capture", return_value="/shared"):
                obj = pipeline.Pipeline(args, config)
            self.assertEqual(obj.env["PATH"].split(":")[0], str(obj.python.parent))
            self.assertEqual(obj.env["PYTHONNOUSERSITE"], "1")
            self.assertEqual(obj.env["PYTHONUNBUFFERED"], "1")

    def test_requested_workers_bound_cuba_and_nested_numeric_threads(self):
        with TemporaryDirectory() as tmp:
            args = pipeline.parse_args(["setup", "--work-dir", str(Path(tmp) / "work"), "--jobs", "3"])
            config = pipeline.configuration(args)
            inherited = {"PATH": "/shared/bin", "CUBACORES": "64", "CUBAACCEL": "8",
                         "OMP_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "32", "MKL_NUM_THREADS": "64"}
            with mock.patch.object(pipeline, "activated_environment", return_value=inherited), \
                    mock.patch.object(pipeline.Pipeline, "capture", return_value="/shared"):
                obj = pipeline.Pipeline(args, config)
            self.assertEqual(obj.env["CUBACORES"], "3")
            self.assertEqual(obj.env["CUBAACCEL"], "0")
            for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
                self.assertEqual(obj.env[name], "1")

    def test_recorded_mg5_ninja_precedes_inherited_same_soname_on_new_invocations(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = pipeline.parse_args(["check", "--work-dir", str(root / "work")])
            config = pipeline.configuration(args)
            chosen = root / "quad-ninja/lib"
            chosen.mkdir(parents=True)
            pipeline.save(args.work_dir / "components.json", {"mg5": {"ninja": str(chosen)}})
            private_config = args.work_dir / "MG5_aMC_v3_5_15/input/mg5_configuration.txt"
            private_config.parent.mkdir(parents=True)
            private_config.write_text("ninja = /unrelated/ordinary-ninja\n")
            inherited = {"PATH": "/shared/bin", "LD_LIBRARY_PATH": "/shared/lib:" + str(chosen)}
            with mock.patch.object(pipeline, "activated_environment", side_effect=lambda _: dict(inherited)), \
                    mock.patch.object(pipeline.Pipeline, "capture", return_value="/shared"), \
                    mock.patch.object(pipeline.sys, "platform", "linux"):
                for _ in range(2):
                    obj = pipeline.Pipeline(args, config)
                    self.assertEqual(obj.env["LD_LIBRARY_PATH"], str(chosen) + ":/shared/lib")
                (args.work_dir / "components.json").unlink()
                private_config.write_text("ninja = " + str(chosen) + "\n")
                obj = pipeline.Pipeline(args, config)
                self.assertEqual(obj.env["LD_LIBRARY_PATH"].split(":")[0], str(chosen))

    def test_initial_mg5_setup_activates_and_preserves_its_selected_ninja(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = self.fixture(root)
            obj.args.workflow = "fourlepton-lo"
            obj.env = {"LD_LIBRARY_PATH": "/shared/lib"}
            obj.capture = mock.Mock(return_value="/shared/bin/tool")
            (obj.mg5 / "bin").mkdir(parents=True)
            (obj.mg5 / "bin/mg5_aMC").write_text("synthetic entrypoint")
            (obj.mg5 / "input").mkdir()
            chosen = obj.prefix / "opt/MG5-selected/HEPTools/ninja/lib"
            chosen.mkdir(parents=True)
            (chosen.parent / "include").mkdir()
            for path in (chosen / "libninja.a", chosen / "libavh_olo.a", chosen.parent / "include/mninja.mod"):
                path.write_text("synthetic library")
            with mock.patch.object(pipeline, "prepare_mg5_archive"), \
                    mock.patch.object(pipeline.sys, "platform", "linux"):
                obj.setup_mg5()
                self.assertEqual(obj.env["LD_LIBRARY_PATH"].split(":")[0], str(chosen))
                self.assertEqual(obj.components()["mg5"]["ninja"], str(chosen))
                # A new earlier-named installation must not silently replace
                # the build recorded by the existing private MG5 installation.
                alternative = obj.prefix / "opt/MG5-a-new/HEPTools/ninja/lib"
                alternative.mkdir(parents=True)
                obj.setup_mg5()
                self.assertEqual(obj.components()["mg5"]["ninja"], str(chosen))

    def test_new_pilot_evidence_excludes_report_presentation(self):
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline.norm, "load_record"):
            root = Path(tmp)
            with self.ho_fixture(root), contextlib.redirect_stdout(io.StringIO()):
                obj, _ = self.fixture(root)
                obj.run()
                report = pipeline.HO / "plots/ho_pilot_01/index.html"
                report.write_text("Updated report presentation")
                pipeline.verify_pilot(obj.work)

    def test_tag_archive_uses_release_packaging_and_keeps_development_utility(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "bin").mkdir()
            (root / "VERSION").write_text("version = 3.5.15\ndate = 2026-04-17\n")
            (root / "bin/create_release.py").write_text("synthetic upstream release utility")
            (root / "bin/mg5_aMC").write_text("import sys\nif sys.version_info < (3, 7):\n    pass\n")
            pipeline.prepare_mg5_archive(root)
            pipeline.prepare_mg5_archive(root)
            self.assertFalse((root / "bin/create_release.py").exists())
            self.assertEqual((root / "source-packaging/create_release.py").read_text(),
                             "synthetic upstream release utility")

    def test_mg5_zero_exit_error_is_detected_without_poisoning_next_attempt(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            obj = pipeline.Pipeline.__new__(pipeline.Pipeline)
            obj.env, obj.logs = dict(os.environ), Path(tmp)
            with self.assertRaisesRegex(RuntimeError, "Error detected in"):
                obj.call("mg5", [sys.executable, "-c", "print('Error detected in synthetic MG5')"], mg5=True)
            obj.call("mg5", [sys.executable, "-c", "print('synthetic success')"], mg5=True)
            self.assertIn("Error detected in", (Path(tmp) / "mg5.log").read_text())

    @unittest.skipUnless(hasattr(tarfile, "data_filter"), "security-patched Python required")
    def test_archive_cannot_write_outside_private_staging(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            archive = root / "bad.tar.gz"
            with tarfile.open(archive, "w:gz") as tar:
                entry = tarfile.TarInfo("../../escaped")
                entry.size = 1
                tar.addfile(entry, io.BytesIO(b"x"))
            obj = pipeline.Pipeline.__new__(pipeline.Pipeline)
            obj.work = root / "work"
            with self.assertRaises(tarfile.FilterError):
                obj.unpack(archive, "bad", obj.work / "package")
            self.assertFalse((root / "escaped").exists())


if __name__ == "__main__":
    unittest.main()
