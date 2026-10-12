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
        with TemporaryDirectory() as tmp, mock.patch.object(pipeline.norm, "load_record"):
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
