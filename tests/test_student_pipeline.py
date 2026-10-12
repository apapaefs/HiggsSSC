"""Public student workflow contracts, without running external generators."""
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest import mock

import run_pipeline as cli


class StudentPipelineTests(unittest.TestCase):
    def args(self, root, *extra):
        return cli.parse_args(["run", "--workflow", "gammagamma-lo", "--run-tag", "pilot",
                               "--work-dir", str(root / "work"), *extra])

    @contextlib.contextmanager
    def repository(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()
            with mock.patch.object(cli, "REPO", root), mock.patch.object(cli.shared, "REPO", root):
                yield root

    def artifacts(self, args):
        directory, report = cli.locations(args)
        for sample in cli.shared.SAMPLES:
            category = "Signal" if sample.startswith("signal") else "Backgrounds"
            folder = directory / category / "events" / sample
            cli.shared.save(folder / "campaign.json", {"sample": sample, "completed": ["generate", "shower", "analyze"]})
            for suffix in (".dat", ".top", "_var.root"):
                (folder / f"{sample}_hwsim_roots-{args.run_tag}{suffix}").write_text("synthetic output")
        report.mkdir(parents=True)
        (report / "index.html").write_text("synthetic report")

    def test_help_and_all_workflow_plans_without_site_packages(self):
        repo = Path(cli.__file__).parent
        for options in (["--help"], *[["run", "--workflow", w, "--run-tag", "pilot", "--dry-run"] for w in cli.WORKFLOWS]):
            result = subprocess.run([sys.executable, "-S", str(repo / "run_pipeline.py"), *options], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_dry_run_does_not_initialize_installation(self):
        with self.repository() as root, mock.patch.object(cli, "StudentPipeline") as constructor:
            work = root / "absent"
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main(["setup", "--workflow", "gammagamma-lo", "--work-dir", str(work), "--dry-run"]), 0)
            constructor.assert_not_called()
            self.assertFalse(work.exists())

    def test_explicit_tags_and_workflow_supported_options(self):
        for options in (["run", "--workflow", "gammagamma-lo"],
                        ["cuts", "--workflow", "fourlepton-lo", "--run-tag", "pilot"],
                        ["run", "--workflow", "gammagamma-lo", "--run-tag", "pilot", "--signal-lhe", "events.lhe"]):
            with self.subTest(options=options), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                cli.parse_args(options)

    def test_missing_input_fails_before_constructing_runtime(self):
        with self.repository() as root, mock.patch.object(cli, "StudentPipeline") as constructor:
            with self.assertRaisesRegex(ValueError, "missing or empty"):
                cli.main(["run", "--workflow", "gammagamma-ho", "--run-tag", "pilot",
                          "--work-dir", str(root / "work"), "--signal-lhe", str(root / "missing")])
            constructor.assert_not_called()
            self.assertFalse((root / "work").exists())

    def test_run_identity_rejects_changed_configuration_before_writes(self):
        with self.repository() as root:
            args = self.args(root)
            config = cli.shared.configuration(args)
            cli.prepare_run(args, config)
            before = cli.state_path(args).read_bytes()
            args.nevents = 31
            with self.assertRaisesRegex(ValueError, "configuration changed"):
                cli.prepare_run(args, config)
            self.assertEqual(cli.state_path(args).read_bytes(), before)

    def test_unmanaged_campaign_is_not_adopted(self):
        with self.repository() as root:
            args = self.args(root)
            directory = cli.locations(args)[0]
            directory.mkdir(parents=True)
            saved = directory / "student.root"
            saved.write_bytes(b"preserve")
            with self.assertRaisesRegex(ValueError, "unmanaged"):
                cli.prepare_run(args, {})
            self.assertEqual(saved.read_bytes(), b"preserve")
            self.assertFalse(cli.state_path(args).exists())

    def test_lock_excludes_setup_while_run_uses_installation(self):
        with self.repository() as root:
            with cli.installation_lock(root / "work"):
                with self.assertRaisesRegex(ValueError, "another command"):
                    with cli.installation_lock(root / "work"):
                        self.fail("lock was acquired twice")

    def test_pilot_evidence_survives_report_changes_but_not_physics_changes(self):
        with self.repository() as root:
            args = self.args(root)
            self.artifacts(args)
            cli.shared.save(cli.pilot_path(args), {"environment": {}, "artifacts": cli.physics_artifacts(args)})
            (cli.locations(args)[1] / "index.html").write_text("new presentation")
            cli.verify_pilot(args, {})
            product = next(cli.locations(args)[0].rglob("*.dat"))
            product.write_text("changed physics")
            with self.assertRaisesRegex(ValueError, "physics artifact changed"):
                cli.verify_pilot(args, {})

    def test_missing_sample_prevents_completion(self):
        with self.repository() as root:
            args = self.args(root)
            self.artifacts(args)
            next(cli.locations(args)[0].rglob("*.top")).unlink()
            with self.assertRaisesRegex(ValueError, "analysis product"):
                cli.physics_artifacts(args)

    def test_report_only_never_calls_campaign_or_generator(self):
        with self.repository() as root:
            args = self.args(root)
            self.artifacts(args)
            pipeline = SimpleNamespace(args=args, python=Path("private-python"), call=mock.Mock())
            with contextlib.redirect_stdout(io.StringIO()):
                cli.report(pipeline)
            pipeline.call.assert_called_once()
            command = [str(v) for v in pipeline.call.call_args.args[1]]
            self.assertTrue(command[1].endswith("make_gammagamma_report.py"))
            self.assertEqual(command[command.index("--analysis-root") + 1], str(cli.locations(args)[0]))

    def test_report_refuses_incomplete_analysis_without_external_call(self):
        with self.repository() as root:
            args = self.args(root)
            pipeline = SimpleNamespace(args=args, call=mock.Mock())
            with self.assertRaises(OSError):
                cli.report(pipeline)
            pipeline.call.assert_not_called()

    def test_campaign_commands_use_private_python_and_inherited_environment(self):
        with self.repository() as root:
            args = self.args(root)
            pipeline = SimpleNamespace(args=args, python=root / "venv/bin/python", mg5=root / "mg5")
            command = cli.campaign_command(pipeline, args.workflow)
            self.assertIn("--no-herwig-env", command)
            self.assertIn("--no-herwig-module", command)
            self.assertIn("--resume", command)
            self.assertEqual(command[command.index("--output-dir") + 1], cli.locations(args)[0])

    def test_cuts_and_ml_receive_selected_campaign_paths(self):
        with self.repository() as root:
            args = self.args(root)
            self.artifacts(args)
            args.action = "cuts"
            pipeline = SimpleNamespace(args=args, python=Path("private-python"), call=mock.Mock())
            cli.downstream_analysis(pipeline)
            command = pipeline.call.call_args.args[1]
            self.assertEqual(command[command.index("--analysis-root") + 1], cli.locations(args)[0])
            self.assertEqual(command[command.index("--run-tag") + 1], "pilot")
            self.assertEqual(command[command.index("--luminosity-fb") + 1], 100.)

    def test_status_missing_tag_is_error_and_running_stage_is_not_liveness_claim(self):
        with self.repository() as root:
            args = self.args(root)
            with self.assertRaisesRegex(ValueError, "no managed campaign"):
                cli.status(args)
            cli.prepare_run(args, {})
            cli.progress(args, state="running", stage="ihixs-calculate")
            text = io.StringIO()
            with contextlib.redirect_stdout(text):
                cli.status(args)
            self.assertIn("ihixs-calculate", text.getvalue())
            self.assertIn("not a liveness check", text.getvalue())


if __name__ == "__main__":
    unittest.main()
