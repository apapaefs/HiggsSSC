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
            folder.mkdir(parents=True, exist_ok=True)
            outputs = []
            for suffix in (".dat", ".top", "_var.root"):
                path = folder / f"{sample}_hwsim_roots-{args.run_tag}{suffix}"
                path.write_text("synthetic output")
                outputs.append({"path": str(path), "sha256": cli.shared.digest(path)})
            manifest = {"sample": sample, "completed": ["generate", "shower", "analyze"],
                        "artifacts": {"analyze": outputs}}
            if args.workflow == "gammagamma-ho":
                # Existing HO manifests have .dat values, not LO artifact records.
                manifest.pop("artifacts")
                manifest["analysis"] = {"events_read": "30", "sum_weight": "1"}
                Path(outputs[0]["path"]).write_text("# summary\nevents_read 30\nsum_weight 1\n")
                (folder / f"normalization-{args.run_tag}.json").write_text("{}")
            cli.shared.save(folder / "campaign.json", manifest)
        report.mkdir(parents=True)
        (report / "index.html").write_text("synthetic report")

    def fourlepton_artifacts(self, args):
        directory, report = cli.locations(args)
        for profile in ("perfect", "ssc"):
            samples = []
            for index in range(9):
                folder = directory / "analysis" / profile / str(index)
                folder.mkdir(parents=True)
                root = folder / "analysis.root"
                summary = folder / "analysis.summary.json"
                root.write_text("synthetic ROOT")
                summary.write_text("synthetic summary")
                samples.append({"name": str(index), "status": "complete", "root_file": str(root),
                                "analysis_summary": str(summary), "output_root_sha256": cli.shared.digest(root),
                                "analysis_summary_sha256": cli.shared.digest(summary)})
            cli.shared.save(directory / f"manifest_{profile}.json", {"samples": samples})
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

    def test_complete_run_rejects_changed_outputs_before_constructing_runtime(self):
        with self.repository() as root, mock.patch.object(cli, "StudentPipeline") as constructor:
            args = self.args(root)
            config = cli.shared.configuration(args)
            cli.shared.save(args.work_dir / "environment.json", config)
            cli.prepare_run(args, config)
            self.artifacts(args)
            cli.progress(args, state="complete", stage="complete")
            next(cli.locations(args)[0].rglob("*_var.root")).write_text("changed nonempty output")
            with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                cli.main(["run", "--workflow", args.workflow, "--run-tag", args.run_tag,
                          "--work-dir", str(args.work_dir)])
            constructor.assert_not_called()

    def test_ho_run_checks_matching_pilot_before_children_even_after_report_failure(self):
        for legacy in (False, True):
            with self.subTest(legacy=legacy), self.repository() as root:
                args = self.args(root, "--workflow", "gammagamma-ho")
                self.artifacts(args)
                original = cli.physics_artifacts(args)
                snapshot = args.work_dir / "pilot-complete.json" if legacy else cli.pilot_path(args)
                cli.shared.save(snapshot, {"environment": {}, "artifacts": original})
                cli.shared.save(cli.state_path(args), {"state": "failed", "stage": "pilot-report"})
                before = snapshot.read_bytes()
                next(cli.locations(args)[0].rglob("*_var.root")).write_text("changed nonempty output")
                pipeline = SimpleNamespace(args=args, config={}, run=mock.Mock())
                with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                    cli.run(pipeline)
                pipeline.run.assert_not_called()
                self.assertEqual(snapshot.read_bytes(), before)

    def test_partial_run_ignores_other_tags_pilot_and_can_resume(self):
        with self.repository() as root:
            args = self.args(root, "--workflow", "gammagamma-ho")
            cli.shared.save(cli.state_path(args), {"state": "failed", "stage": "pilot-powheg"})
            unrelated = str(cli.locations(args)[0].parent / "other_tag" / "missing.root")
            for snapshot in (cli.pilot_path(args), args.work_dir / "pilot-complete.json"):
                cli.shared.save(snapshot, {"artifacts": {unrelated: "0" * 64}})
            calls = []
            pipeline = SimpleNamespace(args=args, config={}, run=lambda: calls.append("run"))
            def completed_artifacts(_):
                calls.append("validate")
                return {}
            with mock.patch.object(cli, "physics_artifacts", side_effect=completed_artifacts), \
                    contextlib.redirect_stdout(io.StringIO()):
                cli.run(pipeline)
            self.assertEqual(calls, ["run", "validate"])

    def test_ho_production_saves_per_campaign_artifact_evidence(self):
        with self.repository() as root:
            args = self.args(root, "--workflow", "gammagamma-ho", "--nevents", "100")
            self.artifacts(args)
            cli.shared.save(cli.state_path(args), {"state": "prepared", "stage": "preflight"})
            original = cli.physics_artifacts(args)
            pipeline = SimpleNamespace(args=args, config={}, run=mock.Mock())
            with contextlib.redirect_stdout(io.StringIO()):
                cli.run(pipeline)
            saved = cli.read_state(args)
            self.assertEqual(saved["state"], "complete")
            self.assertEqual(saved["analysis_artifacts"], original)
            self.assertFalse(cli.pilot_path(args).exists())
            self.assertFalse((args.work_dir / "pilot-complete.json").exists())

    def test_ho_production_evidence_survives_report_failure_and_blocks_changed_outputs(self):
        for action in ("run", "report", "cuts", "xgboost"):
            with self.subTest(action=action), self.repository() as root:
                args = self.args(root, "--workflow", "gammagamma-ho", "--nevents", "100")
                self.artifacts(args)
                original = cli.physics_artifacts(args)
                cli.shared.save(cli.state_path(args), {
                    "state": "complete", "stage": "complete", "analysis_artifacts": original,
                })
                cli.progress(args, state="failed", stage="production-report", error="presentation failed")
                product = next(cli.locations(args)[0].rglob("*_var.root"))
                product.write_text("changed nonempty production output")
                args.action = action
                pipeline = SimpleNamespace(args=args, config={}, run=mock.Mock(), call=mock.Mock())
                operation = cli.run if action == "run" else cli.report if action == "report" else cli.downstream_analysis
                with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                    operation(pipeline)
                pipeline.run.assert_not_called()
                pipeline.call.assert_not_called()
                self.assertEqual(cli.read_state(args)["analysis_artifacts"], original)

    def test_saved_campaign_artifact_inventory_cannot_silently_lose_an_output(self):
        with self.repository() as root:
            args = self.args(root, "--workflow", "gammagamma-ho", "--nevents", "100")
            self.artifacts(args)
            original = cli.physics_artifacts(args)
            original.pop(next(name for name in original if name.endswith("_var.root")))
            cli.shared.save(cli.state_path(args), {"state": "failed", "analysis_artifacts": original})
            with self.assertRaisesRegex(ValueError, "inventory differs"):
                cli.validate_saved_run_artifacts(args)

    def test_lo_report_cuts_and_ml_reject_changed_completed_outputs(self):
        for action in ("report", "cuts", "xgboost"):
            for suffix in (".dat", "_var.root"):
                with self.subTest(action=action, suffix=suffix), self.repository() as root:
                    args = self.args(root)
                    args.action = action
                    self.artifacts(args)
                    product = next(cli.locations(args)[0].rglob("*" + suffix))
                    product.write_text("changed but nonempty")
                    pipeline = SimpleNamespace(args=args, call=mock.Mock())
                    with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                        (cli.report if action == "report" else cli.downstream_analysis)(pipeline)
                    pipeline.call.assert_not_called()

    def test_fourlepton_report_checks_saved_root_and_summary_hashes(self):
        for suffix in (".root", ".summary.json"):
            with self.subTest(suffix=suffix), self.repository() as root:
                args = self.args(root, "--workflow", "fourlepton-lo")
                self.fourlepton_artifacts(args)
                self.assertEqual(len(cli.physics_artifacts(args)), 38)
                product = next(cli.locations(args)[0].rglob("*" + suffix))
                product.write_text("changed but nonempty")
                pipeline = SimpleNamespace(args=args, call=mock.Mock())
                with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                    cli.report(pipeline)
                pipeline.call.assert_not_called()

    def test_ho_legacy_summary_contract_and_matching_pilot_snapshot(self):
        with self.repository() as root:
            args = self.args(root, "--workflow", "gammagamma-ho")
            self.artifacts(args)
            original = cli.physics_artifacts(args)
            # An unrelated pilot may no longer exist; it must not block this report.
            cli.shared.save(cli.pilot_path(args), {"artifacts": {str(root / "other_run/missing.root"): "0" * 64}})
            self.assertEqual(cli.physics_artifacts(args), original)
            cli.shared.save(args.work_dir / "pilot-complete.json", {"artifacts": original})
            product = next(cli.locations(args)[0].rglob("*_var.root"))
            product.write_text("changed but nonempty")
            with self.assertRaisesRegex(ValueError, "completed analysis product changed"):
                cli.physics_artifacts(args)

    def test_ho_legacy_changed_summary_is_rejected_without_historical_hashes(self):
        with self.repository() as root:
            args = self.args(root, "--workflow", "gammagamma-ho")
            self.artifacts(args)
            product = next(cli.locations(args)[0].rglob("*.dat"))
            product.write_text("events_read 30\nsum_weight 2\n")
            with self.assertRaisesRegex(ValueError, "HO analysis summary differs"):
                cli.physics_artifacts(args)

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
