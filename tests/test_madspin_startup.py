"""Real MG5 startup and synthetic recovery checks; no decay generation."""
from copy import deepcopy
from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from hfourlepton import madspin_worker
from hfourlepton import repair_madspin_startup as repair
from hfourlepton import run_four_lepton_campaign as four

REPO = Path(__file__).resolve().parents[1]
MG5 = REPO / "MG5_aMC_v3_5_15"
WORKER = Path(madspin_worker.__file__).resolve()


class MadSpinStartupTests(unittest.TestCase):
    def startup(self, worker, root, *extra):
        return subprocess.run([sys.executable, str(worker), "--mg5-dir", str(MG5), "--cores", "5", *extra],
            cwd=root, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, capture_output=True, text=True, timeout=30)

    def test_real_mg5_startup_and_private_cores_without_events_or_config_edits(self):
        config = MG5 / "input/mg5_configuration.txt"
        before = config.read_bytes(), config.stat().st_mtime_ns
        with tempfile.TemporaryDirectory() as tmp:
            result = self.startup(WORKER, tmp, "--check")
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("5 private cores; no events generated", result.stdout)
            self.assertEqual(list(Path(tmp).iterdir()), [])
        self.assertEqual(before, (config.read_bytes(), config.stat().st_mtime_ns))

    def test_missing_formatter_import_reproduces_the_reported_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            worker = Path(tmp) / "old_startup.py"
            worker.write_text(WORKER.read_text().replace("    import madgraph.interface.coloring_logging  # noqa: F401\n", ""))
            result = self.startup(worker, tmp, "--check")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("No module named 'logging.ColorFormatter'", result.stderr)

    def test_launch_applies_private_cores_after_banner_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            harness = Path(tmp) / "startup.py"
            harness.write_text(f"""import sys
sys.path[:0] = [{str(MG5)!r}, {str(REPO)!r}]
from MadSpin.interface_madspin import MadSpinInterface
from hfourlepton import madspin_worker
def check_launch(self, line):
    assert int(self.mg5cmd.options['nb_core']) == 5
    assert int(self.mg5cmd.options['run_mode']) == 2
    assert self.mg5cmd.options['automatic_html_opening'] is False
    print('PRIVATE_LAUNCH_OK')
def fake_import(self, path):
    self.mg5cmd.options.update(nb_core=40, run_mode=0, automatic_html_opening=True)
    self.do_launch('')
MadSpinInterface.do_launch = check_launch
MadSpinInterface.import_command_file = fake_import
madspin_worker.main(['--mg5-dir', {str(MG5)!r}, '--cores', '5', 'unused-card.dat'])
""")
            result = subprocess.run([sys.executable, harness], cwd=tmp,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("PRIVATE_LAUNCH_OK", result.stdout)

    def test_check_and_card_are_mutually_exclusive(self):
        with tempfile.TemporaryDirectory() as tmp:
            for extra in ((), ("--check", "unused-card.dat")):
                result = self.startup(WORKER, tmp, *extra)
                self.assertEqual(result.returncode, 2)
                self.assertIn("provide either", result.stderr)


class StartupRepairTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.patches = [patch.object(four, "SCRIPT_DIR", self.root),
                       patch.dict(os.environ, {}, clear=True)]
        for item in self.patches:
            item.start()
            self.addCleanup(item.stop)
        self.cfg = replace(four.parse_config(["--order", "ho", "--preset", "minimal", "--run-tag", "startup_test",
            "--run-samples", "signal_gg_h_4e", "--ebeam", "6500", "--no-herwig-env", "--no-herwig-module",
            "--dry-run"]), dry_run=False)
        self.samples = four.selected_samples(self.cfg)
        self.current = {"madspin_worker": {"path": str(WORKER), "sha256": repair.FIXED_WORKER},
                        "runner": {"sha256": "unchanged-runner"}}
        item = patch.object(four, "campaign_tracked_inputs", return_value=self.current)
        item.start()
        self.addCleanup(item.stop)
        previous = deepcopy(self.current)
        previous["madspin_worker"]["sha256"] = repair.OLD_WORKER
        self.versions = dict(python="3.9.25", mg5="3.5.15", herwig="7.3.0", root="6.40.04")
        self.manifest = {"status": "running", "tracked_inputs": previous,
            "software": {name + "_version": value for name, value in self.versions.items()},
            "configuration_fingerprint": four.configuration_fingerprint(self.cfg, self.samples,
                tracked_inputs=previous, software_versions=self.versions),
            "samples": [{"name": self.samples[0].name, "status": "failed", "error": "logging startup"}]}
        self.cfg.manifest_path.parent.mkdir(parents=True)
        four.write_manifest(self.manifest, self.cfg)
        self.workers = self.cfg.run_dir / "workers/perfect"
        self.workers.mkdir(parents=True)
        name = self.samples[0].name
        (self.workers / (name + ".json")).write_text(json.dumps(self.manifest))
        (self.workers / (name + ".log")).write_text("logging.config.fileConfig\nNo module named 'logging.ColorFormatter'\n")
        self.decay = four.sample_decay_dir(self.samples[0], self.cfg)
        self.decay.mkdir(parents=True)
        for filename in ("stable_input_source.lhe", "stable_input.lhe", "madspin_card.dat"):
            (self.decay / filename).write_text("synthetic staged input")
        (self.decay / "madspin_card.dat").write_text(four.madspin_card_text(
            self.decay / "stable_input.lhe", self.samples[0], four.seed_for(name + ":madspin", self.cfg)))
        (self.decay / "fourlepton_madspin_banner_stage.json").write_text(json.dumps({
            "artifact": four.artifact_stage_record(self.decay / "stable_input.lhe"),
            "identity": {"source_sha256": four.sha256_file(self.decay / "stable_input_source.lhe")}}))
        self.state_path = self.root / "managed.json"
        self.state_path.write_text('{"state":"failed"}')

    def test_audit_and_apply_preserve_inputs_and_archive_worker_evidence(self):
        original = self.cfg.manifest_path.read_bytes()
        inputs = {path: path.read_bytes() for path in self.decay.iterdir()}
        with patch.object(repair, "completed_production", return_value={"sha256": "verified-hard-events"}):
            report = repair.repair({"state": "failed"}, self.state_path, self.cfg)
            self.assertFalse(report["generators_launched"])
            self.assertEqual(self.cfg.manifest_path.read_bytes(), original)
            self.assertTrue(self.workers.is_dir())
            repair.repair({"state": "failed"}, self.state_path, self.cfg, apply=True)
        backup = self.cfg.run_dir / "recovery" / repair.MIGRATION
        self.assertEqual((backup / "manifest_perfect.json").read_bytes(), original)
        self.assertTrue((backup / "workers/perfect").is_dir())
        self.assertFalse(self.workers.exists())
        self.assertEqual(inputs, {path: path.read_bytes() for path in inputs})
        updated = json.loads(self.cfg.manifest_path.read_text())
        self.assertNotEqual(updated["configuration_fingerprint"], self.manifest["configuration_fingerprint"])
        self.assertEqual(updated["samples"], self.manifest["samples"])
        before = self.cfg.manifest_path.read_bytes()
        repair.repair({"state": "failed"}, self.state_path, self.cfg, apply=True)
        self.assertEqual(self.cfg.manifest_path.read_bytes(), before)

    def test_repair_rejects_other_software_or_configuration_changes(self):
        with patch.object(repair, "completed_production", side_effect=AssertionError("audit should reject first")):
            self.current["runner"]["sha256"] = "unexpected-runner"
            with self.assertRaisesRegex(ValueError, "beyond the MadSpin"):
                repair.repair({"state": "failed"}, self.state_path, self.cfg)
            self.current["runner"]["sha256"] = "unchanged-runner"
            with self.assertRaisesRegex(ValueError, "beyond the MadSpin"):
                repair.repair({"state": "failed"}, self.state_path, replace(self.cfg, ebeam=20000))

    def test_repair_rejects_decay_results_and_changed_staged_inputs(self):
        output = self.decay / "stable_input_decayed.lhe.gz"
        output.write_bytes(b"existing decay result")
        with self.assertRaisesRegex(ValueError, "progressed beyond startup"):
            repair.repair({"state": "failed"}, self.state_path, self.cfg)
        output.unlink()
        (self.decay / "stable_input.lhe").write_text("changed input")
        with self.assertRaisesRegex(ValueError, "staged decay input changed"):
            repair.repair({"state": "failed"}, self.state_path, self.cfg)

    def test_metadata_write_failure_restores_worker_directory(self):
        original = self.cfg.manifest_path.read_bytes()
        with patch.object(repair, "completed_production", return_value={}), \
             patch.object(four, "write_manifest", side_effect=OSError("test write failure")):
            with self.assertRaisesRegex(OSError, "test write failure"):
                repair.repair({"state": "failed"}, self.state_path, self.cfg, apply=True)
        self.assertTrue(self.workers.is_dir())
        self.assertEqual(self.cfg.manifest_path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
