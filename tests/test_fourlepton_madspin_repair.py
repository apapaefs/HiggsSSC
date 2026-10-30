"""Exercise the bounded migration with synthetic files; no simulation is run."""
from dataclasses import replace
import gzip
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hfourlepton import run_four_lepton_campaign as four, higher_order as ho
from hfourlepton import repair_madspin_banner as repair


class MadSpinRepairTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for target, value in (("SCRIPT_DIR", self.root),):
            patcher = patch.object(four, target, value)
            patcher.start()
            self.addCleanup(patcher.stop)
        patches = [
            # The fixture synthesizes both sides of an allowed transition.
            # Keep production's historical hash allowlist fixed; unrelated new
            # runner features must not broaden the real migration's scope.
            patch.object(repair, "FIXED_RUNNER", four.sha256_file(Path(four.__file__))),
            patch.object(four, "campaign_software_versions", return_value={"python": "3.9.25", "mg5": "3.5.15", "herwig": "7.3", "root": "6.40"}),
            patch.object(four, "run_runtime_command", side_effect=AssertionError("no generator may run")),
            patch.object(ho, "validate_nlo_pole_checks"),
            patch.object(ho, "inspect_production", return_value={"events": 3}),
        ]
        for patcher in patches:
            patcher.start()
            self.addCleanup(patcher.stop)
        with patch.dict("os.environ", {}, clear=True):
            self.cfg = replace(four.parse_config(["--order", "ho", "--signal-modes", "all",
                "--run-tag", "repair", "--ebeam", "6500", "--higgs-mass", "125", "--nevents", "1",
                "--detector-response", "perfect", "--selection-profile", "cms", "--no-herwig-module",
                "--mg5-dir", str(self.root / "mg5"), "--dry-run"]), dry_run=False)
        self.samples = four.selected_samples(self.cfg)
        self.manifest = four.build_manifest(self.cfg, self.samples)
        self.manifest["tracked_inputs"]["runner"]["sha256"] = repair.OLD_RUNNER
        self.manifest["software"]["runner_sha256"] = repair.OLD_RUNNER
        self.manifest["configuration_fingerprint"] = four.configuration_fingerprint(
            self.cfg, self.samples, tracked_inputs=self.manifest["tracked_inputs"])
        for sample in self.samples[:6]:
            entry = four.manifest_sample(self.manifest, sample.name)
            entry["status"] = "complete"
            for key, suffix in (("output_root", ".root"), ("analysis_summary", ".summary.json"),
                                ("source_lhe", ".source.lhe"), ("input_lhe", ".input.lhe")):
                path = four.expected_analysis_output(sample, self.cfg).with_suffix(suffix)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps({"events_read": 1, "generated_raw": {}, "generated_scaled": {},
                    "selected_scaled": {}}) if key == "analysis_summary" else key)
                entry[key], entry[key + "_sha256"] = str(path), four.sha256_file(path)
        self.sample = self.samples[6]
        four.manifest_sample(self.manifest, self.sample.name)["status"] = "running"
        four.write_manifest(self.manifest, self.cfg)
        self.original_manifest = self.cfg.manifest_path.read_bytes()

        group = four.production_dir(self.sample, self.cfg)
        process_dir = group / "mg5_process"
        self.source = process_dir / "Events/repair/events.lhe.gz"
        self.source.parent.mkdir(parents=True)
        self.source.write_bytes(b"synthetic sealed production")
        process, launch = ho.nlo_cards(self.sample, self.cfg, process_dir, 3)
        identity = dict(process_card=process, compile_card=ho.nlo_compile_card(launch), launch_card=launch,
            adapter_sha256=four.sha256_file(Path(ho.__file__)),
            guard_sha256=four.sha256_file(Path(ho.__file__).with_name("mg5_guard.py")))
        four.write_artifact_stage_marker(group / "fourlepton_production_stage.json", identity, self.source)
        slha = """BLOCK SMINPUTS
1 132.507
2 1.16639e-5
3 .118
BLOCK MASS
5 0
6 173
15 1.777
23 91.188
24 80.419
25 125
BLOCK YUKAWA
6 173
15 1.777
DECAY 6 1.4915
DECAY 23 2.441404
DECAY 24 2.0476
DECAY 25 .006382339
"""
        raw = ("<LesHouchesEvents>\n<header>\n<MG5ProcCard>\nimport model loop_sm-no_b_mass\n"
               "generate p p > h w+ [QCD]\n</MG5ProcCard>\n<MGRunCard>\n3 = lhe_version\n</MGRunCard>\n"
               f"<slha>\n{slha}</slha>\n</header>\n<init>\n</init>\n</LesHouchesEvents>\n")
        sliced = group / "decay_inputs/4e.lhe.gz"
        sliced.parent.mkdir()
        with gzip.open(sliced, "wt") as stream:
            stream.write(raw)
        four.write_artifact_stage_marker(sliced.with_suffix(".slice.json"),
            dict(kind="disjoint_hard_event_slice", source=str(self.source.resolve()),
                 source_sha256=four.sha256_file(self.source), start=0, stop=1, sample=self.sample.name,
                 adapter_sha256=four.sha256_file(Path(ho.__file__))), sliced)
        self.failed = four.sample_decay_dir(self.sample, self.cfg)
        self.failed.mkdir(parents=True)
        (self.failed / "stable_input_source.lhe").write_text(raw)
        staged = self.failed / "stable_input.lhe"
        staged.write_text(raw.replace("5 0\n", "5 4.7\n").replace("BLOCK YUKAWA\n", "BLOCK YUKAWA\n5 4.7\n"))
        four.write_artifact_stage_marker(self.failed / "fourlepton_madspin_banner_stage.json",
            dict(runner_sha256=repair.OLD_RUNNER, source_sha256=four.sha256_file(self.failed / "stable_input_source.lhe")), staged)

    def invoke(self, **kwargs):
        with patch("sys.stdout", new=io.StringIO()):
            return repair.repair(self.cfg, **kwargs)

    def test_preview_is_read_only_then_apply_archives_and_preserves_outputs(self):
        self.invoke()
        self.assertEqual(self.cfg.manifest_path.read_bytes(), self.original_manifest)
        self.assertFalse((self.cfg.run_dir / "recovery").exists())
        digest = four.sha256_file(self.source)
        self.invoke(apply=True)
        backup = self.cfg.run_dir / "recovery" / repair.MIGRATION
        self.assertEqual((backup / self.cfg.manifest_path.name).read_bytes(), self.original_manifest)
        self.assertTrue((backup / "failed_madspin/stable_input.lhe").exists())
        self.assertFalse(self.failed.exists())
        self.assertEqual(four.sha256_file(self.source), digest)
        updated = json.loads(self.cfg.manifest_path.read_text())
        self.assertEqual(four.manifest_sample(updated, self.sample.name)["status"], "planned")
        for sample in self.samples[:6]:
            entry = four.manifest_sample(updated, sample.name)
            repair.verify_completed(entry, sample, self.cfg)
            self.assertEqual(entry["completion_software"]["runner_sha256"], repair.OLD_RUNNER)
        merged = four.merge_manifest_samples(updated, four.build_manifest(self.cfg, self.samples))
        self.assertEqual(merged["software_migrations"], updated["software_migrations"])
        self.invoke(apply=True)  # Idempotent: no second archive or replacement.

    def test_changed_completed_output_is_rejected(self):
        Path(self.manifest["samples"][0]["output_root"]).write_text("changed")
        with self.assertRaisesRegex(ValueError, "changed recovery artifact"):
            self.invoke(apply=True)
        self.assertEqual(self.cfg.manifest_path.read_bytes(), self.original_manifest)
        self.assertTrue(self.failed.exists())

    def test_changed_configuration_and_unreviewed_runner_are_rejected(self):
        for setting in ("seed_base", "higgs_mass"):
            with self.subTest(setting=setting):
                original = getattr(self.cfg, setting)
                setattr(self.cfg, setting, original + 1)
                with self.assertRaisesRegex(ValueError, "changed beyond"):
                    self.invoke()
                setattr(self.cfg, setting, original)
        with patch.object(repair, "FIXED_RUNNER", "0" * 64):
            with self.assertRaisesRegex(ValueError, "outside this audited migration"):
                self.invoke()

    def test_existing_decayed_output_and_changed_production_are_rejected(self):
        output = self.failed / "stable_input_decayed.lhe.gz"
        output.write_bytes(b"partial")
        with self.assertRaisesRegex(ValueError, "decayed WH output exists"):
            self.invoke()
        output.unlink()
        self.source.write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "hard production"):
            self.invoke()


if __name__ == "__main__":
    unittest.main()
