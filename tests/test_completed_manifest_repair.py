"""Strict recovery of timing-only resume metadata, without simulation programs."""
from copy import deepcopy
from dataclasses import asdict, replace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_pipeline as cli
from hfourlepton import repair_completed_manifest as repair
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton.minimal_runtime import _encode
from test_fourlepton_minimal import config
import test_student_pipeline as pipeline_fixtures


class CompletedManifestRepairTests(unittest.TestCase):
    def setUp(self):
        helper = pipeline_fixtures.StudentPipelineTests()
        self.context = helper.repository()
        self.root = self.context.__enter__()
        self.addCleanup(self.context.__exit__, None, None, None)
        self.args = helper.args(self.root, "--workflow", "fourlepton-ho", "--energy-tev", "13")
        helper.fourlepton_artifacts(self.args)
        self.directory = cli.locations(self.args)[0]
        self.manifest_path = self.directory / "manifest_perfect.json"
        self.original = cli.physics_artifacts(self.args)
        self.saved = {"state": "failed", "stage": "fourlepton-perfect", "analysis_artifacts": self.original}
        self.state = cli.state_path(self.args)
        cli.shared.save(self.state, self.saved)
        self.manifest = json.loads(self.manifest_path.read_text())
        self.manifest["production_wall_seconds"]["sample"] = 0.01
        cli.shared.save(self.manifest_path, self.manifest)
        self.cfg = SimpleNamespace(run_dir=self.directory, manifest_path=self.manifest_path)

    def report(self):
        products, old, new = repair.frozen_inventory(self.saved, self.manifest_path, self.manifest)
        return {"id": repair.REPAIR, "at": "2026-11-02T23:00:00+00:00", "unchanged_products": products,
                "old_manifest_sha256": old, "new_manifest_sha256": new,
                "completed_production": {}, "samples": [s["name"] for s in self.manifest["samples"]]}, []

    def test_only_manifest_difference_is_allowed_and_all_outputs_remain_pinned(self):
        products, old, new = repair.frozen_inventory(self.saved, self.manifest_path, self.manifest)
        self.assertEqual(len(products), 54)
        self.assertNotEqual(old, new)
        sample = self.manifest["samples"][0]
        for path_key, hash_key in (("root_file", "output_root_sha256"),
                                  ("analysis_summary", "analysis_summary_sha256")):
            product = Path(sample[path_key])
            before = product.read_bytes()
            product.write_bytes(b"changed completed output")
            with self.assertRaisesRegex(ValueError, "ROOT/summary changed"):
                repair.frozen_inventory(self.saved, self.manifest_path, self.manifest)
            # Updating the unfrozen manifest's hash must not defeat the snapshot.
            original_hash = sample[hash_key]
            sample[hash_key] = cli.shared.digest(product)
            with self.assertRaisesRegex(ValueError, "ROOT/summary changed"):
                repair.frozen_inventory(self.saved, self.manifest_path, self.manifest)
            sample[hash_key] = original_hash
            product.write_bytes(before)

    def test_changed_inventory_is_not_adopted(self):
        self.saved["analysis_artifacts"][str(self.root / "unrelated.root")] = "a" * 64
        with self.assertRaisesRegex(ValueError, "inventory differs"):
            repair.frozen_inventory(self.saved, self.manifest_path, self.manifest)

    def test_audit_apply_and_repeat_preserve_all_products_and_archive_old_state(self):
        state_before = self.state.read_bytes()
        before = {p: (Path(p).read_bytes(), Path(p).stat().st_mtime_ns) for p in self.original}
        report = self.report()
        with patch.object(repair, "audit", return_value=report):
            repair.repair(self.saved, self.state, self.cfg, self.args)
            self.assertEqual(self.state.read_bytes(), state_before)
            repair.repair(self.saved, self.state, self.cfg, self.args, apply=True)
        self.assertEqual(before, {p: (Path(p).read_bytes(), Path(p).stat().st_mtime_ns) for p in before})
        backup = self.directory / "recovery" / repair.REPAIR
        self.assertEqual((backup / "managed_campaign.json").read_bytes(), state_before)
        self.assertEqual((backup / "manifest_perfect.json").read_bytes(), self.manifest_path.read_bytes())
        updated = cli.read_state(self.args)
        changed = [p for p in self.original if updated["analysis_artifacts"][p] != self.original[p]]
        self.assertEqual(changed, [str(self.manifest_path)])
        self.assertEqual(updated["state"], "prepared")
        self.assertEqual(updated["stage"], "report")
        self.assertIsNone(updated["error"])
        state_after = self.state.read_bytes()
        repair.repair(updated, self.state, self.cfg, self.args, apply=True)
        self.assertEqual(self.state.read_bytes(), state_after)

    def test_verification_failure_rolls_back_managed_state(self):
        before = self.state.read_bytes()
        with patch.object(repair, "audit", return_value=self.report()), \
                patch.object(cli, "physics_artifacts", side_effect=ValueError("changed during repair")):
            with self.assertRaisesRegex(ValueError, "changed during repair"):
                repair.repair(self.saved, self.state, self.cfg, self.args, apply=True)
        self.assertEqual(self.state.read_bytes(), before)
        self.assertEqual(cli.shared.digest(self.manifest_path), self.report()[0]["new_manifest_sha256"])

    def test_concurrent_metadata_change_is_rejected_before_archive_or_write(self):
        result = self.report()
        def audit(*unused):
            self.state.write_text('{"changed":true}')
            return result
        with patch.object(repair, "audit", side_effect=audit):
            with self.assertRaisesRegex(ValueError, "metadata changed"):
                repair.repair(self.saved, self.state, self.cfg, self.args, apply=True)
        self.assertEqual(self.state.read_text(), '{"changed":true}')
        self.assertFalse((self.directory / "recovery" / repair.REPAIR).exists())


class CompletedWorkerEvidenceTests(unittest.TestCase):
    def test_worker_evidence_allows_only_completed_resume_and_timing_status_changes(self):
        with TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)):
            cfg = replace(config("--run-samples", "signal_gg_h_4e"), dry_run=False)
            samples = four.selected_samples(cfg)
            name = samples[0].name
            directory = cfg.run_dir / "workers/perfect"
            directory.mkdir(parents=True)
            journal = directory / (name + ".json")
            job_path = directory / (name + ".job.json")
            log = directory / (name + ".log")
            current = {"status": "complete", "production_wall_seconds": {"group": .1},
                       "samples": [{"name": name, "status": "complete", "last_attempt_wall_seconds": .02,
                                    "normalization": {"weight_scale": 1.}}]}
            old = deepcopy(current)
            old.update(status="running", production_wall_seconds={"group": 100.})
            old["samples"][0]["last_attempt_wall_seconds"] = 120.
            worker_cfg = replace(cfg, postprocess_only=True, manifest_override=journal)
            job = {"config": asdict(worker_cfg), "sample": name, "manifest": old}
            job_path.write_text(json.dumps(job, default=_encode))
            journal.write_text(json.dumps(old))
            log.write_text(f"==> Resuming: {name} is already complete\n")
            self.assertEqual(len(repair.worker_evidence(cfg, samples, current)), 3)
            old["samples"][0]["normalization"]["weight_scale"] = 2.
            journal.write_text(json.dumps(old))
            with self.assertRaisesRegex(ValueError, "beyond timing/status"):
                repair.worker_evidence(cfg, samples, current)
            old["samples"][0]["normalization"]["weight_scale"] = 1.
            old["samples"][0]["status"] = "planned"
            journal.write_text(json.dumps(old))
            with self.assertRaisesRegex(ValueError, "beyond timing/status"):
                repair.worker_evidence(cfg, samples, current)
            old["samples"][0]["status"] = "complete"
            journal.write_text(json.dumps(old))
            log.write_text("Running a new shower\n")
            with self.assertRaisesRegex(ValueError, "did not only resume"):
                repair.worker_evidence(cfg, samples, current)


if __name__ == "__main__":
    unittest.main()
