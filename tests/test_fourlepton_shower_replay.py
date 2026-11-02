"""Finite-LHE boundary and reversible recovery checks; no simulation programs."""
from copy import deepcopy
from dataclasses import replace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hfourlepton import higher_order as ho, repair_shower_replay as repair
from hfourlepton import run_four_lepton_campaign as four
from test_fourlepton_minimal import config
import test_ho_shower_completion as fixtures


class ShowerReplayTests(unittest.TestCase):
    def test_real_launcher_requests_two_slots_beyond_source_and_records_eof(self):
        sample = ho.all_samples()[0]
        with TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)), \
                patch.object(fixtures, "SAMPLE", sample.name):
            cfg = replace(config(), dry_run=False, nevents=3, sample_nevents={sample.name: 3})
            work = four.shower_work_dir(sample, cfg)
            fixture = fixtures.ShowerCompletionTests().fixture(work, termination="source_exhausted")
            directory, source, metadata = fixture[:3]
            root = Path(metadata[0]["path"])
            root.unlink()
            (directory / (sample.name + ".in")).unlink()
            calls = []

            def command(argv, _cfg, _cwd, log):
                calls.append(list(map(str, argv)))
                if argv[1] == "run":
                    self.assertEqual(argv[-1], "-N5")
                    root.write_bytes(b"synthetic finalized ROOT")
                    metadata[0].update(size=root.stat().st_size, mtime_ns=root.stat().st_mtime_ns)
                    return 1  # The deliberate, fully checked first EOF.
                if log.name == "root-metadata.log":
                    log.write_text(json.dumps(metadata) + "\n")
                return 0

            with patch.object(ho, "logged_command", side_effect=command):
                ho.run_shower(sample, source, cfg)
            record = json.loads((directory / "matched_completion.json").read_text())
            self.assertEqual(record["execution_plan"]["execution_events"], 5)
            self.assertEqual(record["saved_events"], 3)
            self.assertEqual(record["attempted_events"], 3)
            self.assertEqual(record["termination"], "source_exhausted")
            self.assertEqual(len(calls), 3)

    def test_replay_inspection_requires_known_counts_card_and_unmodified_roots(self):
        sample = ho.all_samples()[0]
        with TemporaryDirectory() as tmp, patch.object(four, "SCRIPT_DIR", Path(tmp)), \
                patch.object(fixtures, "SAMPLE", sample.name):
            cfg = replace(config(), dry_run=False, nevents=3, sample_nevents={sample.name: 3})
            fixture = fixtures.ShowerCompletionTests().fixture(four.shower_work_dir(sample, cfg),
                saved=3, attempts=4, requested=3)
            directory, source, metadata = fixture[:3]
            (directory / (sample.name + ".in")).write_text(ho.herwig_card(
                sample, source.resolve(), four.seed_for(sample.name + ":herwig", cfg), 3, cfg))
            (directory / "root-metadata.log").write_text(json.dumps(metadata) + "\n")
            result = repair.inspect_replay(sample, cfg)
            self.assertEqual(result["counts"]["attempted_events"], 4)
            Path(metadata[0]["path"]).write_bytes(b"changed output")
            with self.assertRaisesRegex(ValueError, "ROOT changed"):
                repair.inspect_replay(sample, cfg)


class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        patcher = patch.object(four, "SCRIPT_DIR", self.root)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.cfg = replace(config("--run-samples", "signal_gg_h_4e"), dry_run=False)
        self.samples = four.selected_samples(self.cfg)
        self.name = self.samples[0].name
        current = {"higher_order_adapter": {"sha256": repair.FIXED_ADAPTER}, "runner": {"sha256": "same"}}
        self.current = current
        patcher = patch.object(four, "campaign_tracked_inputs", return_value=current)
        patcher.start()
        self.addCleanup(patcher.stop)
        old = deepcopy(current)
        old["higher_order_adapter"]["sha256"] = repair.OLD_ADAPTER
        versions = dict(python="3.9", mg5="3.5.15", herwig="7.3.0", root="6.40.04")
        self.manifest = {"tracked_inputs": old,
            "software": {name + "_version": value for name, value in versions.items()},
            "configuration_fingerprint": four.configuration_fingerprint(self.cfg, self.samples,
                tracked_inputs=old, software_versions=versions),
            "samples": [{"name": self.name, "status": "failed", "error": "count mismatch"}]}
        self.cfg.run_dir.mkdir(parents=True)
        four.write_manifest(self.manifest, self.cfg)
        self.workers = self.cfg.run_dir / "workers/perfect"
        self.workers.mkdir(parents=True)
        (self.workers / "worker.log").write_text("old worker evidence")
        self.shower = four.shower_work_dir(self.samples[0], self.cfg) / "herwig"
        self.shower.mkdir(parents=True)
        (self.shower / "events.root").write_bytes(b"invalid shower preserved for audit")
        self.decay = four.sample_decay_dir(self.samples[0], self.cfg)
        self.decay.mkdir(parents=True)
        (self.decay / "decays.lhe").write_text("completed decays")
        self.marker = self.cfg.run_dir / "production/producer.json"
        self.marker.parent.mkdir()
        self.marker.write_text(json.dumps({"identity": {"adapter_sha256": repair.OLD_ADAPTER}}))
        self.state = self.root / "managed.json"
        self.state.write_text('{"state":"failed"}')
        for patcher in (
            patch.object(repair, "completed_production", return_value={"verified": True}),
            patch.object(repair, "stage_updates", return_value={self.marker: {
                "identity": {"adapter_sha256": repair.FIXED_ADAPTER}}}),
            patch.object(repair, "check_samples", return_value=([{"sample": self.name}], [])),
            patch.object(ho, "logged_command", side_effect=AssertionError("must not run simulations")),
        ):
            patcher.start()
            self.addCleanup(patcher.stop)

    def test_audit_apply_and_second_apply_preserve_outputs_and_archive_originals(self):
        original = self.cfg.manifest_path.read_bytes()
        marker = self.marker.read_bytes()
        report = repair.repair({"state": "failed"}, self.state, self.cfg)
        self.assertFalse(report["generators_launched"])
        self.assertEqual(self.cfg.manifest_path.read_bytes(), original)
        self.assertEqual(self.marker.read_bytes(), marker)
        repair.repair({"state": "failed"}, self.state, self.cfg, apply=True)
        backup = self.cfg.run_dir / "recovery" / repair.MIGRATION
        self.assertEqual((backup / "metadata/manifest_perfect.json").read_bytes(), original)
        self.assertEqual((backup / "metadata/production/producer.json").read_bytes(), marker)
        self.assertTrue((backup / "workers/perfect/worker.log").exists())
        self.assertTrue((backup / "showers" / self.name / "events.root").exists())
        self.assertFalse(self.shower.exists())
        self.assertEqual((self.decay / "decays.lhe").read_text(), "completed decays")
        updated = json.loads(self.cfg.manifest_path.read_text())
        self.assertEqual(updated["samples"], [{"name": self.name, "status": "planned"}])
        before = self.cfg.manifest_path.read_bytes()
        repair.repair({"state": "failed"}, self.state, self.cfg, apply=True)
        self.assertEqual(self.cfg.manifest_path.read_bytes(), before)

    def test_other_configuration_changes_are_not_migrated(self):
        self.current["runner"]["sha256"] = "unexpected"
        with self.assertRaisesRegex(ValueError, "beyond the shower launcher"):
            repair.repair({"state": "failed"}, self.state, self.cfg, apply=True)
        self.assertTrue(self.shower.exists())

    def test_failed_metadata_write_rolls_back_markers_shower_and_workers(self):
        original = self.cfg.manifest_path.read_bytes()
        marker = self.marker.read_bytes()
        with patch.object(four, "write_manifest", side_effect=OSError("write failed")):
            with self.assertRaisesRegex(OSError, "write failed"):
                repair.repair({"state": "failed"}, self.state, self.cfg, apply=True)
        self.assertEqual(self.cfg.manifest_path.read_bytes(), original)
        self.assertEqual(self.marker.read_bytes(), marker)
        self.assertTrue(self.workers.exists())
        self.assertTrue(self.shower.exists())


if __name__ == "__main__":
    unittest.main()
