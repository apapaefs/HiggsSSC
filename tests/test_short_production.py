"""Count and recovery gates on synthetic LHEs; no integration or showering."""
from contextlib import ExitStack
from dataclasses import replace
import gzip
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from hfourlepton import checked_campaign as checked
from hfourlepton import retry_short_production as retry
from hfourlepton import minimal_runtime, higher_order
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton.repair_loop_production import completed_production


def lhe(events, settings=None):
    header = "\n".join(f" {value} = {key}" for key, value in (settings or {}).items())
    return ('<LesHouchesEvents version="3.0">\n<header>\n' + header + '\n</header>\n<init>\n'
            '2212 2212 6500 6500 0 0 335900 335900 3 1\n'
            '0.001 0.00001 0.001 1\n</init>\n' +
            '<event>\n1 1 0.001 125 0.007 0.118\n13 1 0 0 0 0 0 0 10 10 0 0 9\n</event>\n' * events +
            '</LesHouchesEvents>\n')


def card_text():
    return ('2000 = nevents\n0 = iseed\n6500 = ebeam1\n6500 = ebeam2\n'
            'lhapdf = pdlabel\n335900 = lhaid\n3 = ptl\n3 = etal\n4 = mmll\n'
            '0.02 = drll\n1 = sde_strategy ! strategy\n0 = hard_survey ! survey\n')


class ProductionChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for context in (patch.object(four, "SCRIPT_DIR", self.root), patch.dict(os.environ, {}, clear=True)):
            context.start()
            self.addCleanup(context.stop)
        self.cfg = replace(four.parse_config(["--order", "ho", "--preset", "minimal", "--run-tag", "count_test",
            "--run-samples", "bkg_gg_4e,bkg_gg_4mu", "--ebeam", "6500", "--nb-core", "40",
            "--no-herwig-env", "--no-herwig-module", "--dry-run"]), dry_run=False)
        self.samples = four.selected_samples(self.cfg)
        self.sample = next(s for s in self.samples if s.name == "bkg_gg_4mu")

    def production(self, sample, events):
        group = four.production_dir(sample, self.cfg)
        process = group / "mg5_process"
        path = process / "Events" / self.cfg.run_tag / "unweighted_events.lhe.gz"
        path.parent.mkdir(parents=True)
        with gzip.open(path, "wt") as stream:
            stream.write(lhe(events))
        (process / "Cards").mkdir()
        (process / "Cards/run_card.dat").write_text(card_text())
        (group / "cards").mkdir()
        (group / "cards/process_card.dat").write_text(four.process_card_text(sample, process, self.cfg))
        identity = four.production_stage_identity(sample, sample, self.samples, process, self.cfg)
        four.write_artifact_stage_marker(group / "fourlepton_production_stage.json", identity, path)
        return path

    def manifest(self):
        current = {"madspin_worker": {"sha256": retry.startup.FIXED_WORKER}}
        versions = dict(python="test", mg5="3.5.15", herwig="test", root="test")
        manifest = dict(tracked_inputs=current, software={k + "_version": v for k, v in versions.items()},
            configuration_fingerprint=four.configuration_fingerprint(self.cfg, self.samples,
                tracked_inputs=current, software_versions=versions),
            samples=[dict(name=s.name, status="planned") for s in self.samples])
        self.cfg.run_dir.mkdir(parents=True, exist_ok=True)
        four.write_manifest(manifest, self.cfg)
        return current

    def test_count_gzip_and_reject_truncation_and_short_output(self):
        path = self.production(self.sample, 19)
        self.assertEqual(checked.count_events(path), 19)
        with self.assertRaisesRegex(ValueError, "19 events, expected at least 2000"):
            checked.require_events(path, 2000)
        self.assertEqual(checked.require_events(path, 19), 19)
        for broken in (lhe(2).replace('</LesHouchesEvents>', ''),
                       lhe(2).replace('</event>', '', 1),
                       lhe(2).replace('<event>', '<event>\n<event>', 1)):
            with gzip.open(path, "wt") as stream:
                stream.write(broken)
            with self.assertRaisesRegex(ValueError, "truncated|malformed"):
                checked.count_events(path)

    def test_coordinator_stops_before_worker_creation_and_preserves_short_output(self):
        path = self.production(self.sample, 19)
        before = path.read_bytes()
        with patch.object(four, "obtain_production_lhe", return_value=path) as producer, \
             patch.object(minimal_runtime, "ThreadPoolExecutor", side_effect=AssertionError("workers started")):
            with checked.production_checks(), self.assertRaisesRegex(ValueError, "19 events"):
                minimal_runtime.execute_campaign([self.sample], {}, self.cfg)
            self.assertIs(four.obtain_production_lhe, producer)
        self.assertFalse((self.cfg.run_dir / "workers").exists())
        self.assertFalse((self.cfg.run_dir / "production_validation").exists())
        self.assertEqual(path.read_bytes(), before)

    def test_complete_group_has_count_receipt_without_changing_native_marker(self):
        path = self.production(self.sample, 2000)
        marker = four.production_dir(self.sample, self.cfg) / "fourlepton_production_stage.json"
        before = marker.read_bytes(), path.read_bytes()
        with patch.object(four, "obtain_production_lhe", return_value=path), checked.production_checks():
            four.obtain_production_lhe(self.sample, self.samples, self.cfg)
        receipt = json.loads((self.cfg.run_dir / "production_validation/bkg_gg_4mu.json").read_text())
        self.assertEqual(receipt["events"], 2000)
        self.assertEqual(receipt["required_events"], 2000)
        self.assertEqual(before, (marker.read_bytes(), path.read_bytes()))

    def test_signal_count_gate_checks_native_group_instead_of_channel_slice(self):
        samples = [s for s in higher_order.all_samples() if s.production_mode == "ggf"]
        native, sliced = self.root / "native.lhe", self.root / "slice.lhe"
        native.write_text(lhe(6))
        sliced.write_text(lhe(2))
        cfg = replace(self.cfg, signal_lhe=native, sample_nevents={s.name: 2 for s in samples})
        with patch.object(four, "obtain_production_lhe", return_value=sliced), \
             patch.object(higher_order, "inspect_production", return_value={"events": 6}) as inspect, \
             checked.production_checks():
            self.assertEqual(four.obtain_production_lhe(samples[0], samples, cfg), sliced)
        inspect.assert_called_once_with(native, samples[0], cfg)
        record = json.loads((cfg.run_dir / "production_validation/signal_gg_h_stable.json").read_text())
        self.assertEqual(record["events"], 6)
        self.assertEqual(record["required_events"], 6)

    def test_audit_validates_other_groups_and_rejects_altered_provenance(self):
        good = next(s for s in self.samples if s != self.sample)
        self.production(good, 2000)
        short = self.production(self.sample, 19)
        current = self.manifest()
        before = {p: p.read_bytes() for p in self.cfg.run_dir.rglob('*') if p.is_file()}
        with patch.object(four, "campaign_tracked_inputs", return_value=current):
            result = retry.audit({"state": "failed"}, self.cfg, self.sample.name)
            self.assertEqual(result["short_production"]["events"], 19)
            self.assertEqual(len(result["completed_production"]), 1)
            self.assertFalse(result["generators_launched"])
            self.assertEqual(before, {p: p.read_bytes() for p in before})
            with self.assertRaisesRegex(ValueError, "expected at least 2000"):
                completed_production(self.sample, self.samples, self.cfg)
            with gzip.open(short, "wt") as stream:
                stream.write(lhe(20))
            with self.assertRaisesRegex(ValueError, "provenance"):
                retry.audit({"state": "failed"}, self.cfg, self.sample.name)

    def test_sampling_override_changes_only_two_numerical_settings(self):
        card = self.root / "run_card.dat"
        card.write_text(card_text())
        retry.patch_sampling(card)
        self.assertEqual(retry.numerical_settings(card.read_text()), retry.SETTINGS)
        self.assertEqual(card.read_text(), card_text().replace('1 = sde_strategy', '2 = sde_strategy')
                         .replace('0 = hard_survey', '1 = hard_survey'))
        with self.assertRaisesRegex(ValueError, "unexpected original"):
            retry.patch_sampling(card)

    def test_actual_producer_rejects_short_retry_before_completion_marker(self):
        self.fake_generation(19)

    def test_actual_producer_records_complete_retry_and_checks_receipt(self):
        result = self.fake_generation(2000)
        self.assertEqual(result["events"], 2000)
        marker = four.production_dir(self.sample, self.cfg) / "fourlepton_production_stage.json"
        checked.validate_retry_receipt(marker)
        marker.with_name("numerical_retry.json").write_text('{}')
        with self.assertRaisesRegex(ValueError, "retry provenance"):
            checked.validate_retry_receipt(marker)

    def fake_generation(self, count):
        """Exercise the real LO producer, replacing only export/compile/generation."""
        group = four.production_dir(self.sample, self.cfg)
        process = group / "mg5_process"

        def export(*args, **kwargs):
            (process / "Cards").mkdir(parents=True)
            (process / "Cards/run_card.dat").write_text(card_text())

        def generate(*args, **kwargs):
            settings = retry.numerical_settings((process / "Cards/run_card.dat").read_text())
            self.assertEqual(settings, retry.SETTINGS)
            path = process / "Events" / self.cfg.run_tag / "unweighted_events.lhe.gz"
            path.parent.mkdir(parents=True)
            with gzip.open(path, "wt") as stream:
                stream.write(lhe(count, settings))

        with ExitStack() as stack:
            stack.enter_context(patch.object(four, "run_runtime_command", side_effect=export))
            stack.enter_context(patch.object(four, "run_generate_events", side_effect=generate))
            stack.enter_context(patch.object(minimal_runtime, "configure_export"))
            stack.enter_context(patch.object(higher_order, "inspect_production", return_value={"cross_section_pb": .001}))
            for function in ("patch_mg5_toolchain", "patch_mg5_hidden_pdf_defaults", "patch_mg5_pdf_label_sync",
                             "patch_param_card_higgs_mass", "generate_and_patch_run_card_include", "restore_mg5_symmetry_factors"):
                stack.enter_context(patch.object(four, function, return_value=None))
            if count < 2000:
                with self.assertRaisesRegex(ValueError, "19 events"):
                    retry.execute_retry(self.sample, self.samples, self.cfg, {"id": retry.RETRY})
                self.assertFalse((group / "fourlepton_production_stage.json").exists())
                self.assertTrue((process / "Events" / self.cfg.run_tag / "unweighted_events.lhe.gz").exists())
                return
            return retry.execute_retry(self.sample, self.samples, self.cfg, {"id": retry.RETRY})


if __name__ == "__main__":
    unittest.main()
