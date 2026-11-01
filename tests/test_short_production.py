"""Count and recovery gates on synthetic LHEs; no integration or showering."""
from contextlib import ExitStack
from dataclasses import replace
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
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


def card_text(*, expanded=False):
    return ('2000 = nevents\n0 = iseed\n6500 = ebeam1\n6500 = ebeam2\n'
            'lhapdf = pdlabel\n335900 = lhaid\n3 = ptl\n3 = etal\n4 = mmll\n'
            '0.02 = drll\n1 = sde_strategy ! strategy\n' +
            ('0 = hard_survey ! survey\n' if expanded else ''))


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
        (process / "Cards/run_card.dat").write_text(card_text(expanded=True))
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
        for expanded in (False, True):
            card.write_text(card_text(expanded=expanded))
            retry.patch_sampling(card)
            self.assertEqual(retry.numerical_settings(card.read_text()), retry.SETTINGS)
            expected = card_text(expanded=expanded).replace('1 = sde_strategy', '2 = sde_strategy')
            expected = (expected.replace('0 = hard_survey', '1 = hard_survey') if expanded else
                        expected + ' 1 = hard_survey ! explicit numerical retry setting\n')
            self.assertEqual(card.read_text(), expected)
            with self.assertRaisesRegex(ValueError, "unexpected original"):
                retry.patch_sampling(card)

    def test_missing_hidden_setting_is_allowed_only_when_patching_fresh_cards(self):
        self.assertEqual(retry.numerical_settings(card_text(), fresh_card=True),
                         {"sde_strategy": 1, "hard_survey": 0})
        with self.assertRaisesRegex(ValueError, "explicit hard_survey"):
            retry.numerical_settings(card_text())
        card = self.root / "run_card.dat"
        for text in (card_text() + 'bad = hard_survey\n',
                     card_text() + '0 = hard_survey\n1 = hard_survey\n',
                     card_text() + '0 = hard_survey\nbad = hard_survey\n',
                     card_text().replace('1 = sde_strategy', 'bad = sde_strategy')):
            card.write_text(text)
            with self.assertRaisesRegex(ValueError, "exactly one explicit"):
                retry.patch_sampling(card)
            self.assertEqual(card.read_text(), text)

    def test_actual_mg5_fresh_card_and_parser_accept_the_explicit_override(self):
        repo = Path(__file__).resolve().parents[1]
        card = self.root / "fresh.dat"
        script = f'''import sys
from pathlib import Path
sys.path[:0] = [{str(repo / "MG5_aMC_v3_5_15")!r}, {str(repo)!r}]
from madgraph.various.banner import RunCardLO
from hfourlepton.retry_short_production import patch_sampling, numerical_settings
card = Path({str(card)!r})
original = RunCardLO()
assert original['hard_survey'] == 0
assert 'hard_survey' in original.hidden_param
original['sde_strategy'] = 1
original.write(str(card))
assert not any('= hard_survey' in line.replace('  ', ' ') for line in card.read_text().splitlines())
patch_sampling(card)
assert numerical_settings(card.read_text()) == {{'sde_strategy': 2, 'hard_survey': 1}}
parsed = RunCardLO(str(card))
assert parsed['sde_strategy'] == 2 and parsed['hard_survey'] == 1
for key, value in original.items():
    if key.lower() not in ('sde_strategy', 'hard_survey'):
        assert parsed[key] == value, key
print('MG5_FRESH_CARD_OK')
'''
        result = subprocess.run([sys.executable, '-c', script], cwd=self.root,
            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}, capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn('MG5_FRESH_CARD_OK', result.stdout)

    def startup_fixture(self):
        good = next(s for s in self.samples if s != self.sample)
        self.production(good, 2000)
        self.production(self.sample, 19)
        current = self.manifest()
        group = four.production_dir(self.sample, self.cfg)
        destination = self.cfg.run_dir / 'recovery' / retry.RETRY / self.sample.name
        saved = dict(state='failed', identity=dict(environment=dict(work_dir=str(self.root / 'work'))))
        with patch.object(four, 'campaign_tracked_inputs', return_value=current):
            report = retry.audit(saved, self.cfg, self.sample.name)
        report['retry_code']['sha256'] = retry.OLD_RETRY_SHA256
        retry.archive(group, destination, report, [self.cfg.manifest_path])
        process = group / 'mg5_process'
        (process / 'Events').mkdir(parents=True)
        (process / 'Events/.keep').touch()
        (process / 'Cards').mkdir()
        (process / 'Cards/run_card.dat').write_text(card_text())
        (group / 'cards').mkdir()
        (group / 'cards/process_card.dat').write_text(four.process_card_text(self.sample, process, self.cfg))
        log = self.root / 'work/logs' / (self.cfg.run_tag + '-retry-gg4mu.log')
        log.parent.mkdir(parents=True)
        log.write_text('ERROR: expected exactly one explicit hard_survey setting\n')
        return current, saved, destination, report

    def test_card_startup_audit_preserves_both_attempts_and_requires_no_integration(self):
        current, saved, destination, report = self.startup_fixture()
        group = four.production_dir(self.sample, self.cfg)
        before = {p: p.read_bytes() for p in self.cfg.run_dir.rglob('*') if p.is_file()}
        with patch.object(four, 'campaign_tracked_inputs', return_value=current):
            checked_report = retry.audit_card_startup(saved, self.cfg, self.sample, destination, report)
            self.assertIn('card_startup_recovery', checked_report)
            self.assertFalse(checked_report['generators_launched'])
            self.assertEqual(before, {p: p.read_bytes() for p in before})
            for name in ('numerical_retry.json', 'mg5_process/Source/run_card.inc',
                         'mg5_process/Events/count_test/events.lhe'):
                progress = group / name
                progress.parent.mkdir(parents=True, exist_ok=True)
                progress.write_text('existing output')
                with self.assertRaisesRegex(ValueError, 'progressed beyond'):
                    retry.audit_card_startup(saved, self.cfg, self.sample, destination, report)
                progress.unlink()
                if name.endswith('events.lhe'):
                    progress.parent.rmdir()
            retry.archive(group, destination / retry.CARD_RECOVERY, checked_report, [self.cfg.manifest_path])
        self.assertFalse(group.exists())
        self.assertTrue((destination / 'failed_production/mg5_process/Events/count_test/unweighted_events.lhe.gz').is_file())
        self.assertTrue((destination / retry.CARD_RECOVERY / 'failed_production/mg5_process/Cards/run_card.dat').is_file())

    def test_card_startup_recovery_rejects_changed_archive_or_other_production(self):
        current, saved, destination, report = self.startup_fixture()
        source = destination / 'failed_production/mg5_process/Events/count_test/unweighted_events.lhe.gz'
        original = source.read_bytes()
        with patch.object(four, 'campaign_tracked_inputs', return_value=current):
            source.write_bytes(b'changed archive')
            with self.assertRaisesRegex(ValueError, 'short-sample archive changed'):
                retry.audit_card_startup(saved, self.cfg, self.sample, destination, report)
            source.write_bytes(original)
            good = next(s for s in self.samples if s != self.sample)
            marker = four.production_dir(good, self.cfg) / 'fourlepton_production_stage.json'
            marker.write_text('{}')
            with self.assertRaisesRegex(ValueError, 'provenance'):
                retry.audit_card_startup(saved, self.cfg, self.sample, destination, report)

    def test_cli_audits_without_running_and_archives_only_when_operator_requests_run(self):
        current, saved, destination, _report = self.startup_fixture()
        state_path = self.root / 'managed.json'
        state_path.write_text(json.dumps(saved))
        work = self.root / 'work'
        command = ['--work-dir', str(work), '--run-tag', self.cfg.run_tag, '--sample', self.sample.name]
        before = {p: p.read_bytes() for p in self.cfg.run_dir.rglob('*') if p.is_file()}

        def simulate_after_archive(*args):
            self.assertTrue((destination / retry.CARD_RECOVERY / 'failed_production').is_dir())
            self.assertFalse(four.production_dir(self.sample, self.cfg).exists())
            self.assertTrue((destination / 'failed_production').is_dir())
            return {'events': 2000}

        with patch.object(four, 'campaign_tracked_inputs', return_value=current), \
             patch.object(retry, 'saved_configuration', return_value=(saved, state_path, self.cfg)), \
             patch.object(retry, 'saved_options', return_value=(saved, state_path, None, {})), \
             patch.object(retry.cli, 'StudentPipeline', return_value=SimpleNamespace(env=dict(os.environ), mg5=self.root)), \
             patch.object(minimal_runtime, 'validate_libraries'), \
             patch('hfourlepton.mg5_loop_recovery.prepare_campaign'), \
             patch.object(retry, 'execute_retry', side_effect=simulate_after_archive) as generate:
            retry.main(command)
            generate.assert_not_called()
            self.assertEqual(before, {p: p.read_bytes() for p in before})
            retry.main([*command, '--run'])
            generate.assert_called_once()

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
