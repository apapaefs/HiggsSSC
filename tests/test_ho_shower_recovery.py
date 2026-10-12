"""Preparing a failed-shower restart preserves completed generation and evidence."""

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hgammagamma import run_gammagamma_ho_campaign as ho


def fingerprint(configuration):
    return hashlib.sha256(json.dumps(configuration, sort_keys=True).encode()).hexdigest()


class FailedShowerRecoveryTests(unittest.TestCase):
    sample = ho.SAMPLES[1]

    def arguments(self, directory, *extra):
        return ho.parse_args(["--output-dir", str(directory), "--no-herwig-module",
                              "--nevents", "1000", *extra])

    def fixture(self, directory):
        original_args = self.arguments(directory)
        old = ho.prepare(original_args, self.sample)
        sample_dir = ho.sample_dir(original_args, self.sample)
        source = ho.lhe_path(original_args, self.sample)
        source.parent.mkdir(parents=True)
        source.write_bytes(b"unchanged complete generated LHE fixture")
        old["completed"] = ["generate"]
        old["lhe"] = {"events": 1000, "source": str(source), "cross_section_pb": 1.0}
        manifest = sample_dir / "campaign.json"
        ho.write_json(manifest, old, original_args)
        herwig = sample_dir / "herwig"
        (herwig / "events").mkdir()
        root = herwig / "events/partial.root"
        root.write_bytes(b"finalized partial ROOT fixture")
        for name in ("run.log", f"{self.sample.name}.log", f"{self.sample.name}.out"):
            (herwig / name).write_text("original error-limit diagnostic\n")
        counts = dict.fromkeys(ho.shower_completion.SEVERITIES, 0)
        counts["eventerror"] = 100
        parsed = {
            "generated_events": 378, "saved_events": 378, "attempted_events": 478,
            "exception_counts": counts,
            "exception_classes": [{"class": "ThePEG::Exception", "severity": "eventerror", "count": 100}],
            "root_writers": [str(root.resolve())],
            "artifacts": {"run_log": ho.shower_completion._artifact(herwig / "run.log")},
        }
        args = self.arguments(directory, "--resume", "--restart-failed-shower",
                              "--herwig-max-errors", "1000")
        configuration = copy.deepcopy(old["configuration"])
        configuration["herwig_card"] = ho.herwig_card(args, self.sample)
        return args, old, configuration, parsed, sample_dir, source, manifest, root

    def invoke(self, fixture, **options):
        args, old, configuration, parsed, *_ = fixture
        with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])) as inspect, \
                patch.object(ho.shower_completion, "_parse_logs", return_value=copy.deepcopy(parsed)) as logs, \
                patch.object(ho, "command") as command, contextlib.redirect_stdout(io.StringIO()):
            updated = ho.recover_failed_shower(args, self.sample, old, configuration,
                                              fingerprint(configuration))
            command.assert_not_called()
            return updated, inspect, logs

    def test_default_budget_preserves_existing_card_and_fingerprint(self):
        with TemporaryDirectory() as directory:
            args = self.arguments(directory)
            self.assertEqual(args.herwig_max_errors, 100)
            card = ho.herwig_card(args, self.sample)
            self.assertEqual(card.count("set EventGenerator:MaxErrors 100\n"), 1)
            old = ho.prepare(args, self.sample)
            explicit = self.arguments(directory, "--herwig-max-errors", "100")
            self.assertEqual(ho.herwig_card(explicit, self.sample), card)
            self.assertEqual(ho.prepare(explicit, self.sample), old)
            custom = self.arguments(directory, "--herwig-max-errors", "1000")
            self.assertEqual(ho.herwig_card(custom, self.sample).replace(
                "MaxErrors 1000\n", "MaxErrors 100\n"), card)

    def test_restart_cli_requires_prepare_resume_and_positive_budget(self):
        cases = (("--herwig-max-errors", "0"), ("--herwig-max-errors", "-1"),
                 ("--herwig-max-errors", "1.5"), ("--restart-failed-shower",),
                 ("--restart-failed-shower", "--resume", "--stage", "shower"),
                 ("--restart-failed-shower", "--resume", "--stage", "all"))
        for extra in cases:
            with self.subTest(extra=extra), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    self.arguments("/tmp/recovery-cli", *extra)

    def test_recovery_archives_failed_shower_and_preserves_generation_bytes(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            args, old, configuration, parsed, sample_dir, source, manifest, root = fixture
            manifest_bytes = manifest.read_bytes()
            source_stat = source.stat()
            source_bytes = source.read_bytes()
            root_bytes = root.read_bytes()
            updated, inspect, logs = self.invoke(fixture)
            archive = Path(updated["shower_restarts"][-1]["archive"])
            self.assertEqual((archive / "campaign.json").read_bytes(), manifest_bytes)
            self.assertEqual((archive / "herwig/events/partial.root").read_bytes(), root_bytes)
            self.assertEqual((archive / "herwig" / f"{self.sample.name}.in").read_text(),
                             old["configuration"]["herwig_card"])
            self.assertEqual((sample_dir / "herwig" / f"{self.sample.name}.in").read_text(),
                             configuration["herwig_card"])
            self.assertFalse((sample_dir / "herwig/events").exists())
            self.assertEqual(source.read_bytes(), source_bytes)
            self.assertEqual(source.stat().st_mtime_ns, source_stat.st_mtime_ns)
            self.assertEqual(updated["lhe"], old["lhe"])
            self.assertEqual(updated["completed"], ["generate"])
            self.assertEqual(json.loads(manifest.read_text()), updated)
            self.assertEqual(updated["fingerprint"], fingerprint(configuration))
            self.assertNotIn("shower_restarts", old)
            inspect.assert_called_once_with(source, self.sample, args)
            logs.assert_called_once_with(sample_dir / "herwig", self.sample.name, "error_limit")

    def test_identical_recovery_request_is_a_noop(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            updated, _, _ = self.invoke(fixture)
            repeat = (fixture[0], updated, fixture[2], *fixture[3:])
            result, inspect, logs = self.invoke(repeat)
            self.assertIs(result, updated)
            inspect.assert_not_called()
            logs.assert_not_called()
            self.assertEqual(len(updated["shower_restarts"]), 1)

    def test_prepare_recovery_publishes_identity_then_ordinary_resume_succeeds(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            args, old, configuration, parsed, sample_dir, *_ = fixture
            with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                    patch.object(ho.shower_completion, "_parse_logs", return_value=copy.deepcopy(parsed)), \
                    contextlib.redirect_stdout(io.StringIO()):
                updated = ho.prepare(args, self.sample)
            self.assertEqual(updated["configuration"], configuration)
            self.assertEqual(updated["fingerprint"], fingerprint(configuration))
            self.assertEqual(len(updated["shower_restarts"]), 1)
            self.assertEqual(json.loads((sample_dir / "campaign.json").read_text()), updated)
            args.restart_failed_shower = False
            with patch.object(ho, "inspect_lhe") as inspect, \
                    patch.object(ho.shower_completion, "_parse_logs") as logs:
                self.assertEqual(ho.prepare(args, self.sample), updated)
            inspect.assert_not_called()
            logs.assert_not_called()

    def test_prepare_recovery_dry_run_leaves_original_root_and_manifest(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            args, old, configuration, parsed, sample_dir, *_ = fixture
            args.dry_run = True
            before = {str(path): path.read_bytes() for path in Path(directory).rglob("*") if path.is_file()}
            with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                    patch.object(ho.shower_completion, "_parse_logs", return_value=copy.deepcopy(parsed)), \
                    contextlib.redirect_stdout(io.StringIO()):
                updated = ho.prepare(args, self.sample)
            after = {str(path): path.read_bytes() for path in Path(directory).rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            self.assertEqual(updated["configuration"], configuration)
            self.assertEqual(updated["fingerprint"], fingerprint(configuration))
            self.assertFalse((sample_dir / "shower-recovery").exists())

    def test_dry_run_checks_inputs_but_does_not_archive_or_write(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            fixture[0].dry_run = True
            before = {str(path): path.read_bytes() for path in Path(directory).rglob("*") if path.is_file()}
            updated, inspect, logs = self.invoke(fixture)
            after = {str(path): path.read_bytes() for path in Path(directory).rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            self.assertFalse((fixture[4] / "shower-recovery").exists())
            self.assertEqual(updated["configuration"], fixture[2])
            inspect.assert_called_once()
            logs.assert_called_once()

    def test_changed_generation_or_nonincreased_budget_is_rejected(self):
        for change in ("seed", "runtime", "budget", "duplicate_budget"):
            with self.subTest(change=change), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                configuration = fixture[2]
                if change == "seed":
                    configuration["launch_card"] += "set iseed 9999\n"
                elif change == "runtime":
                    configuration["herwig"] = "OtherHerwig"
                elif change == "budget":
                    configuration["herwig_card"] = configuration["herwig_card"].replace("MaxErrors 1000", "MaxErrors 10")
                else:
                    configuration["herwig_card"] += "set EventGenerator:MaxErrors 1000\n"
                with self.assertRaisesRegex(ValueError, "MaxErrors"):
                    self.invoke(fixture)
                self.assertTrue(fixture[7].exists())
                self.assertFalse((fixture[4] / "shower-recovery").exists())

    def test_completed_or_analyzed_campaign_is_rejected(self):
        for key in ("shower", "analyze", "shower_completion", "analysis", "rate_normalization", "root_files"):
            with self.subTest(key=key), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                old = fixture[1]
                if key in ("shower", "analyze"):
                    old["completed"].append(key)
                else:
                    old[key] = {}
                with self.assertRaisesRegex(ValueError, "incomplete, unanalyzed"):
                    self.invoke(fixture)
                self.assertTrue(fixture[7].exists())

    def test_failed_evidence_must_match_old_limit_and_partial_population(self):
        for change in ("under_limit", "over_limit", "fatal", "runerror", "source_exhausted", "requested_complete", "writers"):
            with self.subTest(change=change), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                parsed = fixture[3]
                if change in ("under_limit", "over_limit"):
                    count = 99 if change == "under_limit" else 101
                    parsed["exception_classes"][0]["count"] = count
                    parsed["exception_counts"]["eventerror"] = count
                elif change in ("fatal", "runerror"):
                    parsed["exception_counts"]["abortnow" if change == "fatal" else "runerror"] = 1
                elif change == "source_exhausted":
                    parsed["attempted_events"] = 1000
                elif change == "requested_complete":
                    parsed["generated_events"] = parsed["saved_events"] = 1000
                else:
                    parsed["root_writers"] = []
                with self.assertRaises(ValueError):
                    self.invoke(fixture)
                self.assertTrue(fixture[7].exists())
                self.assertFalse((fixture[4] / "shower-recovery").exists())

    def test_changed_card_or_manifest_fingerprint_is_rejected(self):
        for change in ("card", "fingerprint", "beam"):
            with self.subTest(change=change), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                if change == "card":
                    (fixture[4] / "herwig" / f"{self.sample.name}.in").write_text("changed card\n")
                elif change == "fingerprint":
                    fixture[1]["fingerprint"] = "0" * 64
                else:
                    fixture[0].ebeam = 6500
                with self.assertRaises(ValueError):
                    self.invoke(fixture)
                self.assertFalse((fixture[4] / "shower-recovery").exists())

    def test_changed_lhe_is_rejected_before_mutation(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            with patch.object(ho, "inspect_lhe", return_value={"events": 500}):
                with self.assertRaisesRegex(ValueError, "LHE changed"):
                    ho.recover_failed_shower(fixture[0], self.sample, fixture[1], fixture[2], fingerprint(fixture[2]))
            self.assertTrue(fixture[7].exists())
            self.assertFalse((fixture[4] / "shower-recovery").exists())

    def test_manifest_publication_failure_restores_originals_and_archives_replacement(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            manifest_bytes = fixture[6].read_bytes()
            root_bytes = fixture[7].read_bytes()
            original_writer = ho.write_json

            def fail_after_write(path, data, args):
                original_writer(path, data, args)
                raise OSError("injected publication failure")

            with patch.object(ho, "write_json", side_effect=fail_after_write):
                with self.assertRaisesRegex(OSError, "injected publication"):
                    self.invoke(fixture)
            self.assertEqual(fixture[6].read_bytes(), manifest_bytes)
            self.assertEqual(fixture[7].read_bytes(), root_bytes)
            self.assertEqual((fixture[4] / "herwig" / f"{self.sample.name}.in").read_text(),
                             fixture[1]["configuration"]["herwig_card"])
            archives = list((fixture[4] / "shower-recovery").iterdir())
            self.assertEqual(len(archives), 1)
            self.assertTrue((archives[0] / "replacement-herwig" / f"{self.sample.name}.in").is_file())
            self.assertTrue((archives[0] / "replacement-campaign.json").is_file())


if __name__ == "__main__":
    unittest.main()
