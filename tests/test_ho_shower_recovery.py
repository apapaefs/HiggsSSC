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

    def write_error_limit_logs(self, sample_dir, *, eventerrors=100):
        herwig = sample_dir / "herwig"
        (herwig / f"{self.sample.name}.out").write_text(
            "Statistics for Les Houches event handler 'LesHouchesHandler':\n"
            f"Total: 378 {378 + eventerrors} 0.214(1)e+00\n"
            "Per Les Houches Reader breakdown:\n")
        (herwig / f"{self.sample.name}.log").write_text(
            "The following exception classes were reported in this run:\n"
            f"ThePEG::Exception eventerror ({eventerrors} times)\n\n")
        (herwig / "run.log").write_text(
            "A root tree has been written to the file:events/partial.root\n"
            "Number of events that pass basic cuts: 378\n"
            "Weight of events that pass basic cuts: 378\n"
            "Herwig: ThePEG::Exception caught.\n"
            "Too many tries in a hadron decay.\nSee logfile for details.\n")

    def corrupt_error_limit_logs(self, sample_dir, failure):
        herwig = sample_dir / "herwig"
        if failure == "missing":
            (herwig / "run.log").unlink()
        elif failure == "unrelated":
            (herwig / "run.log").write_text("Segmentation fault\n")
        elif failure == "inconsistent":
            path = herwig / f"{self.sample.name}.out"
            path.write_text(path.read_text().replace("Total: 378 478", "Total: 378 479"))
        elif failure == "malformed":
            path = herwig / f"{self.sample.name}.out"
            path.write_text(path.read_text().replace("0.214(1)e+00", " "))
        else:
            path = herwig / f"{self.sample.name}.log"
            path.write_text(path.read_text().replace("eventerror (100 times)", "abortnow (100 times)"))

    def test_error_limit_hint_recommends_scaled_or_doubled_budget_without_writes(self):
        for nevents, previous_budget, suggested in ((1000, 100, 200), (5000000, 100, 1000),
                                                     (5000000, 750, 1500)):
            with self.subTest(nevents=nevents, previous_budget=previous_budget), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                old, sample_dir = fixture[1], fixture[4]
                old["configuration"]["herwig_card"] = old["configuration"]["herwig_card"].replace(
                    "MaxErrors 100\n", f"MaxErrors {previous_budget}\n")
                old["nevents_requested"] = old["lhe"]["events"] = nevents
                args = self.arguments(directory, "--nevents", str(nevents), "--resume")
                self.write_error_limit_logs(sample_dir, eventerrors=previous_budget)
                before = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                hint = ho.error_limit_recovery_hint(args, self.sample, old)
                self.assertIn(f"--herwig-max-errors {suggested}", hint)
                self.assertIn("--stage prepare", hint)
                self.assertIn("--restart-failed-shower", hint)
                self.assertIn("--run-samples bkg_prompt_aa", hint)
                self.assertIn("from the beginning", hint)
                self.assertIn("completed generation is reused", hint)
                after = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                self.assertEqual(before, after)

    def test_resume_of_error_limit_abort_explains_recovery_and_preserves_outputs(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            old, sample_dir = fixture[1], fixture[4]
            args = self.arguments(directory, "--resume")
            self.write_error_limit_logs(sample_dir)
            before = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
            with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                    patch.object(ho, "command") as command, \
                    patch.object(ho, "record_shower_completion") as completion, \
                    patch.object(ho, "recover_failed_shower") as restart:
                with self.assertRaisesRegex(ValueError, "per-class event-error limit"):
                    ho.shower(args, self.sample, old)
            command.assert_not_called()
            completion.assert_not_called()
            restart.assert_not_called()
            after = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            self.assertFalse((sample_dir / "shower-recovery").exists())

    def test_new_error_limit_failure_explains_recovery_without_adopting_or_archiving(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            old, sample_dir, root = fixture[1], fixture[4], fixture[7]
            args = self.arguments(directory, "--resume")
            root.unlink()
            failure = RuntimeError("original generator error")
            stopped_files = {}

            def fail_run(_args, argv, _cwd, _log):
                if argv[1] == "run":
                    root.write_bytes(b"newly finalized partial ROOT")
                    self.write_error_limit_logs(sample_dir)
                    stopped_files.update({str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()})
                    raise failure

            with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                    patch.object(ho, "command", side_effect=fail_run), \
                    patch.object(ho, "record_shower_completion") as completion, \
                    patch.object(ho, "recover_failed_shower") as restart:
                with self.assertRaisesRegex(RuntimeError, "per-class event-error limit") as caught:
                    ho.shower(args, self.sample, old)
            self.assertIs(caught.exception.__cause__, failure)
            completion.assert_not_called()
            restart.assert_not_called()
            after = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
            self.assertEqual(stopped_files, after)
            self.assertFalse((sample_dir / "shower-recovery").exists())

    def test_unrecognized_failure_keeps_original_error_and_never_suggests_restart(self):
        for reason in ("unrelated", "missing", "inconsistent", "malformed", "fatal"):
            with self.subTest(reason=reason), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                old, sample_dir, root = fixture[1], fixture[4], fixture[7]
                args = self.arguments(directory, "--resume")
                root.unlink()
                failure = RuntimeError(f"original {reason} generator failure")

                def fail_run(_args, argv, _cwd, _log):
                    if argv[1] == "run":
                        root.write_bytes(b"unrecognized failure ROOT evidence")
                        self.write_error_limit_logs(sample_dir)
                        self.corrupt_error_limit_logs(sample_dir, reason)
                        raise failure

                with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                        patch.object(ho, "command", side_effect=fail_run), \
                        patch.object(ho, "record_shower_completion") as completion, \
                        patch.object(ho, "recover_failed_shower") as restart:
                    with self.assertRaises(RuntimeError) as caught:
                        ho.shower(args, self.sample, old)
                self.assertIs(caught.exception, failure)
                self.assertIsNone(ho.error_limit_recovery_hint(args, self.sample, old))
                completion.assert_not_called()
                restart.assert_not_called()
                self.assertEqual(root.read_bytes(), b"unrecognized failure ROOT evidence")
                self.assertFalse((sample_dir / "shower-recovery").exists())

    def test_unrecognized_resume_keeps_existing_completion_audit_error(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            old, sample_dir = fixture[1], fixture[4]
            args = self.arguments(directory, "--resume")
            self.write_error_limit_logs(sample_dir)
            self.corrupt_error_limit_logs(sample_dir, "inconsistent")
            failure = ValueError("original strict completion audit failure")
            with patch.object(ho, "inspect_lhe", return_value=copy.deepcopy(old["lhe"])), \
                    patch.object(ho, "record_shower_completion", side_effect=failure) as completion, \
                    patch.object(ho, "command") as command:
                with self.assertRaises(ValueError) as caught:
                    ho.shower(args, self.sample, old)
            self.assertIs(caught.exception, failure)
            completion.assert_called_once_with(args, self.sample, old, "source_exhausted")
            command.assert_not_called()

    def test_default_budget_preserves_existing_card_and_fingerprint(self):
        with TemporaryDirectory() as directory:
            args = self.arguments(directory)
            self.assertIsNone(args.herwig_max_errors)
            card = ho.herwig_card(args, self.sample)
            self.assertEqual(card.count("set EventGenerator:MaxErrors 100\n"), 1)
            old = ho.prepare(args, self.sample)
            explicit = self.arguments(directory, "--herwig-max-errors", "100")
            self.assertEqual(ho.herwig_card(explicit, self.sample), card)
            self.assertEqual(ho.prepare(explicit, self.sample), old)
            custom = self.arguments(directory, "--herwig-max-errors", "1000")
            self.assertEqual(ho.herwig_card(custom, self.sample).replace(
                "MaxErrors 1000\n", "MaxErrors 100\n"), card)

    def test_new_campaign_error_budget_scales_with_event_request_and_rounds_up(self):
        cases = ((100000, 100), (500000, 100), (500001, 101),
                 (1000000, 200), (1000001, 201), (5000000, 1000))
        for nevents, expected in cases:
            with self.subTest(nevents=nevents), TemporaryDirectory() as directory:
                args = self.arguments(directory, "--nevents", str(nevents))
                manifest = ho.prepare(args, self.sample)
                card = manifest["configuration"]["herwig_card"]
                self.assertEqual(card.count(f"set EventGenerator:MaxErrors {expected}\n"), 1)
                self.assertEqual((ho.sample_dir(args, self.sample) / "herwig" /
                                  f"{self.sample.name}.in").read_text(), card)
                self.assertEqual(manifest["fingerprint"], fingerprint(manifest["configuration"]))
                self.assertIsNone(args.herwig_max_errors)

    def test_explicit_budget_overrides_scaled_default_for_new_campaign(self):
        with TemporaryDirectory() as directory:
            args = self.arguments(directory, "--nevents", "5000000", "--herwig-max-errors", "237")
            manifest = ho.prepare(args, self.sample)
            self.assertIn("set EventGenerator:MaxErrors 237\n", manifest["configuration"]["herwig_card"])
            self.assertEqual(args.herwig_max_errors, 237)

    def test_existing_large_campaign_keeps_legacy_and_custom_budgets_on_resume(self):
        for previous_budget in (100, 777):
            with self.subTest(previous_budget=previous_budget), TemporaryDirectory() as directory:
                original = self.arguments(directory, "--nevents", "5000000",
                                          "--herwig-max-errors", str(previous_budget))
                previous = ho.prepare(original, self.sample)
                sample_dir = ho.sample_dir(original, self.sample)
                root = sample_dir / "herwig/events/retained.root"
                root.parent.mkdir()
                root.write_bytes(b"existing campaign ROOT must remain unchanged")
                previous["completed"] = ["generate", "shower"]
                ho.write_json(sample_dir / "campaign.json", previous, original)
                before = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                args = self.arguments(directory, "--nevents", "5000000", "--resume")
                resumed = ho.prepare(args, self.sample)
                self.assertEqual(resumed, previous)
                self.assertIn(f"set EventGenerator:MaxErrors {previous_budget}\n",
                              ho.herwig_card(args, self.sample))
                after = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                self.assertEqual(before, after)
                self.assertIsNone(args.herwig_max_errors)

    def test_common_resume_resolves_budget_per_sample_and_scales_only_new_samples(self):
        with TemporaryDirectory() as directory:
            existing = {}
            for sample, budget in ((ho.SAMPLES[1], 100), (ho.SAMPLES[2], 777)):
                original = self.arguments(directory, "--nevents", "5000000",
                                          "--herwig-max-errors", str(budget))
                existing[sample.name] = ho.prepare(original, sample)
            args = self.arguments(directory, "--nevents", "5000000", "--resume")
            for sample in (ho.SAMPLES[1], ho.SAMPLES[2]):
                self.assertEqual(ho.prepare(args, sample), existing[sample.name])
            fresh = ho.prepare(args, ho.SAMPLES[3])
            self.assertIn("set EventGenerator:MaxErrors 1000\n", fresh["configuration"]["herwig_card"])
            self.assertIsNone(args.herwig_max_errors)

    def test_invalid_persisted_error_budget_is_rejected_without_overwriting_campaign(self):
        settings = ("", "set EventGenerator:MaxErrors 0\n", "set EventGenerator:MaxErrors -1\n",
                    "set EventGenerator:MaxErrors 100.0\n", "set EventGenerator:MaxErrors nan\n",
                    "set EventGenerator:MaxErrors 100\nset EventGenerator:MaxErrors 100\n")
        for setting in settings:
            with self.subTest(setting=setting), TemporaryDirectory() as directory:
                args = self.arguments(directory)
                manifest = ho.prepare(args, self.sample)
                sample_dir = ho.sample_dir(args, self.sample)
                card = manifest["configuration"]["herwig_card"].replace(
                    "set EventGenerator:MaxErrors 100\n", setting)
                manifest["configuration"]["herwig_card"] = card
                manifest["fingerprint"] = fingerprint(manifest["configuration"])
                (sample_dir / "herwig" / f"{self.sample.name}.in").write_text(card)
                ho.write_json(sample_dir / "campaign.json", manifest, args)
                before = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                with self.assertRaises(ValueError):
                    ho.prepare(self.arguments(directory, "--resume"), self.sample)
                after = {str(path): path.read_bytes() for path in sample_dir.rglob("*") if path.is_file()}
                self.assertEqual(before, after)

    def test_restart_cli_requires_prepare_resume_and_positive_budget(self):
        cases = (("--herwig-max-errors", "0"), ("--herwig-max-errors", "-1"),
                 ("--herwig-max-errors", "1.5"), ("--restart-failed-shower",),
                 ("--restart-failed-shower", "--resume"),
                 ("--restart-failed-shower", "--resume", "--herwig-max-errors", "1000", "--stage", "shower"),
                 ("--restart-failed-shower", "--resume", "--herwig-max-errors", "1000", "--stage", "all"))
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
