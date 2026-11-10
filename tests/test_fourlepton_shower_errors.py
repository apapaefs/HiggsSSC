"""Recovery must preserve production and change only the shower error limit."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from hfourlepton import recover_shower_errors as fix


class RecoveryTests(unittest.TestCase):
    def fixture(self, directory):
        run = Path(directory)
        source = run / "events.lhe.gz"
        source.write_bytes(b"completed production, never regenerate")
        marker = run / "production.json"
        marker.write_text("completed provenance")
        card = "set EventGenerator:MaxErrors 200\nset Random:Seed 321\n"
        rows = {name: dict(source=fix.four.artifact_stage_record(source),
            production={**fix.four.artifact_stage_record(source),
                        "marker": fix.four.artifact_stage_record(marker)},
            old_card=card, card=fix.change_budget(card, 1000), max_errors=1000)
            for name in fix.SAMPLES}
        receipt = dict(run_dir=str(run), samples=rows)
        fix.write_json(run / fix.RECEIPT, receipt)
        return run, source, card, receipt

    def test_only_error_budget_changes(self):
        card = "set EventGenerator:MaxErrors 200\nset X:Seed 300\n"
        self.assertEqual(fix.change_budget(card, 1000), card.replace("Errors 200", "Errors 1000"))
        for bad, budget in ((card, 200), (card, 199), (card + card, 1000), ("", 1000)):
            with self.assertRaises(ValueError):
                fix.change_budget(bad, budget)

    def test_context_prevents_generation_and_routes_isolated_workers(self):
        with TemporaryDirectory() as directory:
            run, source, card, receipt = self.fixture(directory)
            sample = SimpleNamespace(name="bkg_qq_4e")
            with patch.object(fix.ho, "herwig_card", return_value=card), \
                    patch.object(fix.four, "obtain_production_lhe") as generator, \
                    patch.object(fix.runtime.subprocess, "Popen") as popen:
                with fix.recovery_runtime(receipt):
                    self.assertEqual(fix.four.obtain_production_lhe(sample, [], None), source.resolve())
                    self.assertEqual(fix.ho.herwig_card(sample, source, 321, 1000000, None),
                                     card.replace("Errors 200", "Errors 1000"))
                    fix.runtime.subprocess.Popen(["python", str(Path(fix.runtime.__file__).resolve()),
                                                  "--worker", "job.json"])
                    self.assertEqual(popen.call_args.args[0],
                        ["python", str(Path(fix.__file__).resolve()), "worker", "job.json"])
                    fix.runtime.subprocess.Popen(["Herwig", "run", "sample.run"])
                    self.assertEqual(popen.call_args.args[0], ["Herwig", "run", "sample.run"])
                generator.assert_not_called()
                self.assertIs(fix.four.obtain_production_lhe, generator)
            self.assertEqual(source.read_bytes(), b"completed production, never regenerate")

    def test_context_refuses_changed_hard_input_or_card(self):
        for corruption in ("source", "marker", "card", "new_card"):
            with self.subTest(corruption=corruption), TemporaryDirectory() as directory:
                run, source, card, receipt = self.fixture(directory)
                sample = SimpleNamespace(name="bkg_qq_4e")
                if corruption in ("source", "marker"):
                    (source if corruption == "source" else run / "production.json").write_text("changed")
                if corruption == "new_card":
                    receipt["samples"][sample.name]["card"] += "set Random:Seed 999\n"
                with patch.object(fix.ho, "herwig_card", return_value=card + ("changed" if corruption == "card" else "")):
                    with fix.recovery_runtime(receipt), self.assertRaises(ValueError):
                        if corruption in ("source", "marker"):
                            fix.four.obtain_production_lhe(sample, [], None)
                        else:
                            fix.ho.herwig_card(sample, source, 321, 1000000, None)

    def test_prepare_archives_all_failed_files_without_touching_generation(self):
        with TemporaryDirectory() as directory:
            run, source, _, receipt = self.fixture(directory)
            (run / fix.RECEIPT).unlink()
            (run / "manifest_perfect.json").write_text("original manifest")
            (run / "workers").mkdir()
            (run / "workers/failed.log").write_text("original traceback")
            for name in fix.SAMPLES:
                path = run / "samples/shared" / name / "herwig/events"
                path.mkdir(parents=True)
                (path / "partial.root").write_bytes(b"failed output")
            with patch.object(fix, "audit", return_value=receipt):
                result = fix.prepare(run, 1000)
            for name in fix.SAMPLES:
                self.assertFalse((run / "samples/shared" / name / "herwig").exists())
                self.assertEqual((Path(result["archive"]) / name / "events/partial.root").read_bytes(), b"failed output")
            self.assertEqual((run / "manifest_perfect.json").read_text(), "original manifest")
            self.assertEqual(source.read_bytes(), b"completed production, never regenerate")
            with self.assertRaisesRegex(ValueError, "already prepared"):
                fix.prepare(run, 1000)

    def test_prepare_failure_rolls_back_moved_showers(self):
        with TemporaryDirectory() as directory:
            run, _, _, receipt = self.fixture(directory)
            (run / fix.RECEIPT).unlink()
            (run / "manifest_perfect.json").write_text("original")
            (run / "workers").mkdir()
            for name in fix.SAMPLES:
                (run / "samples/shared" / name / "herwig").mkdir(parents=True)
            with patch.object(fix, "audit", return_value=receipt), \
                    patch.object(fix, "write_json", side_effect=OSError("disk failure")):
                with self.assertRaises(OSError):
                    fix.prepare(run, 1000)
            for name in fix.SAMPLES:
                self.assertTrue((run / "samples/shared" / name / "herwig").is_dir())


if __name__ == "__main__":
    unittest.main()
