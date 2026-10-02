"""Persisted shower evidence and signed source populations; no generator needed."""

import copy
import gzip
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from hgammagamma import ho_shower_completion as completion


SAMPLE = "signal_gg_h_aa"
CARD = """create ThePEG::LesHouchesFileReader LesHouchesReader
set LesHouchesReader:AllowedToReOpen No
set LesHouchesHandler:WeightOption VarNegWeight
set LesHouchesHandler:Cuts /Herwig/Cuts/NoCuts
set LesHouchesReader:Cuts /Herwig/Cuts/NoCuts
insert LesHouchesHandler:LesHouchesReaders 0 LesHouchesReader
set /Herwig/Analysis/HwSim:OnTheFlyAnalysis No
"""


def lhe_text(weights=(10., -2., 8.), idwtup=-4):
    events = "".join(
        f"<event>\n1 1 {weight} 100 .007 .118\n"
        "25 1 0 0 0 0 0 0 0 125 125 0 9\n"
        "<rwgt><wgt id='1'>1.0</wgt></rwgt>\n</event>\n"
        for weight in weights)
    return ('<LesHouchesEvents version="3.0">\n<header>test input</header>\n'
            f"<init>\n2212 2212 20000 20000 -1 -1 336100 336100 {idwtup} 1\n"
            "5.333333333333333 .2 9999 1\n</init>\n" + events + "</LesHouchesEvents>\n")


class ShowerCompletionTests(unittest.TestCase):
    def fixture(self, directory, *, weights=(10., -2., 8.), saved=None,
                requested=None, termination="requested_events", attempts=None):
        directory = Path(directory)
        herwig = directory / "herwig"
        (herwig / "events").mkdir(parents=True)
        source = directory / "source.lhe"
        source.write_text(lhe_text(weights))
        attempts = len(weights) if attempts is None else attempts
        saved = attempts if saved is None else saved
        requested = saved if requested is None else requested
        discarded = attempts - saved
        out = ("Statistics for Les Houches event handler 'LesHouchesHandler':\n"
               "generated number of attempts Cross-section (nb)\n"
               f"Total: {saved} {attempts} 0.214(1)e+00\n"
               "Per Les Houches Reader breakdown:\n"
               f"LesHouchesReader {saved} {attempts} 0.214(1)e+00\n"
               "Statistics for unrelated UE process:\nTotal: 100291 100291 1.7e+08\n")
        (herwig / f"{SAMPLE}.out").write_text(out)
        root = herwig / "events" / f"{SAMPLE}.root"
        root.write_bytes(b"finalized ROOT fixture")
        normalized = [weight / max(map(abs, weights)) for weight in weights[:saved]]
        metadata = [{"path": str(root.resolve()), "entries": saved,
                     "size": root.stat().st_size, "mtime_ns": root.stat().st_mtime_ns,
                     "sum_weight": math.fsum(normalized),
                     "sum_abs_weight": math.fsum(map(abs, normalized)),
                     "sum_weight_squared": math.fsum(weight * weight for weight in normalized)}]
        run = ("HwSim Loaded.\nBasicConsistency: maximum 4-momentum violation: 4120632.432 MeV\n"
               f"A root tree has been written to the file:events/./{SAMPLE}.root\n"
               f"Number of events that pass basic cuts: {saved}\n"
               f"Weight of events that pass basic cuts: {metadata[0]['sum_weight']:.6g}\n")
        classes = ["ThePEG::EventHandler::ConsistencyException warning (206 times)"]
        if discarded:
            classes.append(f"ThePEG::Exception eventerror ({discarded} times)")
        if termination == "source_exhausted":
            run += ("Herwig: ThePEG::Exception caught.\n" + completion.EOF_MESSAGE
                    + " LesHouchesReader\nSee logfile for details.\n")
            classes.append("ThePEG::Exception runerror (1 times)")
        (herwig / "run.log").write_text(run)
        (herwig / f"{SAMPLE}.log").write_text(
            "Some detailed diagnostics\nThe following exception classes were reported in this run:\n"
            + "\n".join(classes) + "\n\nMiscellaneous output from modules:\n")
        card = CARD + f"set LesHouchesReader:FileName {source.resolve()}\n"
        (herwig / f"{SAMPLE}.in").write_text(card)
        manifest = {"sample": SAMPLE, "nevents_requested": requested,
                    "configuration": {"herwig_card": card},
                    "root_files": [{key: entry[key] for key in ("path", "size", "mtime_ns")}
                                   for entry in metadata],
                    "lhe_file": str(source.resolve()),
                    "lhe": {"events": len(weights), "source": str(source.resolve()),
                            "source_size": source.stat().st_size,
                            "source_mtime_ns": source.stat().st_mtime_ns}}
        return herwig, source, metadata, manifest, termination, attempts

    def build(self, fixture):
        herwig, source, metadata, manifest, termination, attempts = fixture
        weights = completion.source_weight_summary(source, attempts)
        record = completion.build_completion(
            herwig, SAMPLE, manifest["nevents_requested"], manifest["lhe"]["events"],
            metadata, termination, source_weights=weights)
        manifest["shower_completion"] = record
        return record, manifest

    def test_normal_completion_is_sealed_and_reuses_saved_count(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            record, manifest = self.build(fixture)
            self.assertEqual(record["saved_events"], 3)
            self.assertEqual(record["discarded_events"], 0)
            self.assertEqual(record["max_momentum_violation_mev"], 4120632.432)
            self.assertEqual(record["exception_counts"]["warning"], 206)
            self.assertEqual(completion.expected_analysis_events(manifest), 3)
            self.assertIs(completion.validate_completion(record, manifest, directory), record)
            self.assertAlmostEqual(completion.normalization_denominator(manifest, 1.6 * .00227, .00227),
                                   1.6 * .00227)

    def test_exhausted_source_retains_failed_signed_weights_in_denominator(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory, weights=(10., 10., -10., 10., 10., 10.),
                                   saved=4, requested=6, termination="source_exhausted")
            record, manifest = self.build(fixture)
            self.assertEqual(record["discarded_events"], 2)
            self.assertEqual(completion.expected_analysis_events(manifest), 4)
            self.assertEqual(record["source_weights"]["sum_weight"], 4.)
            self.assertEqual(record["root_summary"]["sum_weight"], 2.)
            self.assertAlmostEqual(completion.normalization_denominator(manifest, .00454, .00227), .00908)
            self.assertIs(completion.validate_completion(record, manifest, directory), record)

    def test_discarded_positive_and_negative_weights_do_not_use_count_correction(self):
        with TemporaryDirectory() as directory:
            weights = (10., 10., 10., 10., 10., -10., 10., 10., 10., 10., 10., -10.)
            fixture = self.fixture(directory, weights=weights, saved=6, requested=12,
                                   termination="source_exhausted")
            record, manifest = self.build(fixture)
            self.assertEqual(record["root_summary"]["sum_weight"], 4.)
            self.assertEqual(record["source_weights"]["sum_weight"], 8.)
            self.assertEqual(record["discarded_events"], 6)
            self.assertEqual(completion.normalization_denominator(manifest, 4., 1.), 8.)
            self.assertEqual(record["source_weights"]["negative_events"], 2)

    def test_source_prefix_uses_full_file_maximum_and_gzip_hash(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.lhe.gz"
            with gzip.open(path, "wt") as stream:
                stream.write(lhe_text((10., -2., 8., 20.)))
            summary = completion.source_weight_summary(path, 3)
            self.assertEqual(summary["lhe_events"], 4)
            self.assertEqual(summary["effective_max_weight"], 20.)
            self.assertAlmostEqual(summary["sum_weight"], .8)
            self.assertAlmostEqual(summary["sum_abs_weight"], 1.)
            self.assertAlmostEqual(summary["sum_weight_squared"], .42)
            self.assertEqual(summary["artifact"]["path"], str(path.resolve()))

    def test_source_population_supports_three_and_four_only(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.lhe"
            for idwtup in (-4, -3, 3, 4):
                with self.subTest(idwtup=idwtup):
                    path.write_text(lhe_text(idwtup=idwtup))
                    self.assertEqual(completion.source_weight_summary(path, 3)["idwtup"], idwtup)
            path.write_text(lhe_text(idwtup=-1))
            with self.assertRaisesRegex(ValueError, "IDWTUP"):
                completion.source_weight_summary(path, 3)

    def test_source_rejects_unclosed_nonfinite_zero_and_invalid_events(self):
        cases = {
            "unclosed": lhe_text().replace("</LesHouchesEvents>", ""),
            "zero": lhe_text((10., 0., 8.)),
            "nan_weight": lhe_text((10., float("nan"), 8.)),
            "nonfinite_beam": lhe_text().replace("20000 20000", "inf 20000"),
            "wrong_nup": lhe_text().replace("1 1 10.0", "2 1 10.0"),
            "undeclared": lhe_text().replace("1 1 10.0", "1 2 10.0"),
            "bad_particle": lhe_text().replace("25 1 0 0 0 0 0 0 0 125 125 0 9", "25 1 0 0 0 0 0"),
            "fractional_nup": lhe_text().replace("1 1 10.0", "1.5 1 10.0"),
            "multiple_processes": lhe_text().replace("-4 1\n5.333", "-4 2\n5.333"),
        }
        with TemporaryDirectory() as directory:
            path = Path(directory) / "source.lhe"
            for name, text in cases.items():
                with self.subTest(name=name):
                    path.write_text(text)
                    with self.assertRaises(ValueError):
                        completion.source_weight_summary(path, 3)
            path.write_text(lhe_text())
            with self.assertRaises(ValueError):
                completion.source_weight_summary(path, 4)

    def test_provisional_exhausted_record_requires_attached_source_proof(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory, weights=(10., 10., -10., 10., 10., 10.),
                                   saved=4, requested=6, termination="source_exhausted")
            herwig, source, metadata, manifest, termination, attempts = fixture
            record = completion.build_completion(herwig, SAMPLE, 6, 6, metadata, termination)
            with self.assertRaisesRegex(ValueError, "source-weight proof"):
                completion.validate_completion(record, manifest)
            record["source_weights"] = completion.source_weight_summary(source, attempts)
            record["fingerprint"] = completion.completion_fingerprint(record)
            self.assertIs(completion.validate_completion(record, manifest, directory), record)

    def test_normal_completion_can_account_for_retried_event_errors(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory, weights=(10., 10., -10., 10., 10., 10.),
                                   saved=4, requested=4, termination="requested_events", attempts=6)
            record, _ = self.build(fixture)
            self.assertEqual(record["attempted_events"], 6)
            self.assertEqual(record["generated_events"], 4)
            self.assertEqual(record["exception_counts"]["eventerror"], 2)

    def test_exhaustion_rejects_other_fatal_errors_and_unaccounted_loss(self):
        changes = {
            "extra_runerror": (f"{SAMPLE}.log", "runerror (1 times)", "runerror (2 times)"),
            "fatal_setup": (f"{SAMPLE}.log", "warning (206 times)", "setuperror (206 times)"),
            "wrong_discard_count": (f"{SAMPLE}.log", "eventerror (2 times)", "eventerror (1 times)"),
            "missing_summary": (f"{SAMPLE}.log", "The following exception classes were reported in this run:", "Summary absent:"),
            "non_eof_error": ("run.log", completion.EOF_MESSAGE, "Unrelated Herwig run failure"),
            "post_eof_error": ("run.log", "See logfile for details.", "See logfile for details.\nSegmentation fault"),
            "missing_writer": ("run.log", "A root tree has been written to the file:", "Writer unfinished:"),
            "footer_mismatch": ("run.log", "basic cuts: 4", "basic cuts: 3"),
            "not_exhausted": (f"{SAMPLE}.out", "Total: 4 6", "Total: 4 5"),
            "replayed": (f"{SAMPLE}.log", "Some detailed diagnostics", "Reopening LesHouchesReader"),
            "nonfinite_xsec": (f"{SAMPLE}.out", "0.214(1)e+00", "nan"),
        }
        for name, (filename, old, new) in changes.items():
            with self.subTest(name=name), TemporaryDirectory() as directory:
                fixture = self.fixture(directory, weights=(10., 10., -10., 10., 10., 10.),
                                       saved=4, requested=6, termination="source_exhausted")
                path = fixture[0] / filename
                path.write_text(path.read_text().replace(old, new))
                with self.assertRaises(ValueError):
                    self.build(fixture)

    def test_root_count_and_signed_moment_mismatches_are_rejected(self):
        for name in ("entries", "sum_weight", "sum_abs_weight", "sum_weight_squared"):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                fixture[2][0][name] += 1
                with self.assertRaises(ValueError):
                    self.build(fixture)

    def test_analysis_closure_and_artifact_changes_are_detected(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            record, manifest = self.build(fixture)
            with self.assertRaisesRegex(ValueError, "analysis weights"):
                completion.normalization_denominator(manifest, .00227, .00227)
            root = Path(record["root_metadata"][0]["path"])
            root.write_bytes(b"changed ROOT fixture")
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                completion.validate_completion(record, manifest, directory)

    def test_source_file_provenance_and_campaign_identity_are_checked(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            record, manifest = self.build(fixture)
            altered = copy.deepcopy(manifest)
            altered["nevents_requested"] = 4
            with self.assertRaisesRegex(ValueError, "campaign identity"):
                completion.validate_completion(record, altered)
            altered = copy.deepcopy(manifest)
            altered["lhe"]["events"] = 4
            with self.assertRaisesRegex(ValueError, "LHE event count"):
                completion.validate_completion(record, altered)
            altered = copy.deepcopy(manifest)
            altered["lhe_file"] = str(Path(directory) / "different.lhe")
            with self.assertRaisesRegex(ValueError, "different LHE"):
                completion.validate_completion(record, altered)
            fixture[1].write_text(lhe_text((10., -3., 8.)))
            with self.assertRaisesRegex(ValueError, "artifact changed"):
                completion.validate_completion(record, manifest, directory)

    def test_reweighting_card_blocks_source_weight_proof(self):
        changes = ("set LesHouchesHandler:NormalizeWeights Yes\n",
                   "set LesHouchesReader:MaxScan 100\n",
                   "set LesHouchesReader:MaxFactor 2\n",
                   "set LesHouchesReader:ReweightPDF Yes\n",
                   "set EventGenerator:WeightNormalization CrossSection\n",
                   "insert LesHouchesReader:Reweights 0 /Herwig/Reweights/Foo\n",
                   "create ThePEG::LesHouchesFileReader AnotherReader\n")
        for change in changes:
            with self.subTest(change=change), TemporaryDirectory() as directory:
                fixture = self.fixture(directory)
                card = fixture[0] / f"{SAMPLE}.in"
                card.write_text(card.read_text() + change)
                with self.assertRaises(ValueError):
                    self.build(fixture)

    def test_reader_filename_must_match_weight_source(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            card = fixture[0] / f"{SAMPLE}.in"
            card.write_text(card.read_text().replace(str(fixture[1].resolve()), "/tmp/different.lhe"))
            with self.assertRaisesRegex(ValueError, "FileName differs"):
                self.build(fixture)

    def test_campaign_card_and_root_inventory_remain_bound_to_completion(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            record, manifest = self.build(fixture)
            changed = copy.deepcopy(manifest)
            changed["configuration"]["herwig_card"] += "# replacement card\n"
            with self.assertRaisesRegex(ValueError, "recorded campaign card"):
                completion.validate_completion(record, changed)
            changed = copy.deepcopy(manifest)
            changed["root_files"][0]["size"] += 1
            with self.assertRaisesRegex(ValueError, "ROOT inventory"):
                completion.validate_completion(record, changed)

    def test_no_exceptions_footer_and_legacy_fallback(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            (fixture[0] / f"{SAMPLE}.log").write_text("No exceptions reported in this run.\n")
            record, manifest = self.build(fixture)
            self.assertEqual(sum(record["exception_counts"].values()), 0)
            del manifest["shower_completion"]
            self.assertEqual(completion.expected_analysis_events(manifest), 3)
            self.assertEqual(completion.normalization_denominator(manifest, .91, .00227), .91)

    def test_tampered_record_and_nonfinite_denominator_are_rejected(self):
        with TemporaryDirectory() as directory:
            fixture = self.fixture(directory)
            record, manifest = self.build(fixture)
            record["saved_events"] = 2
            with self.assertRaisesRegex(ValueError, "fingerprint changed"):
                completion.expected_analysis_events(manifest)
            for value in (0., float("nan"), float("inf")):
                with self.subTest(value=value), self.assertRaises(ValueError):
                    completion.normalization_denominator({"nevents_requested": 3}, value, 1.)


if __name__ == "__main__":
    unittest.main()
