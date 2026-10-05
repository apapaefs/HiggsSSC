"""Synthetic fresh-build provenance checks; no executable or integration runs."""

import copy
import unittest

from hgammagamma import ho_signal_normalization as norm
from test_ihixs_normalization import synthetic_record
from test_ihixs_parser_repair_record import synthetic_repair_record


def add_primary_parser_provenance(record):
    """Describe a synthetic parser-fixed primary build, separate from recovery."""
    descriptor = {
        "replacements": {"src/tools/user_interface.cpp": [
            ["synthetic undersized string", "synthetic terminated string"],
            ["synthetic out-of-bounds sentinel", "synthetic bounded sentinel"]]},
        "counts": {"src/tools/user_interface.cpp": [1, 4]},
    }
    digest = norm.canonical_sha256(descriptor)
    record["provenance"].update(parser_patch_sha256=digest, parser_patches=descriptor)
    record["validations"]["benchmark"]["parser_patch_sha256"] = digest
    for run in record["runs"]:
        # Existing mixed-build fixtures already bind repaired comparison runs.
        if "build_variant" not in run:
            run.update(build_variant="lhapdf",
                       executable_sha256=record["provenance"]["executable_sha256"])
    return norm.seal_record(record)


def synthetic_student_record():
    return add_primary_parser_provenance(synthetic_record())


class StudentParserRecordTests(unittest.TestCase):
    def assert_invalid(self, record):
        with self.assertRaises(ValueError):
            norm.validate_record(norm.seal_record(record))

    def test_fresh_primary_build_has_one_benchmark_and_no_legacy_recovery(self):
        record = synthetic_student_record()
        self.assertNotIn("parser_repair", record["provenance"])
        self.assertNotIn("parser_repair_benchmark", record["validations"])
        self.assertNotIn("parser_repair_central", record["validations"])
        self.assertEqual(len(record["runs"]), 109)
        self.assertIs(norm.validate_record(record), record)

    def test_historical_unpatched_and_mixed_recovery_records_remain_valid(self):
        for record in (synthetic_record(), synthetic_repair_record()):
            with self.subTest(mixed="parser_repair" in record["provenance"]):
                self.assertNotIn("parser_patch_sha256", record["provenance"])
                self.assertIs(norm.validate_record(record), record)

    def test_primary_patch_digest_descriptor_and_benchmark_are_all_required(self):
        for container, key in (("provenance", "parser_patch_sha256"),
                               ("provenance", "parser_patches"),
                               ("benchmark", "parser_patch_sha256")):
            with self.subTest(container=container, missing=key):
                record = synthetic_student_record()
                target = (record["validations"]["benchmark"] if container == "benchmark"
                          else record["provenance"])
                target.pop(key)
                self.assert_invalid(record)

    def test_patch_metadata_cannot_be_orphaned_on_a_historical_primary_build(self):
        descriptor = synthetic_student_record()["provenance"]["parser_patches"]
        for container, key, value in (
                ("provenance", "parser_patches", descriptor),
                ("benchmark", "parser_patch_sha256", norm.canonical_sha256(descriptor))):
            with self.subTest(container=container, orphan=key):
                record = synthetic_record()
                target = (record["validations"]["benchmark"] if container == "benchmark"
                          else record["provenance"])
                target[key] = copy.deepcopy(value)
                self.assert_invalid(record)

    def test_primary_and_benchmark_patch_digests_must_be_valid_and_match(self):
        for container in ("provenance", "benchmark"):
            for value in ("invalid-digest", "f" * 64, None):
                with self.subTest(container=container, value=value):
                    record = synthetic_student_record()
                    target = (record["validations"]["benchmark"] if container == "benchmark"
                              else record["provenance"])
                    target["parser_patch_sha256"] = value
                    self.assert_invalid(record)

    def test_descriptor_contents_are_covered_by_the_patch_digest(self):
        for field in ("replacements", "counts"):
            with self.subTest(field=field):
                record = synthetic_student_record()
                descriptor = record["provenance"]["parser_patches"]
                descriptor[field]["src/tools/user_interface.cpp"] = ["changed synthetic patch"]
                self.assert_invalid(record)

    def test_descriptor_contains_both_nonempty_replacement_and_count_mappings(self):
        for descriptor in ({}, {"replacements": {}, "counts": {}},
                           {"replacements": {"synthetic": []}},
                           {"replacements": [], "counts": {"synthetic": []}},
                           {"replacements": {"synthetic": []}, "counts": None}):
            with self.subTest(descriptor=descriptor):
                record = synthetic_student_record()
                digest = norm.canonical_sha256(descriptor)
                record["provenance"].update(parser_patches=descriptor, parser_patch_sha256=digest)
                record["validations"]["benchmark"]["parser_patch_sha256"] = digest
                self.assert_invalid(record)

    def test_fresh_primary_build_binds_every_run_to_its_executable(self):
        for index in (0, 6, 7, 106, 107, 108):
            for field in ("build_variant", "executable_sha256"):
                with self.subTest(index=index, missing=field):
                    record = synthetic_student_record()
                    record["runs"][index].pop(field)
                    self.assert_invalid(record)
        for field, value in (("build_variant", "lhapdf-parser-repair"),
                             ("executable_sha256", "f" * 64)):
            record = synthetic_student_record()
            record["runs"][-1][field] = value
            self.assert_invalid(record)

    def test_primary_benchmark_still_requires_the_correct_executable(self):
        record = synthetic_student_record()
        record["validations"]["benchmark"]["production_executable_sha256"] = "f" * 64
        self.assert_invalid(record)

    def test_primary_patch_metadata_does_not_waive_separate_recovery_validation(self):
        record = add_primary_parser_provenance(synthetic_repair_record())
        self.assertIs(norm.validate_record(record), record)
        for field in ("parser_repair_benchmark", "parser_repair_central"):
            with self.subTest(missing=field):
                incomplete = copy.deepcopy(record)
                incomplete["validations"].pop(field)
                self.assert_invalid(incomplete)


if __name__ == "__main__":
    unittest.main()
