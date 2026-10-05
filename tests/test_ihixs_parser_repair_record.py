"""Synthetic mixed-build records only; no executable or physics calculation."""

import copy
import unittest

from hgammagamma import ho_signal_normalization as norm
from test_ihixs_normalization import synthetic_record


def synthetic_repair_record():
    record = synthetic_record()
    provenance = record["provenance"]
    integration = {"epsrel": 5e-5, "epsabs": 0., "mineval": 200000,
                   "maxeval": 200000000, "nstart": 40000, "nincrease": 4000,
                   "cuba_verbose": 0}
    provenance["settings"].update(integration=copy.deepcopy(integration),
                                  numerical_relative_tolerance=.0005)
    provenance["settings_sha256"] = norm.canonical_sha256(provenance["settings"])
    repair = {"ihixs_commit": norm.IHIXS_COMMIT,
              "source_sha256": provenance["source_sha256"],
              "adapter_sha256": provenance["adapter_sha256"],
              "parser_patch_sha256": "b" * 64, "executable_sha256": "c" * 64,
              "build_manifest_sha256": "d" * 64}
    provenance["parser_repair"] = repair
    for run in record["runs"]:
        run.update(executable_sha256=provenance["executable_sha256"], build_variant="lhapdf")
    central = record["runs"][0]
    central["integration"] = copy.deepcopy(integration)
    for run in record["runs"][-2:]:
        run.update(executable_sha256=repair["executable_sha256"], build_variant="lhapdf-parser-repair")
    benchmark = copy.deepcopy(record["validations"]["benchmark"])
    benchmark.update(executable_sha256="e" * 64,
                     production_executable_sha256=repair["executable_sha256"],
                     parser_patch_sha256=repair["parser_patch_sha256"])
    record["validations"]["parser_repair_benchmark"] = benchmark
    sigma, error = 100.001, .01
    multiplier = record["prefactor_correction"]["multiplier"]
    record["validations"]["parser_repair_central"] = {
        "passed": True, "relative_tolerance": .0005,
        "relative_difference": abs(sigma / central["cross_section_pb"] - 1.),
        "original_cross_section_pb": central["cross_section_pb"],
        "original_cross_section_error_pb": central["cross_section_error_pb"],
        "cross_section_pb": sigma, "cross_section_error_pb": error,
        "result_key": "eftn3lo", "qcd_order": "N3LO", "pdf_set": norm.PDF_SET,
        "pdf_member": 0, "mur_gev": 62.5, "muf_gev": 62.5,
        "alpha_s_mur": .125, "alpha_s_at_91_1876": .118, "pdf_alpha_s_mur": .125,
        "input_sha256": "1" * 64, "output_sha256": "2" * 64,
        "pdf_info_sha256": central["pdf_info_sha256"],
        "pdf_member_sha256": central["pdf_member_sha256"],
        "executable_sha256": repair["executable_sha256"], "build_variant": "lhapdf-parser-repair",
        "run_directory": "parser-repair/runs/central_check", "precision_refinements": 0,
        "integration": copy.deepcopy(integration), "raw_cross_section_pb": sigma / multiplier,
        "raw_cross_section_error_pb": error / multiplier}
    return norm.seal_record(record)


class ParserRepairRecordTests(unittest.TestCase):
    def assert_invalid(self, record):
        with self.assertRaises(ValueError):
            norm.validate_record(norm.seal_record(record))

    def test_historical_records_without_per_run_executables_remain_valid(self):
        record = synthetic_record()
        for run in record["runs"]:
            run.pop("executable_sha256", None)
            run.pop("build_variant", None)
        record = norm.seal_record(record)
        self.assertIs(norm.validate_record(record), record)

    def test_original_executable_metadata_is_valid_without_repair(self):
        record = synthetic_record()
        for run in record["runs"]:
            run.update(executable_sha256=record["provenance"]["executable_sha256"], build_variant="lhapdf")
        record = norm.seal_record(record)
        self.assertIs(norm.validate_record(record), record)

    def test_mixed_build_record_keeps_original_central_and_benchmark(self):
        record = synthetic_repair_record()
        original_benchmark = copy.deepcopy(record["validations"]["benchmark"])
        self.assertIs(norm.validate_record(record), record)
        self.assertEqual(record["cross_section_pb"], 100.)
        self.assertEqual(record["validations"]["parser_repair_central"]["cross_section_pb"], 100.001)
        self.assertEqual(record["validations"]["benchmark"], original_benchmark)
        self.assertEqual(record["runs"][0]["build_variant"], "lhapdf")
        self.assertEqual(record["runs"][-2]["build_variant"], "lhapdf-parser-repair")

    def test_repaired_central_cannot_replace_the_original_build_baseline(self):
        record = synthetic_repair_record()
        central = record["runs"][0]
        central.update(build_variant="lhapdf-parser-repair",
                       executable_sha256=record["provenance"]["parser_repair"]["executable_sha256"])
        self.assert_invalid(record)

    def test_repair_requires_complete_hashes_and_same_source_adapter_commit(self):
        for key in ("source_sha256", "adapter_sha256", "parser_patch_sha256",
                    "executable_sha256", "build_manifest_sha256"):
            for missing in (False, True):
                with self.subTest(key=key, missing=missing):
                    record = synthetic_repair_record()
                    if missing:
                        record["provenance"]["parser_repair"].pop(key)
                    else:
                        record["provenance"]["parser_repair"][key] = "invalid-digest"
                    self.assert_invalid(record)
        for key, value in (("ihixs_commit", "different-commit"),
                           ("source_sha256", "f" * 64), ("adapter_sha256", "f" * 64)):
            with self.subTest(key=key):
                record = synthetic_repair_record()
                record["provenance"]["parser_repair"][key] = value
                self.assert_invalid(record)

    def test_each_run_is_bound_to_its_explicit_build_variant(self):
        for index, key, value in ((0, "build_variant", "unrecognized"),
                                  (0, "build_variant", "lhapdf-parser-repair"),
                                  (0, "executable_sha256", "c" * 64),
                                  (-2, "build_variant", "lhapdf"),
                                  (-2, "executable_sha256", "a" * 64),
                                  (-2, "executable_sha256", "invalid-digest")):
            with self.subTest(index=index, key=key, value=value):
                record = synthetic_repair_record()
                record["runs"][index][key] = value
                self.assert_invalid(record)
        for index in (0, -2):
            for key in ("executable_sha256", "build_variant"):
                with self.subTest(index=index, missing=key):
                    record = synthetic_repair_record()
                    record["runs"][index].pop(key)
                    self.assert_invalid(record)

    def test_partial_or_repaired_run_metadata_without_repair_is_rejected(self):
        for metadata in ({"build_variant": "lhapdf"}, {"executable_sha256": "a" * 64},
                         {"build_variant": "lhapdf", "executable_sha256": "c" * 64},
                         {"build_variant": "lhapdf-parser-repair", "executable_sha256": "c" * 64}):
            with self.subTest(metadata=metadata):
                record = synthetic_record()
                record["runs"][0].update(metadata)
                self.assert_invalid(record)

    def test_repair_benchmark_is_required_and_uses_same_published_accuracy_gates(self):
        record = synthetic_repair_record()
        record["validations"].pop("parser_repair_benchmark")
        self.assert_invalid(record)
        record = synthetic_repair_record()
        record["validations"]["parser_repair_benchmark"].pop("parser_patch_sha256")
        self.assert_invalid(record)
        for key, value in (("passed", False), ("ihixs_commit", "different-commit"),
                           ("production_executable_sha256", "a" * 64),
                           ("parser_patch_sha256", "f" * 64),
                           ("parser_patch_sha256", "invalid-digest"),
                           ("executable_sha256", "invalid-digest"), ("input_sha256", "invalid-digest"),
                           ("output_sha256", "invalid-digest"), ("result_key", "Higgs XS"),
                           ("pdf_set", norm.PDF_SET), ("expected_pb", 45.),
                           ("relative_tolerance", .01), ("raw_cross_section_pb", 50.),
                           ("raw_cross_section_error_pb", .03)):
            with self.subTest(key=key):
                record = synthetic_repair_record()
                record["validations"]["parser_repair_benchmark"][key] = value
                self.assert_invalid(record)

    def test_repair_validations_cannot_be_orphaned_from_build_provenance(self):
        for key in ("parser_repair_benchmark", "parser_repair_central"):
            with self.subTest(key=key):
                record = synthetic_record()
                record["validations"][key] = synthetic_repair_record()["validations"][key]
                self.assert_invalid(record)

    def test_repaired_central_check_is_required_and_consistent(self):
        record = synthetic_repair_record()
        record["validations"].pop("parser_repair_central")
        self.assert_invalid(record)
        for key, value in (("passed", False), ("original_cross_section_pb", 101.),
                           ("original_cross_section_error_pb", .02), ("relative_tolerance", .001),
                           ("relative_difference", 0.), ("cross_section_error_pb", .1),
                           ("result_key", "eftnnlo"), ("qcd_order", "NNLO"),
                           ("pdf_set", norm.NATIVE_PDF_SET), ("pdf_member", 1),
                           ("mur_gev", 31.25), ("muf_gev", 125.),
                           ("alpha_s_mur", .126), ("alpha_s_at_91_1876", .119),
                           ("pdf_alpha_s_mur", .126), ("executable_sha256", "a" * 64),
                           ("build_variant", "lhapdf"), ("pdf_info_sha256", "f" * 64),
                           ("pdf_member_sha256", "f" * 64), ("input_sha256", "invalid-digest"),
                           ("output_sha256", "invalid-digest"), ("raw_cross_section_pb", 100.),
                           ("raw_cross_section_error_pb", .02),
                           ("run_directory", "runs/central"), ("precision_refinements", True)):
            with self.subTest(key=key):
                record = synthetic_repair_record()
                record["validations"]["parser_repair_central"][key] = value
                self.assert_invalid(record)
        record = synthetic_repair_record()
        check = record["validations"]["parser_repair_central"]
        check.update(cross_section_pb=100.1, relative_difference=abs(100.1 / 100. - 1.),
                     raw_cross_section_pb=100.1 / record["prefactor_correction"]["multiplier"])
        self.assert_invalid(record)

    def test_central_check_refinement_metadata_matches_original_selected_precision(self):
        record = synthetic_repair_record()
        central = record["runs"][0]
        # The original central may itself already have been numerically refined.
        central["integration"]["epsrel"] = 1e-6
        check = record["validations"]["parser_repair_central"]
        check.update(precision_refinements=2,
                     run_directory="parser-repair/runs/central_check__precision_2")
        check["integration"] = copy.deepcopy(central["integration"])
        check["integration"]["epsrel"] = min(1e-6 / 10., 1e-6) / 10.
        record = norm.seal_record(record)
        self.assertIs(norm.validate_record(record), record)
        for key, value in (("precision_refinements", 4), ("run_directory", "parser-repair/runs/central_check")):
            bad = copy.deepcopy(record)
            bad["validations"]["parser_repair_central"][key] = value
            self.assert_invalid(bad)
        for key, value in (("epsrel", 5e-5), ("maxeval", 400000000)):
            bad = copy.deepcopy(record)
            bad["validations"]["parser_repair_central"]["integration"][key] = value
            self.assert_invalid(bad)


if __name__ == "__main__":
    unittest.main()
