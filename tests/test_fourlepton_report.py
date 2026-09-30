from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

from hfourlepton.make_fourlepton_report import (
    DEFAULT_CUT_STAGES,
    RootRecords,
    WeightStats,
    _comparison_identity,
    cutflow_stats,
    generate_report,
    grouped_weight_stats,
    mass_spectrum,
    optional_weight_audit,
    response_closure_audit,
)


SELECTED_MASK = (1 << 10) - 1


class FourLeptonReportHelpersTest(unittest.TestCase):
    def test_signed_weight_stats_keep_negative_and_zero_totals(self) -> None:
        stats = WeightStats()
        for value in (2.0, -3.0, 1.0):
            stats.add(value)

        result = stats.as_dict()
        self.assertEqual(result["sumw"], 0.0)
        self.assertEqual(result["sumabsw"], 6.0)
        self.assertEqual(result["sumw2"], 14.0)
        self.assertEqual(result["sumw_positive"], 3.0)
        self.assertEqual(result["sumw_negative"], -3.0)
        self.assertAlmostEqual(result["negative_weight_fraction"], 1.0 / 3.0)
        self.assertEqual(result["negative_absolute_weight_fraction"], 0.5)

    def test_sumw2_groups_exclusive_hypotheses_by_source(self) -> None:
        records = [
            {"source_index": 10, "event_weight": 0.4},
            {"source_index": 10, "event_weight": 0.6},
            {"source_index": 11, "event_weight": -0.25},
        ]
        stats = grouped_weight_stats(records)
        self.assertAlmostEqual(stats.sumw, 0.75)
        self.assertAlmostEqual(stats.sumw2, 1.0 + 0.25**2)
        self.assertEqual(stats.entries, 2)

    def test_cut_masks_are_cumulative_and_pretrigger_weight_is_used(self) -> None:
        self.assertEqual(DEFAULT_CUT_STAGES[-1].required_mask, SELECTED_MASK)
        records = [
            {
                "source_index": 1,
                "cut_mask": SELECTED_MASK,
                "pretrigger_event_weight": 1.0,
                "event_weight": 0.98,
            },
            {
                "source_index": 2,
                "cut_mask": (1 << 6) - 1,
                "pretrigger_event_weight": -0.5,
                "event_weight": -0.495,
            },
        ]
        cutflow = cutflow_stats(records)
        self.assertAlmostEqual(cutflow[0]["sumw"], 0.5)
        self.assertAlmostEqual(cutflow[5]["sumw"], 0.5)
        self.assertAlmostEqual(cutflow[6]["sumw"], 1.0)
        self.assertAlmostEqual(cutflow[-1]["sumw"], 0.98)
        self.assertEqual([row["bit"] for row in cutflow], list(range(10)))

    def test_response_closure_uses_failure_plus_hypothesis_probabilities(self) -> None:
        sources = [
            {
                "source_index": 1,
                "response_failure_probability": 0.2,
                "response_closure_delta": 0.0,
            },
            {
                "source_index": 2,
                "response_failure_probability": 1.0,
                "response_closure_delta": 0.0,
            },
        ]
        rows = [
            {"source_index": 1, "hypothesis_probability": 0.3},
            {"source_index": 1, "hypothesis_probability": 0.5},
        ]
        closure = response_closure_audit(sources, rows)
        self.assertTrue(closure["ok"])
        self.assertEqual(closure["evaluated"], 2)
        self.assertEqual(closure["max_abs_delta"], 0.0)

    def test_mass_spectrum_preserves_signed_bins_and_source_grouped_sumw2(self) -> None:
        records = [
            {"source_index": 1, "event_weight": 0.6, "m4l": 124.0},
            {"source_index": 1, "event_weight": -0.1, "m4l": 124.5},
            {"source_index": 2, "event_weight": -0.3, "m4l": 123.0},
            {"source_index": 3, "event_weight": 2.0, "m4l": 250.0},
        ]
        normalization = {"status": "ok", "scale_pb_per_weight": 1.0}
        spectrum = mass_spectrum(records, normalization, 100.0)
        mass_bin = next(
            value
            for value in spectrum["bins"]
            if value["low_gev"] == 120.0 and value["high_gev"] == 125.0
        )
        self.assertAlmostEqual(mass_bin["sumw"], 0.2)
        self.assertAlmostEqual(mass_bin["sumw2"], 0.5**2 + 0.3**2)
        self.assertAlmostEqual(mass_bin["events"], 20_000.0)
        self.assertAlmostEqual(spectrum["bins"][-1]["events"], 200_000.0)

    def test_optional_variations_keep_negative_raw_sums_and_check_names(self) -> None:
        sources = [
            {
                "source_index": 1,
                "optional_weight_names": ["scale_up", "scale_down"],
                "optional_weights": [1.2, 0.8],
                "optional_selected_event_weights": [0.6, 0.4],
            },
            {
                "source_index": 2,
                "optional_weight_names": ["scale_up", "scale_down"],
                "optional_weights": [-0.6, -0.4],
                "optional_selected_event_weights": [-0.3, -0.2],
            },
        ]
        audit = optional_weight_audit(sources)
        self.assertEqual(audit["status"], "ok_raw_unnormalized")
        self.assertFalse(audit["normalized"])
        self.assertAlmostEqual(audit["variations"][0]["generated_raw"]["sumw"], 0.6)
        self.assertAlmostEqual(
            audit["variations"][0]["generated_raw"]["sumw2"], 1.2**2 + 0.6**2
        )
        self.assertAlmostEqual(audit["variations"][0]["selected_raw"]["sumw"], 0.3)

        inconsistent = optional_weight_audit(
            [
                sources[0],
                {
                    **sources[1],
                    "optional_weight_names": ["scale_down", "scale_up"],
                },
            ]
        )
        self.assertEqual(inconsistent["status"], "inconsistent")
        self.assertTrue(inconsistent["errors"])

    def test_manifest_defined_optional_variations_reject_synthesized_names(self) -> None:
        audit = optional_weight_audit(
            [
                {
                    "source_index": 1,
                    "optional_weight_names": ["optional_0", "optional_1"],
                    "optional_weights": [1.2, 0.8],
                    "optional_selected_event_weights": [0.6, 0.4],
                }
            ],
            expected_names=["lhe_id=1001:scale up", "lhe_id=1002:scale down"],
            synthesized_name_events=1,
        )
        self.assertEqual(audit["status"], "inconsistent")
        self.assertEqual(audit["anonymous_variation_names"], ["optional_0", "optional_1"])
        self.assertTrue(
            any("synthesized positional" in message for message in audit["errors"])
        )
        self.assertTrue(
            any("despite manifest definitions" in message for message in audit["errors"])
        )

    def test_profile_identity_handles_background_na_and_detects_mismatch(self) -> None:
        def sample(source_checksum: str | None) -> dict:
            return {
                "category": "irreducible",
                "comparison_identity": {
                    "input_lhe_sha256": "input-checksum",
                    "source_lhe_sha256": source_checksum,
                    "production_group_id": "qq4l",
                    "nevents_requested": 100,
                    "source_events_read": 100,
                    # MadSpin is intentionally absent for continuum samples.
                    "seeds": {
                        "production": 11,
                        "madspin": None,
                        "herwig": 12,
                        "analysis": 13,
                    },
                },
            }

        identical = _comparison_identity(sample("source"), sample("source"))
        self.assertEqual(identical["status"], "verified_identical_source")
        self.assertNotIn("seed.madspin", identical["missing_fields"])

        mismatch = _comparison_identity(sample("source-a"), sample("source-b"))
        self.assertEqual(mismatch["status"], "mismatch")
        self.assertEqual(mismatch["mismatched_fields"], ["source_lhe_sha256"])

        missing = _comparison_identity(sample(None), sample(None))
        self.assertEqual(missing["status"], "unverified_missing_provenance")
        self.assertIn("source_lhe_sha256", missing["missing_fields"])


class FourLeptonReportOutputTest(unittest.TestCase):
    def _manifest(self, directory: Path) -> Path:
        manifest = {
            "schema_version": 1,
            "run_tag": "unit",
            "collider": {"sqrt_s_tev": 40},
            "luminosity_fb": 100,
            "detector_profile": "perfect",
            "samples": [
                {
                    "name": "signal_4e",
                    "label": "Higgs 4e",
                    "category": "signal",
                    "channel": "4e",
                    "perturbative_order": "LO",
                    "source_kind": "mg5",
                    "normalization": {
                        "native_cross_section_pb": 2.0,
                        "target_cross_section_pb": 1.0,
                        "weight_scale": 0.5,
                        "generated_sumw": 0.5,
                        "generated_sumw2": 1.25,
                    },
                    "output_root": "signal.root",
                    "production_group_id": "signal",
                },
                {
                    "name": "qq4l_4e",
                    "label": "qq to 4e",
                    "category": "irreducible",
                    "channel": "4e",
                    "perturbative_order": "LO",
                    "source_kind": "mg5",
                    "normalization": {
                        "native_cross_section_pb": 4.0,
                        "target_cross_section_pb": 4.0,
                        "weight_scale": 1.0,
                        "generated_sumw": 2.0,
                        "generated_sumw2": 2.0,
                    },
                    "output_root": "background.root",
                    "production_group_id": "qq4l",
                },
            ],
        }
        path = directory / "manifest.json"
        path.write_text(json.dumps(manifest))
        return path

    @staticmethod
    def _loader(path: Path) -> RootRecords:
        if path.name == "signal.root":
            sources = [
                {
                    "source_index": 1,
                    "generator_weight": 1.0,
                    "response_failure_probability": 0.0,
                    "optional_weight_names": ["scale_up"],
                    "optional_weights": [1.2],
                    "optional_selected_event_weights": [0.6],
                },
                {
                    "source_index": 2,
                    "generator_weight": -0.5,
                    "response_failure_probability": 0.0,
                    "optional_weight_names": ["scale_up"],
                    "optional_weights": [-0.6],
                    "optional_selected_event_weights": [-0.3],
                },
            ]
            rows = [
                {
                    "source_index": 1,
                    "hypothesis_probability": 1.0,
                    "cut_mask": SELECTED_MASK,
                    "channel": 1,
                    "pretrigger_event_weight": 0.5,
                    "event_weight": 0.5,
                    "m4l": 125.0,
                },
                {
                    "source_index": 2,
                    "hypothesis_probability": 1.0,
                    "cut_mask": SELECTED_MASK,
                    "channel": 1,
                    "pretrigger_event_weight": -0.25,
                    "event_weight": -0.25,
                    "m4l": 126.0,
                },
            ]
            return RootRecords(sources, rows)
        sources = [
            {
                "source_index": 1,
                "generator_weight": 1.0,
                "response_failure_probability": 0.0,
            },
            {
                "source_index": 2,
                "generator_weight": 1.0,
                "response_failure_probability": 0.0,
            },
        ]
        rows = [
            {
                "source_index": 1,
                "hypothesis_probability": 1.0,
                "cut_mask": SELECTED_MASK,
                "channel": 1,
                "pretrigger_event_weight": 1.0,
                "event_weight": 1.0,
                "m4l": 124.0,
            },
            {
                "source_index": 2,
                "hypothesis_probability": 1.0,
                "cut_mask": SELECTED_MASK,
                "channel": 1,
                "pretrigger_event_weight": 1.0,
                "event_weight": 1.0,
                "m4l": 140.0,
            },
        ]
        return RootRecords(sources, rows)

    @classmethod
    def _strict_loader(
        cls,
        path: Path,
        *,
        root_muon_efficiency: float = 1.0,
        summary_trigger_efficiency: float = 1.0,
    ) -> RootRecords:
        records = cls._loader(path)
        root_metadata = {
            "schema_version": 1,
            "response_profile": "perfect",
            "sample": "qq4l_4e",
            "category": "irreducible",
            "requested_channel": "4e",
            "seed": 101,
            "weight_scale": 1.0,
            "muon_resolution_scale": 1.0,
            "detector_parameter_source": "perfect_identity_response",
            "electron_efficiency": 1.0,
            "muon_efficiency": root_muon_efficiency,
            "trigger_4e_efficiency": 1.0,
            "trigger_other_efficiency": 1.0,
            "em_pileup_noise_enabled": False,
            "optional_weight_names_source": "not_available",
        }
        summary = {
            "schema_version": 1,
            "response_profile": "perfect",
            "sample": "qq4l_4e",
            "category": "irreducible",
            "requested_channel": "4e",
            "events_read": 2,
            "four_lepton_rows": 2,
            "seed": 101,
            "weight_scale": 1.0,
            "muon_resolution_scale": 1.0,
            "output_root": str(path.resolve()),
            "diagnostics": {
                "accepted_lepton_overflow_events": 0,
                "optional_weight_name_syntheses": 0,
                "optional_weight_name_mismatches": 0,
            },
            "detector": {
                "parameter_source": "perfect_identity_response",
                "electron_efficiency": 1.0,
                "muon_efficiency": 1.0,
                "trigger_4e_efficiency": summary_trigger_efficiency,
                "trigger_other_efficiency": 1.0,
                "em_pileup_noise_enabled": False,
            },
        }
        records.metadata = {
            **root_metadata,
            "analysis_metadata_present": True,
            "analysis_metadata": root_metadata,
            "summary": summary,
        }
        return records

    def test_normalization_json_html_and_signal_region(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            output = directory / "report"
            report = generate_report(
                [manifest],
                output,
                record_loader=self._loader,
            )

            signal = report["samples"][0]
            # native 2 pb * selected signed 0.25 / generated signed 0.5
            self.assertAlmostEqual(signal["selected"]["cross_section_pb"], 1.0)
            self.assertAlmostEqual(signal["signal_region"]["events"], 100_000.0)
            # source-grouped sumw2: (0.5)^2 + (-0.25)^2, scaled by (2/0.5)^2
            self.assertAlmostEqual(signal["selected"]["cross_section_sumw2_pb2"], 5.0)

            background = report["samples"][1]
            # Only one of two equal positive events is in the open 120--130 window.
            self.assertAlmostEqual(background["signal_region"]["events"], 200_000.0)
            self.assertAlmostEqual(
                report["summary"]["background_composition"][0][
                    "fraction_of_total_background"
                ],
                1.0,
            )
            self.assertTrue(report["audit"]["ok"])
            optional = signal["optional_weights"]
            self.assertEqual(optional["status"], "ok_raw_unnormalized")
            self.assertAlmostEqual(
                optional["variations"][0]["generated_raw"]["sumw"], 0.6
            )

            json_path = output / "report.json"
            html_path = output / "index.html"
            self.assertTrue(json_path.exists())
            self.assertTrue(html_path.exists())
            on_disk = json.loads(json_path.read_text())
            self.assertEqual(on_disk["schema_version"], 1)
            self.assertIn("native_cross_section_pb * sum(event_weight)", json_path.read_text())
            html_text = html_path.read_text()
            self.assertIn("Background composition", html_text)
            self.assertIn("Detector-response closure", html_text)
            self.assertIn("Optional-weight audit", html_text)
            self.assertIn("report.json", html_text)
            self.assertNotIn("{baseline", html_text)
            signal_spectrum = next(
                value
                for value in report["summary"]["mass_spectra"]
                if value["category"] == "signal"
            )
            self.assertAlmostEqual(
                sum(
                    value["cross_section_pb"]
                    for value in signal_spectrum["bins"]
                ),
                signal["selected"]["cross_section_pb"],
            )

    def test_strict_detector_metadata_and_summary_are_cross_checked(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            raw["samples"] = [raw["samples"][1]]
            raw["detector"] = {
                "profile": "perfect",
                "muon_resolution_scale": 1.0,
                "pileup_noise_enabled": False,
                "active_efficiencies": {
                    "electron": 1.0,
                    "muon": 1.0,
                    "trigger": {"all": 1.0},
                },
            }
            raw["samples"][0]["nevents_requested"] = 2
            raw["samples"][0]["seed"] = 101
            raw["samples"][0]["seeds"] = {"analysis": 101}
            raw["samples"][0]["analysis_summary"] = "background.summary.json"
            manifest.write_text(json.dumps(raw))

            valid = generate_report(
                [manifest],
                directory / "valid_report",
                record_loader=self._strict_loader,
            )
            self.assertTrue(valid["audit"]["ok"])
            detector_audit = valid["samples"][0]["detector_audit"]
            self.assertTrue(detector_audit["strict_contract"])
            self.assertTrue(detector_audit["ok"])
            self.assertTrue(
                all(
                    check["root_status"] == "match"
                    and check["summary_status"] == "match"
                    for check in detector_audit["checks"]
                )
            )

            mismatched = generate_report(
                [manifest],
                directory / "mismatched_report",
                record_loader=lambda path: self._strict_loader(
                    path,
                    root_muon_efficiency=0.85,
                    summary_trigger_efficiency=0.98,
                ),
            )
            self.assertFalse(mismatched["audit"]["ok"])
            errors = mismatched["samples"][0]["detector_audit"]["errors"]
            self.assertTrue(
                any(
                    "muon_efficiency" in message
                    and "ROOT AnalysisMetadata" in message
                    for message in errors
                )
            )
            self.assertTrue(
                any(
                    "trigger_4e_efficiency" in message
                    and "analyzer summary" in message
                    for message in errors
                )
            )

    def test_accepted_lepton_overflow_is_an_error_and_forces_exclusion(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            raw["samples"] = [raw["samples"][1]]
            manifest.write_text(json.dumps(raw))

            def overflow_loader(path: Path) -> RootRecords:
                records = self._loader(path)
                records.source_events[0]["accepted_lepton_count"] = 11
                records.source_events[1]["accepted_lepton_count"] = 4
                records.metadata["summary"] = {
                    "diagnostics": {"accepted_lepton_overflow_events": 1}
                }
                return records

            report = generate_report(
                [manifest],
                directory / "overflow_report",
                record_loader=overflow_loader,
                allow_unvalidated_bias=True,
            )
            sample = report["samples"][0]
            self.assertFalse(sample["included_in_physics_totals"])
            self.assertEqual(sample["exclusion_reason"], "accepted_lepton_overflow")
            self.assertTrue(sample["detector_audit"]["accepted_lepton_overflow"])
            self.assertEqual(
                report["summary"]["signal_region"]["background_events"], 0
            )
            self.assertTrue(
                any(
                    "accepted-lepton enumeration overflow" in message
                    for message in report["audit"]["errors"]
                )
            )

    def test_zero_signed_denominator_is_reported_not_replaced_by_entries(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            raw["samples"] = [raw["samples"][0]]
            raw["samples"][0]["normalization"]["generated_sumw"] = 0.0
            raw["samples"][0]["normalization"]["generated_sumw2"] = 2.0
            manifest.write_text(json.dumps(raw))

            def loader(_path: Path) -> RootRecords:
                sources = [
                    {"source_index": 1, "generator_weight": 1.0, "response_failure_probability": 0.0},
                    {"source_index": 2, "generator_weight": -1.0, "response_failure_probability": 0.0},
                ]
                rows = [
                    {
                        "source_index": 1,
                        "hypothesis_probability": 1.0,
                        "cut_mask": SELECTED_MASK,
                        "event_weight": 0.5,
                        "m4l": 125.0,
                    },
                    {
                        "source_index": 2,
                        "hypothesis_probability": 1.0,
                        "cut_mask": SELECTED_MASK,
                        "event_weight": -0.5,
                        "m4l": 125.0,
                    },
                ]
                return RootRecords(sources, rows)

            report = generate_report([manifest], directory / "report", record_loader=loader)
            sample = report["samples"][0]
            self.assertEqual(
                sample["normalization"]["status"], "unavailable_zero_signed_denominator"
            )
            self.assertIsNone(sample["selected"]["events"])
            self.assertTrue(
                any("zero" in message for message in report["audit"]["errors"])
            )

    def test_unvalidated_decay_bias_is_audited_and_excluded_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            raw["samples"] = [raw["samples"][1]]
            raw["samples"][0]["decay_bias"] = {
                "enabled": True,
                "scheme": "herwig_branching_ratio_reweighter",
                "closure_validated": False,
            }
            raw["samples"][0]["normalization"]["denominator_sumw"] = 2.0
            raw["samples"][0]["normalization"]["denominator_sumw2"] = 2.0
            manifest.write_text(json.dumps(raw))

            report = generate_report([manifest], directory / "report", record_loader=self._loader)
            sample = report["samples"][0]
            self.assertFalse(sample["included_in_physics_totals"])
            self.assertEqual(sample["exclusion_reason"], "unvalidated_decay_bias")
            self.assertEqual(report["summary"]["signal_region"]["background_events"], 0)
            self.assertEqual(
                report["summary"]["excluded_samples"][0]["name"], "qq4l_4e"
            )
            self.assertTrue(
                any("decay-bias closure" in message for message in report["audit"]["errors"])
            )

            included = generate_report(
                [manifest],
                directory / "included_report",
                record_loader=self._loader,
                allow_unvalidated_bias=True,
            )
            self.assertTrue(included["samples"][0]["included_in_physics_totals"])
            self.assertAlmostEqual(
                included["summary"]["signal_region"]["background_events"], 200_000.0
            )

    def test_heavy_flavour_bias_control_comparison_is_diagnostic_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            control = raw["samples"][1]
            biased = json.loads(json.dumps(control))
            biased["name"] = "qq4l_4e_hfbiased"
            biased["label"] = "qq to 4e (HF biased)"
            biased["output_root"] = "background_hfbiased.root"
            biased["decay_bias"] = {
                "enabled": True,
                "scheme": "restricted_direct_b_semileptonic_importance_sampling",
                "bias_of": control["name"],
                "full_tail_support": False,
                "closure_validated": False,
            }
            biased["normalization"]["denominator_sumw"] = 2.0
            biased["normalization"][
                "denominator_kind"
            ] = "predecay_lhe_handler_weight_sum"
            raw["samples"] = [control, biased]
            manifest.write_text(json.dumps(raw))

            report = generate_report(
                [manifest],
                directory / "closure_report",
                record_loader=self._loader,
            )
            diagnostic = report["summary"][
                "heavy_flavour_closure_diagnostics"
            ][0]
            self.assertEqual(
                diagnostic["status"], "diagnostic_available_not_validation"
            )
            self.assertEqual(diagnostic["control_sample"], "qq4l_4e")
            self.assertFalse(diagnostic["closure_validation_claimed_by_report"])
            self.assertFalse(diagnostic["full_tail_support"])
            self.assertAlmostEqual(
                diagnostic["selected"]["ratio_biased_over_control"], 1.0
            )
            biased_report = next(
                sample
                for sample in report["samples"]
                if sample["name"] == "qq4l_4e_hfbiased"
            )
            self.assertFalse(biased_report["included_in_physics_totals"])
            self.assertFalse(
                biased_report["decay_bias"]["closure_validated"]
            )
            self.assertIn(
                "transparent weighted comparisons only",
                (directory / "closure_report/index.html").read_text(),
            )

    def test_forced_decay_requires_explicit_predecay_denominator(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            manifest = self._manifest(directory)
            raw = json.loads(manifest.read_text())
            raw["samples"] = [raw["samples"][1]]
            raw["samples"][0]["prompt_decay_profile"] = "w_emu"
            manifest.write_text(json.dumps(raw))

            report = generate_report([manifest], directory / "report", record_loader=self._loader)
            sample = report["samples"][0]
            self.assertEqual(
                sample["normalization"]["status"],
                "unavailable_missing_predecay_denominator",
            )
            self.assertIsNone(sample["signal_region"]["events"])
            self.assertTrue(
                any(
                    "would cancel the decay branching weight" in message
                    for message in report["audit"]["errors"]
                )
            )

            raw["samples"][0]["normalization"]["denominator_sumw"] = 2.0
            raw["samples"][0]["normalization"]["denominator_sumw2"] = 2.0
            raw["samples"][0]["normalization"][
                "denominator_kind"
            ] = "predecay_lhe_handler_weight_sum"
            manifest.write_text(json.dumps(raw))
            valid = generate_report(
                [manifest], directory / "valid_report", record_loader=self._loader
            )
            self.assertEqual(valid["samples"][0]["normalization"]["status"], "ok")
            self.assertEqual(
                valid["samples"][0]["normalization"]["denominator_source"],
                "manifest_predecay_sumw",
            )
            self.assertFalse(
                any(
                    "unrecognized normalization denominator kind" in message
                    for message in valid["audit"]["warnings"]
                )
            )

    def test_matching_perfect_and_ssc_samples_are_compared(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)

            def write_manifest(profile: str) -> Path:
                path = directory / f"manifest_{profile}.json"
                path.write_text(
                    json.dumps(
                        {
                            "schema_version": 1,
                            "run_tag": "compare",
                            "collider": {"sqrt_s_tev": 40},
                            "detector_profile": profile,
                            "samples": [
                                {
                                    "name": "signal_4e",
                                    "label": "Higgs 4e",
                                    "category": "signal",
                                    "channel": "4e",
                                    "production_group_id": "signal_stable",
                                    "nevents_requested": 2,
                                    "input_lhe_sha256": "decayed-checksum",
                                    "source_lhe_sha256": "stable-checksum",
                                    "seeds": {
                                        "production": 11,
                                        "madspin": 12,
                                        "herwig": 13,
                                        "analysis": 14,
                                    },
                                    "normalization": {
                                        "native_cross_section_pb": 1.0,
                                        "target_cross_section_pb": 1.0,
                                        "weight_scale": 1.0,
                                        "generated_sumw": 2.0,
                                        "generated_sumw2": 2.0,
                                    },
                                    "output_root": f"signal_{profile}.root",
                                }
                            ],
                        }
                    )
                )
                return path

            def loader(path: Path) -> RootRecords:
                sources = [
                    {
                        "source_index": 1,
                        "generator_weight": 1.0,
                        "response_failure_probability": 0.0,
                    },
                    {
                        "source_index": 2,
                        "generator_weight": 1.0,
                        "response_failure_probability": 0.0,
                    },
                ]
                if "ssc" in path.name:
                    masses, weights = (123.0, 127.0), (0.8, 0.8)
                else:
                    masses, weights = (124.0, 126.0), (1.0, 1.0)
                rows = [
                    {
                        "source_index": index,
                        "hypothesis_probability": 1.0,
                        "cut_mask": SELECTED_MASK,
                        "channel": 1,
                        "pretrigger_event_weight": weight,
                        "event_weight": weight,
                        "m4l": mass,
                    }
                    for index, (mass, weight) in enumerate(
                        zip(masses, weights), start=1
                    )
                ]
                return RootRecords(sources, rows)

            perfect_manifest = write_manifest("perfect")
            ssc_manifest = write_manifest("ssc")
            report = generate_report(
                [perfect_manifest, ssc_manifest],
                directory / "report",
                record_loader=loader,
            )
            comparison = report["summary"]["perfect_ssc_comparisons"][0]
            self.assertEqual(
                comparison["source_identity"]["status"],
                "verified_identical_source",
            )
            self.assertAlmostEqual(
                comparison["selected_events"]["ratio_ssc_over_perfect"], 0.8
            )
            self.assertAlmostEqual(comparison["mass_sigma_gev"]["perfect"], 1.0)
            self.assertAlmostEqual(comparison["mass_sigma_gev"]["ssc"], 2.0)
            self.assertAlmostEqual(
                comparison["mass_sigma_gev"]["difference_ssc_minus_perfect"], 1.0
            )
            # The two response profiles describe the same hard samples and
            # must never be added into one nominal physics yield.
            self.assertIsNone(
                report["summary"]["signal_region"]["signal_events"]
            )
            self.assertAlmostEqual(
                report["summary"]["physics_by_profile"]["perfect"][
                    "signal_region"
                ]["signal_events"],
                100_000.0,
            )
            self.assertAlmostEqual(
                report["summary"]["physics_by_profile"]["ssc"][
                    "signal_region"
                ]["signal_events"],
                80_000.0,
            )
            html_text = (directory / "report/index.html").read_text()
            self.assertIn("Perfect versus SSC", html_text)
            self.assertEqual(html_text.count("Single-profile background yield"), 1)

            # Duplicate sample/profile entries are ambiguous and must not be
            # silently resolved according to manifest input order.
            perfect_copy = directory / "manifest_perfect_copy.json"
            perfect_copy.write_text(perfect_manifest.read_text())
            duplicate = generate_report(
                [perfect_manifest, perfect_copy, ssc_manifest],
                directory / "duplicate_report",
                record_loader=loader,
            )
            self.assertTrue(
                any(
                    "duplicate sample/profile identity" in message
                    for message in duplicate["audit"]["errors"]
                )
            )
            self.assertEqual(duplicate["summary"]["perfect_ssc_comparisons"], [])

            # A checksum disagreement invalidates a detector-only comparison
            # and is promoted from comparison metadata to the global audit.
            ssc_raw = json.loads(ssc_manifest.read_text())
            ssc_raw["samples"][0]["source_lhe_sha256"] = "different-stable-checksum"
            ssc_manifest.write_text(json.dumps(ssc_raw))
            mismatched = generate_report(
                [perfect_manifest, ssc_manifest],
                directory / "mismatch_report",
                record_loader=loader,
            )
            identity = mismatched["summary"]["perfect_ssc_comparisons"][0][
                "source_identity"
            ]
            self.assertEqual(identity["status"], "mismatch")
            self.assertFalse(identity["detector_only_interpretation_valid"])
            self.assertTrue(
                any(
                    "perfect/SSC source identity differs" in message
                    for message in mismatched["audit"]["errors"]
                )
            )


if __name__ == "__main__":
    unittest.main()
