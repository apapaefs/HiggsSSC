import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from hfourlepton import evaluate_reducible_backgrounds as evaluator


def manifests():
    base = dict(status="complete", collider={"sqrt_s_tev": 40.}, order="ho",
        detector_profile="ssc", selection_profile="ssc", detector={"electron_efficiency": .9},
        cut_stages=[{"name": "selected", "bit": 9}],
        samples=[dict(name="signal", selection_version="ssc_v1", input_lhe_sha256="old-input",
                      source_lhe_sha256="old-source", output_root_sha256="old-root")])
    extra = {**copy.deepcopy(base), "order": "lo", "run_tag": "control", "samples": []}
    extra["samples"] = [dict(name="bkg_ttbar", status="complete", category="reducible",
        perturbative_order="LO", selection_version="ssc_v1", process="p p > t t~",
        decay_bias={"enabled": False}, input_lhe_sha256="new-input", source_lhe_sha256="new-source",
        output_root_sha256="new-root", seeds={"production": 23},
        normalization={"denominator_event_count": 100, "native_cross_section_pb": 10.})]
    return base, extra


def row(weight, shard="a", source=0, sample="bkg_ttbar", channel="4e", mass=125., **features):
    values = {name: 1. for name in evaluator.analysis.VARIABLES}
    values.update(m4l=mass, **features)
    return dict(sample=sample, shard=shard, group=f"{shard}/{sample}", source=source,
                channel=channel, weight=weight, **values)


def meta(shard="a", exposure=100., sample="bkg_ttbar", **kwargs):
    return dict(sample=sample, shard=shard, generated_exposure=exposure, **kwargs)


def weighted_pair():
    base, control = manifests()
    sample = control["samples"][0]
    sample.update(heavy_flavour_decay_backend="herwig_native_table",
                  heavy_flavour_decay_table_sha256="decay-table", heavy_flavour_decay_table="/example/HerwigBDecays.in")
    weighted = copy.deepcopy(control)
    weighted["run_tag"] = "weighted"
    sample = weighted["samples"][0]
    sample["name"] += evaluator.BIAS_SUFFIX
    sample["decay_bias"] = dict(enabled=True, bias_of="bkg_ttbar", full_tail_support=True,
        closure_validated=False, scheme="full_support_direct_b_semileptonic_importance_sampling")
    sample["branching_ratio_reweighter"] = True
    sample["normalization"]["denominator_kind"] = "predecay_lhe_handler_weight_sum"
    sample["seeds"]["production"] += 1
    for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
        sample[field] += "-weighted"
    return base, control, weighted


class PredictOnly:
    def __init__(self):
        self.calls = []

    def predict_proba(self, matrix):
        self.calls.append(matrix.copy())
        score = matrix[:, 0] / 100.
        return np.column_stack((1. - score, score))

    def fit(self, *args, **kwargs):
        raise AssertionError("the frozen model must never be fitted")


class ReducibleEvaluationTests(unittest.TestCase):
    def test_weighted_inputs_require_explicit_proposal_and_matching_decay_backend(self):
        base, control, weighted = weighted_pair()
        with self.assertRaises(ValueError):
            evaluator.validate_manifests(base, [control, weighted])
        evaluator.validate_manifests(base, [control, weighted], study_decay_bias=True)
        for case in ("support", "denominator", "backend"):
            changed = copy.deepcopy(weighted)
            sample = changed["samples"][0]
            if case == "support": sample["decay_bias"]["full_tail_support"] = False
            elif case == "denominator": sample["normalization"]["denominator_kind"] = "SourceEvents.generator_weight"
            else: sample["heavy_flavour_decay_table_sha256"] = "different-table"
            with self.subTest(case=case), self.assertRaises(ValueError):
                evaluator.validate_manifests(base, [control, changed], study_decay_bias=True)
        control["samples"][0]["heavy_flavour_decay_backend"] = "evtgen"
        with self.assertRaisesRegex(ValueError, "physics differs"):
            evaluator.validate_manifests(base, [control, weighted], study_decay_bias=True)

    def test_weighted_and_unbiased_can_share_only_disjoint_parent_ranges(self):
        base, control, weighted = weighted_pair()
        def proof(manifest, start):
            sample = manifest["samples"][0]
            return {sample["name"]: dict(parent_sha256="parent", lhe_sha256=sample["source_lhe_sha256"],
                run_tag=manifest["run_tag"], start=start, stop=start+100, events=100,
                parent_cross_section_pb=10., production_seed=321)}
        proofs = [proof(control, 0), proof(weighted, 100)]
        evaluator.validate_manifests(base, [control, weighted], proofs, study_decay_bias=True)
        proofs[1]["bkg_ttbar_hfbiased"]["start"] = 0
        with self.assertRaisesRegex(ValueError, "overlap"):
            evaluator.validate_manifests(base, [control, weighted], proofs, study_decay_bias=True)

    def test_weighted_yields_never_enter_unbiased_totals_or_nominal_promotion(self):
        rows = [row(2., mZ1=90.), row(8., shard="weighted", sample="bkg_ttbar_hfbiased", mZ1=90.)]
        metadata = [meta(), meta("weighted", sample="bkg_ttbar_hfbiased")]
        data, receipt = evaluator.combine_shards(rows, metadata)
        card = dict(features=["mZ1"])
        result, *_ = evaluator.evaluate_data(data, receipt, card, card["features"], .5, PredictOnly(), study_decay_bias=True)
        self.assertEqual(result["cuts"]["total"]["sum_weight"], 2.)
        self.assertEqual(result["weighted_diagnostics"]["cuts"]["total"]["sum_weight"], 8.)
        self.assertFalse(result["weighted_diagnostics"]["included_in_nominal_totals"])
        closure = evaluator.closure_diagnostics(result, {}, 100.)["bkg_ttbar_hfbiased"]
        self.assertEqual(closure["cuts"]["ratio_biased_over_control"], 4.)
        self.assertFalse(closure["cuts"]["statistically_supported"])
        self.assertIsNone(closure["cuts"]["difference_over_combined_stat_uncertainty"])
        self.assertFalse(closure["closure_validated"])

    def test_zero_zero_is_inconclusive_and_supported_closure_stays_diagnostic(self):
        empty = dict(sum_weight=0., variance=0., effective_count=0., selected_source_events=0)
        zero = evaluator.closure_metric(empty, empty, 100.)
        self.assertIsNone(zero["ratio_biased_over_control"])
        self.assertIsNone(zero["difference_over_combined_stat_uncertainty"])
        self.assertIsNone(zero["difference_95pct_gaussian_interval_pb"])
        self.assertFalse(zero["statistically_supported"])
        supported = dict(sum_weight=100., variance=100., effective_count=100., selected_source_events=100)
        result = evaluator.closure_metric(supported, supported, 100.)
        self.assertTrue(result["statistically_supported"])
        self.assertFalse(result["closure_validated"])
        lo, hi = result["difference_95pct_gaussian_interval_pb"]
        self.assertLess(lo, 0.)
        self.assertAlmostEqual(lo, -hi)

    def test_sidebands_use_predecay_scale_and_separate_strata(self):
        metadata = [meta(manifest="control", native_cross_section_pb=10., combination_coefficient=1.,
                         normalization={"scale_pb_per_weight": .1}),
                    meta("weighted", sample="bkg_ttbar_hfbiased", manifest="weighted", native_cross_section_pb=10.,
                         combination_coefficient=1., normalization={"scale_pb_per_weight": .1})]
        def reader(path):
            return ("bkg_ttbar" + ("_hfbiased" if path == "weighted" else ""), 100, 10.,
                    {"open/3P1F/all": [2., 2., 2] if path == "control" else [1., .5, 2]})
        result = evaluator.sideband_diagnostics(metadata, 1., reader=reader)
        self.assertEqual(result["strata"]["bkg_ttbar"]["open/3P1F/all"]["rate_fb"], 200.)
        weighted = result["strata"]["bkg_ttbar_hfbiased"]["open/3P1F/all"]
        self.assertEqual(weighted["rate_fb"], 100.)
        self.assertEqual(weighted["variance"], 5000.)
        self.assertEqual(weighted["effective_count"], 2.)
        self.assertEqual(result["weighted_control_closure"]["bkg_ttbar_hfbiased"]["open/3P1F/all"]["ratio_biased_over_control"], .5)

    def test_lo_controls_can_differ_in_hard_physics_but_not_detector(self):
        base, extra = manifests()
        evaluator.validate_manifests(base, [extra])
        for field in evaluator.DETECTOR_FIELDS:
            changed = copy.deepcopy(extra)
            changed[field] = "changed"
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "response/selection"):
                evaluator.validate_manifests(base, [changed])

    def test_rejects_overlap_missing_provenance_bias_and_unexpected_inventory(self):
        for case in ("old-source", "missing-hash", "bias", "incomplete", "ttz", "old-name"):
            base, extra = manifests()
            sample = extra["samples"][0]
            if case == "old-source": sample["source_lhe_sha256"] = "old-source"
            elif case == "missing-hash": del sample["input_lhe_sha256"]
            elif case == "bias": sample["decay_bias"]["enabled"] = True
            elif case == "incomplete": sample["status"] = "planned"
            elif case == "ttz": sample["name"] = "bkg_ttz"
            else: base["samples"][0]["name"] = sample["name"]
            with self.subTest(case=case), self.assertRaises(ValueError):
                evaluator.validate_manifests(base, [extra])

    def test_repeated_shards_require_independence_and_same_process_settings(self):
        base, extra = manifests()
        other = copy.deepcopy(extra)
        for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
            other["samples"][0][field] += "-2"
        with self.assertRaisesRegex(ValueError, "production seed"):
            evaluator.validate_manifests(base, [extra, other])
        other["samples"][0]["seeds"]["production"] += 1
        evaluator.validate_manifests(base, [extra, other])
        other["samples"][0]["process"] = "changed"
        with self.assertRaisesRegex(ValueError, "shards change"):
            evaluator.validate_manifests(base, [extra, other])

    def test_common_parent_seed_requires_disjoint_verified_ranges(self):
        base, extra = manifests()
        other = copy.deepcopy(extra)
        other["run_tag"] = "control2"
        for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
            other["samples"][0][field] += "-2"
        def proof(manifest, start):
            return {"bkg_ttbar": dict(parent_sha256="parent", lhe_sha256=manifest["samples"][0]["source_lhe_sha256"],
                run_tag=manifest["run_tag"], start=start, stop=start+100, events=100, parent_cross_section_pb=10.)}
        proofs = [proof(extra, 0), proof(other, 100)]
        evaluator.validate_manifests(base, [extra, other], proofs)
        proofs[1]["bkg_ttbar"]["start"] = 90
        with self.assertRaisesRegex(ValueError, "overlap"):
            evaluator.validate_manifests(base, [extra, other], proofs)

    def test_generation_domains_cannot_mix_and_partition_requires_index(self):
        base, extra = manifests()
        other = copy.deepcopy(extra)
        other["run_tag"] = "control2"
        other["samples"][0]["seeds"]["production"] = 24
        for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
            other["samples"][0][field] += "-2"
        policy = {"name": "inclusive_heavy_flavour_v1", "cut_overrides": {"ptj": 0.}}
        other["reducible_control_partition"] = {"generation_policy": policy}
        proof = {"bkg_ttbar": dict(parent_sha256="parent", lhe_sha256="new-source-2", run_tag="control2",
                 start=0, stop=100, events=100, parent_cross_section_pb=10., generation_policy=policy)}
        with self.assertRaisesRegex(ValueError, "shard-index"):
            evaluator.validate_manifests(base, [other])
        evaluator.validate_manifests(base, [other], [proof])
        with self.assertRaisesRegex(ValueError, "generation_policy"):
            evaluator.validate_manifests(base, [extra, other], [{}, proof])
        proof["bkg_ttbar"]["generation_policy"] = evaluator.LEGACY_GENERATION_POLICY
        with self.assertRaisesRegex(ValueError, "generation policies differ"):
            evaluator.validate_manifests(base, [other], [proof])

    def test_independent_parent_runs_use_generation_seeds_not_worker_seeds(self):
        base, extra = manifests()
        other = copy.deepcopy(extra)
        other["run_tag"] = "control2"
        other["samples"][0]["seeds"]["production"] = 999
        for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256"):
            other["samples"][0][field] += "-2"
        def proof(manifest, parent):
            return {"bkg_ttbar": dict(parent_sha256=parent, lhe_sha256=manifest["samples"][0]["source_lhe_sha256"],
                run_tag=manifest["run_tag"], start=0, stop=100, events=100,
                parent_cross_section_pb=10., production_seed=123, actual_generation_cuts={"ptj": 0.})}
        proofs = [proof(extra, "parent-a"), proof(other, "parent-b")]
        with self.assertRaisesRegex(ValueError, "production seed"):
            evaluator.validate_manifests(base, [extra, other], proofs)
        proofs[1]["bkg_ttbar"]["production_seed"] = 124
        evaluator.validate_manifests(base, [extra, other], proofs)
        proofs[1]["bkg_ttbar"]["actual_generation_cuts"]["ptj"] = 20.
        with self.assertRaisesRegex(ValueError, "actual generation cuts"):
            evaluator.validate_manifests(base, [extra, other], proofs)

    def test_index_checks_content_hashes_and_half_open_intervals(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            parent, shard, manifest = root/"parent.lhe", root/"shard.lhe", root/"manifest.json"
            parent.write_text("parent"); shard.write_text("shard")
            item = dict(index=0, start=0, stop=100, events=100, lhe=str(shard),
                        sha256=evaluator.analysis.digest(shard), run_tag="control", manifest=str(manifest))
            payload = dict(samples={"bkg_ttbar": dict(source_lhe=str(parent),
                source_sha256=evaluator.analysis.digest(parent), source_events=100, cross_section_pb=10., shards=[item])})
            index = root/"index.json"
            index.write_text(json.dumps(payload))
            proof = evaluator.shard_proofs([index], [manifest])
            self.assertEqual(proof[str(manifest.resolve()), "bkg_ttbar"]["events"], 100)
            shard.write_text("changed")
            with self.assertRaises(ValueError):
                evaluator.shard_proofs([index], [manifest])

    def test_independent_rates_average_and_distinct_processes_add(self):
        rows = [row(3., source=0), row(-1., source=1), row(1., shard="b", source=0),
                row(2., shard="b", source=0), row(-1., shard="b", source=1),
                row(5., shard="c", sample="bkg_zbb_ee")]
        data, receipt = evaluator.combine_shards(rows, [meta(), meta("b", 300.), meta("c", sample="bkg_zbb_ee")])
        result = evaluator.selected_yields(data, np.ones(len(rows), bool), {"bkg_ttbar", "bkg_zbb_ee"})
        self.assertEqual(result["by_process"]["bkg_ttbar"]["sum_weight"], 2.)
        self.assertEqual(result["total"]["sum_weight"], 7.)
        self.assertAlmostEqual(result["total"]["variance"], (.25**2 + .75**2)*10. + 25.)
        self.assertEqual(result["total"]["selected_source_events"], 5)
        self.assertEqual([r["combination_coefficient"] for r in receipt], [.25, .75, 1.])
        self.assertEqual(rows[0]["weight"], 3.)

    def test_parent_slices_pool_signed_denominators_then_parent_exposure(self):
        common = {"parent_sha256": "parent"}
        rows = [row(6.), row(2., shard="b"), row(3., shard="c")]
        metadata = [meta(parent_slice=common, normalization={"denominator_sumw": 1.}),
                    meta("b", parent_slice=common, normalization={"denominator_sumw": 3.}), meta("c", 200.)]
        data, receipt = evaluator.combine_shards(rows, metadata)
        np.testing.assert_allclose([r["combination_coefficient"] for r in receipt], [.125, .375, .5])
        self.assertEqual(data["weight"].sum(), 3.)

    def test_channel_totals_keep_correlated_hypotheses_in_one_source(self):
        data, _ = evaluator.combine_shards([row(2., channel="4e"), row(3., channel="4mu")], [meta()])
        result = evaluator.selected_yields(data, np.ones(2, bool), {"bkg_ttbar"})
        self.assertEqual(result["total"]["variance"], 25.)
        self.assertEqual(result["by_channel"]["4e"]["variance"], 4.)
        self.assertEqual(result["by_channel"]["4mu"]["variance"], 9.)
        self.assertEqual(result["total"]["selected_source_events"], 1)

    def test_frozen_features_threshold_and_strict_mass_windows_without_projection(self):
        card = dict(features=["mZ2", "mZ1"], preselection={"m4l": {"min": 105., "max": 140.}},
                    cuts={"m4l": {"min": 120., "max": 130.}})
        rows = [row(2., mZ2=70., mZ1=90.), row(-1., source=1, mZ2=60., mZ1=91.),
                row(4., source=2, mZ2=40., mZ1=92.), row(8., source=3, mass=120., mZ2=90.),
                row(16., source=4, mass=140., mZ2=90.)]
        data, receipt = evaluator.combine_shards(rows, [meta()])
        model = PredictOnly()
        result, _, eligible, baseline = evaluator.evaluate_data(data, receipt, card, card["features"], .5, model)
        self.assertEqual(len(model.calls), 1)
        np.testing.assert_array_equal(model.calls[0][:3], [[70.,90.], [60.,91.], [40.,92.]])
        self.assertEqual(result["frozen_classifier"]["total"]["sum_weight"], 1.)
        self.assertEqual(result["frozen_classifier"]["total"]["variance"], 5.)
        self.assertEqual(result["cuts"]["total"]["sum_weight"], 5.)
        np.testing.assert_array_equal(eligible, [True, True, True, True, False])
        np.testing.assert_array_equal(baseline, [True, True, True, False, False])
        with self.assertRaisesRegex(ValueError, "feature order"):
            evaluator.evaluate_data(data, receipt, card, list(reversed(card["features"])), .5, model)

    def test_empty_selection_reports_zero_counts_without_predicting(self):
        data, receipt = evaluator.combine_shards([], [meta()])
        model = PredictOnly()
        card = dict(features=["mZ1"])
        result, _, _, _ = evaluator.evaluate_data(data, receipt, card, card["features"], .5, model)
        self.assertEqual(model.calls, [])
        self.assertEqual(result["frozen_classifier"]["total"]["selected_source_events"], 0)
        self.assertEqual(result["cuts"]["by_process"]["bkg_ttbar"]["sum_weight"], 0.)

    def test_row_conversion_uses_all_source_provenance_and_signed_scale(self):
        sample = SimpleNamespace(name="bkg_ttbar", raw={}, production_group_id="ttbar", channel="inclusive")
        stage = evaluator.analysis.report.CutStage("selected", "Selected", 9, 1023)
        record = dict(source_index=2, channel=3, cut_mask=1023, event_weight=-.5,
            m4l=125., mZ1=90., mZ2=25., pt4l=30., min_delta_r=.2, y4l=-1., lepton_pt=[10.,40.,20.,30.])
        records = SimpleNamespace(source_events=[{"source_index": 1}, {"source_index": 2}], four_lepton=[record])
        rows = evaluator.rows_from_records(sample, records, [stage], .002, 100., "shard")
        self.assertEqual(rows[0]["weight"], -100.)
        self.assertEqual(rows[0]["channel"], "2e2mu")
        self.assertEqual([rows[0][f"pt{i}"] for i in range(1,5)], [40.,30.,20.,10.])
        records.source_events = [{"source_index": 1}]
        with self.assertRaisesRegex(ValueError, "no SourceEvents"):
            evaluator.rows_from_records(sample, records, [stage], .002, 100., "shard")

    def test_cutflow_averages_rates_squares_coefficients_and_adds_sources(self):
        flow = [dict(name="isolation", label="Isolation", bit=3, yield_events=4., yield_sumw2=10.,
                     independent_source_events=2)]
        result = evaluator.pooled_cutflow([meta(cutflow=flow, combination_coefficient=.25),
                                          meta("b", cutflow=flow, combination_coefficient=.75)])
        row = result["bkg_ttbar"][0]
        self.assertEqual(row["yield_events"], 4.)
        self.assertEqual(row["yield_sumw2"], 6.25)
        self.assertEqual(row["independent_source_events"], 4)

    def test_zero_selection_run_writes_auditable_results_without_modifying_reference(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            frozen_dir, output = root/"frozen", root/"result"
            frozen_dir.mkdir()
            for filename in ("completion.json", "summary.json", "model.json"):
                (frozen_dir/filename).write_text("{}")
            before = {path.name: path.read_bytes() for path in frozen_dir.iterdir()}
            base, extra = manifests()
            manifest = root/"manifest.json"
            manifest.write_text(json.dumps(extra))
            card = dict(features=["mZ1"], preselection={"m4l": {"min":105., "max":140.}},
                        cuts={"m4l": {"min":120., "max":130.}})
            summary = dict(luminosity_fb=100., card=card, features=card["features"],
                           threshold={"threshold": .5})
            model = PredictOnly()
            model.load_model = mock.Mock()
            metadata = [meta(cutflow=[])]
            with mock.patch.object(evaluator.frozen, "frozen_analysis", return_value=(summary, base)), \
                    mock.patch.object(evaluator, "load_shard", return_value=([], metadata)), \
                    mock.patch.dict("sys.modules", {"xgboost": SimpleNamespace(XGBClassifier=lambda: model)}):
                evaluator.run([manifest], frozen_dir, output)
            result = json.loads((output/"summary.json").read_text())
            self.assertEqual(result["frozen_classifier"]["total"]["sum_weight"], 0.)
            self.assertEqual(result["frozen_classifier"]["total"]["selected_source_events"], 0)
            self.assertEqual(model.calls, [])
            self.assertEqual(len((output/"scores.csv").read_text().splitlines()), 1)
            self.assertEqual(result["missing_planned_processes"], ["bkg_zbb_ee", "bkg_zbb_mumu"])
            self.assertIn("Zero selected events", result["limitations"][1])
            self.assertTrue((output/"cutflow.csv").is_file())
            self.assertEqual({path.name: path.read_bytes() for path in frozen_dir.iterdir()}, before)
            completion = json.loads((output/"completion.json").read_text())
            for name, checksum in completion["files"].items():
                self.assertEqual(evaluator.analysis.digest(output/name), checksum)


if __name__ == "__main__":
    unittest.main()
