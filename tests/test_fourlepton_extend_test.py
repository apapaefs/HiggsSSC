import copy
import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from hfourlepton import extend_xgboost_test as extension
from hfourlepton import analyze_fourlepton as analysis


def fixture():
    base = dict(status="complete", collider={"sqrt_s_tev":40.}, detector_profile="ssc",
                selection_profile="ssc", order="ho", sampling={"continuum_m4l_gev":[90.,180.]}, samples=[])
    for i, name in enumerate(sorted(extension.QQ)):
        base["samples"].append(dict(name=name, status="complete", input_lhe_sha256=name+"current",
            source_lhe_sha256=name+"current", output_root_sha256=name+"currentroot",
            seeds={"production":10+i}, normalization={"denominator_event_count":1000000}))
    extra = copy.deepcopy(base)
    for i, row in enumerate(extra["samples"]):
        row.update(input_lhe_sha256=row["name"]+"independent", source_lhe_sha256=row["name"]+"independent",
                   output_root_sha256=row["name"]+"independentroot", seeds={"production":100+i},
                   normalization={"denominator_event_count":200000})
    summary = dict(samples=[dict(name=s["name"], group=s["name"], root_sha256=s["output_root_sha256"],
                                normalization=s["normalization"]) for s in base["samples"]],
                   split_fractions={s["name"]:{"test":.2} for s in base["samples"]})
    return base, extra, summary


class ExtendTestTests(unittest.TestCase):
    def test_independent_matching_population_is_subset_without_edits(self):
        base, extra, summary = fixture()
        before = copy.deepcopy((base, extra, summary))
        subset = extension.compatible_subset(base, extra, summary)
        self.assertEqual(set(s["name"] for s in subset["samples"]), extension.QQ)
        self.assertEqual((base, extra, summary), before)
        coefficients = extension.exposure_coefficients(summary, subset)
        for value in coefficients.values():
            self.assertEqual(value["original_fraction"], .5)
            self.assertEqual(value["additional_fraction"], .5)

    def test_rejects_training_overlap_changed_physics_and_window(self):
        for field in ("input_lhe_sha256", "source_lhe_sha256", "output_root_sha256", "seeds", "status", "model", "window"):
            with self.subTest(field=field):
                base, extra, summary = fixture()
                if field == "window":
                    extra["sampling"]["continuum_m4l_gev"] = [100.,140.]
                elif field in ("model", "status"):
                    extra["samples"][0][field] = "changed"
                else:
                    extra["samples"][0][field] = copy.deepcopy(base["samples"][0][field])
                with self.assertRaises(ValueError):
                    extension.compatible_subset(base, extra, summary)

    def test_assembled_generation_window_is_verified_in_parent(self):
        with TemporaryDirectory() as directory:
            base, extra, summary = fixture()
            path = Path(directory) / "parent.json"
            path.write_text(json.dumps(base))
            for row in base["samples"]:
                row["source_manifest"] = extension.revision.artifact(path)
            base["sampling"] = {"preset":"custom"}
            extension.compatible_subset(base, extra, summary)
            extra["sampling"]["continuum_m4l_gev"] = [100.,140.]
            with self.assertRaisesRegex(ValueError, "generation window"):
                extension.compatible_subset(base, extra, summary)

    def test_pool_preserves_one_rate_and_signed_source_variance(self):
        original = dict(sample=np.array(["bkg_qq_4e", "bkg_qq_4e", "signal", "bkg_gg"]),
            group=np.array(["original/qq", "original/qq", "original/signal", "original/gg"]),
            source=np.array([0,1,0,0]), label=np.array([0,0,1,0]), weight=np.array([3.,-1.,10.,1.]))
        extra = dict(sample=np.array(["bkg_qq_4e", "bkg_qq_4e", "bkg_qq_4e"]),
            group=np.array(["additional/qq"]*3), source=np.array([0,0,1]),
            label=np.array([0,0,0]), weight=np.array([1.,2.,-1.]))
        before = original["weight"].copy(), extra["weight"].copy()
        coefficients = {"bkg_qq_4e":{"original_fraction":.5,"additional_fraction":.5}}
        pooled = extension.pool(original, extra, coefficients)
        totals = analysis.yields(pooled, pooled["weight"], np.ones(len(pooled["weight"]), bool))
        # Two independent estimates of qq=2 yield qq=2, not qq=4.
        self.assertEqual(totals["background"]["sum_weight"], 3.)
        self.assertEqual(totals["signal"]["sum_weight"], 10.)
        self.assertEqual(totals["background"]["variance"], .25*(9+1)+.25*(9+1)+1)
        np.testing.assert_array_equal(before[0], original["weight"])
        np.testing.assert_array_equal(before[1], extra["weight"])
        extra["group"] = np.array(["original/qq"]*3)
        with self.assertRaisesRegex(ValueError, "namespaces overlap"):
            extension.pool(original, extra, coefficients)

    def test_exposure_uses_all_hard_events_not_selected_counts(self):
        base, extra, summary = fixture()
        name = extra["samples"][0]["name"]
        summary["split_fractions"][name]["test"] = .199999
        c = extension.exposure_coefficients(summary, extra)[name]
        self.assertAlmostEqual(c["original_test_exposure"], 199999.)
        self.assertAlmostEqual(c["original_fraction"], 199999./399999.)

    def test_frozen_scores_reproduce_old_test_and_ignore_training_rows(self):
        with TemporaryDirectory() as directory:
            path = Path(directory)
            data = dict(group=np.array(["original/qq", "original/qq", "original/signal"]),
                source=np.array([1,2,1]), label=np.array([0,0,1]), weight=np.array([3.,-1.,10.]))
            card = {"min_background_effective_events":100., "background_fractional_uncertainty":.1}
            summary = dict(samples=[dict(name="qq",group="qq",category="irreducible"),
                                    dict(name="signal",group="signal",category="signal")],
                card=card, split_fractions={"qq":{"test":.2},"signal":{"test":.2}},
                threshold={"threshold":.488}, xgboost_test={"totals":analysis.totals(data,data["weight"],np.ones(3,bool),card)})
            with (path/"scores.csv").open("w",newline="") as stream:
                w=csv.writer(stream)
                w.writerow(["sample","source_group","source","partition","m4l","weight","score","baseline"])
                w.writerows([["qq","qq",1,"test",125,.6,.7,True],["qq","qq",2,"test",125,-.2,.6,True],
                             ["signal","signal",1,"test",125,2.,.8,True],["qq","qq",99,"train",125,10.,.99,True]])
            result=extension.original_test(path,summary)
            np.testing.assert_allclose(result["weight"],[3.,-1.,10.])
            summary["threshold"]["threshold"] = .75
            with self.assertRaisesRegex(ValueError,"do not reproduce"):
                extension.original_test(path,summary)


if __name__ == "__main__":
    unittest.main()
