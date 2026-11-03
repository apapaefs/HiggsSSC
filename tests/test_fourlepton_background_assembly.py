import copy
import unittest
from hfourlepton import replace_qq_backgrounds as assembly


def fixtures():
    base = dict(samples=[], collider={"sqrt_s_tev": 13.}, selection_profile="cms",
                detector_profile="perfect", order="ho", generation_cuts={"mass_window": [90,180]})
    for name in ["signal", *sorted(assembly.QQ)]:
        base["samples"].append(dict(name=name, selection_version="cms_fiducial_v3", ancestry_schema_version=1,
            nevents_requested=20000, input_lhe_sha256=name + "old", normalization={"target_cross_section_pb": 1.}))
    replacement = copy.deepcopy(base)
    replacement["samples"] = [s for s in replacement["samples"] if s["name"] in assembly.QQ]
    for s in replacement["samples"]:
        s.update(nevents_requested=200000, input_lhe_sha256=s["name"]+"new")
        s["normalization"]["target_cross_section_pb"] = 1.01
    return base, replacement


class BackgroundAssemblyTest(unittest.TestCase):
    def test_replacement_does_not_double_count_rates_or_mutate_inputs(self):
        base, replacement = fixtures()
        original = copy.deepcopy((base, replacement))
        result = assembly.assemble(base, replacement, {"base": "base", "replacement": "new"})
        self.assertEqual((base, replacement), original)
        self.assertEqual(len(result["samples"]), 4)
        for s in result["samples"]:
            self.assertEqual(s["normalization"]["target_cross_section_pb"], 1.01 if s["name"] in assembly.QQ else 1.)
            self.assertEqual(s["source_manifest"], "new" if s["name"] in assembly.QQ else "base")
        self.assertEqual(result["sampling"]["events_by_sample"]["bkg_qq_4e"], 200000)

    def test_rejects_incompatible_or_reused_population(self):
        for mutation in (
            lambda r: r["samples"].pop(),
            lambda r: r.update(generation_cuts={"mass_window": [100,150]}),
            lambda r: r["samples"][0].update(selection_version="cms_fiducial_v2"),
            lambda r: r["samples"][0].update(pdf_lhaid=999),
            lambda r: r["samples"][0].update(nevents_requested=20000),
            lambda r: r["samples"][0].update(input_lhe_sha256=r["samples"][0]["name"]+"old"),
        ):
            base, replacement = fixtures()
            mutation(replacement)
            with self.assertRaises(ValueError):
                assembly.assemble(base, replacement, {"base": "base", "replacement": "new"})


if __name__ == "__main__":
    unittest.main()
