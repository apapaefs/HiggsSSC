import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
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


def ssc_fixtures(profile="ssc"):
    base, replacement = fixtures()
    for manifest in (base, replacement):
        manifest.update(collider={"sqrt_s_tev": 40.}, selection_profile="ssc", detector_profile=profile,
                        detector={"profile": profile, "muon_resolution_scale": 1.},
                        sampling={"continuum_m4l_gev": [90., 180.]})
        for sample in manifest["samples"]:
            sample.update(selection_version="ssc_v1", ancestry_schema_version=0,
                          source_lhe_sha256=sample["input_lhe_sha256"],
                          seeds={"production": 1 if manifest is base else 2, "analysis": 3},
                          shower_completion={"artifacts": {"roots": [{"path": sample["input_lhe_sha256"]+".root",
                                                                    "sha256": sample["input_lhe_sha256"]}]}})
    return base, replacement


class BackgroundAssemblyTest(unittest.TestCase):
    def test_finishing_training_does_not_imply_supported_statistics(self):
        cuts = {"cuts_full_sample": {"totals": {"statistically_supported": True,
                 "background": {"effective_count": 1000.}}}}
        ml = {"threshold": {"threshold": .5}, "xgboost_test": {"totals": {
              "statistically_supported": False, "background": {"effective_count": 90.}}}}
        self.assertFalse(assembly.statistical_support(cuts, ml)["ml_supported"])
        ml["xgboost_test"]["totals"]["statistically_supported"] = True
        self.assertTrue(assembly.statistical_support(cuts, ml)["ml_supported"])
        self.assertFalse(assembly.statistical_support(cuts, {"threshold": None})["ml_supported"])

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

    def test_ssc_replacement_preserves_signal_and_ggf_normalization(self):
        base, replacement = ssc_fixtures()
        before = copy.deepcopy((base, replacement))
        result = assembly.assemble(base, replacement, {"base": "old", "replacement": "new"})
        self.assertEqual((base, replacement), before)
        signal = result["samples"][0]
        self.assertEqual({k: v for k, v in signal.items() if k != "source_manifest"}, base["samples"][0])
        self.assertEqual(result["detector_profile"], "ssc")
        self.assertEqual(len(result["samples"]), len(base["samples"]))
        self.assertIn("SSC 40 TeV", result["sampling"]["description"])

    def test_ssc_rejects_changed_detector_selection_window_and_seed(self):
        for mutation in (
            lambda r: r["detector"].update(muon_resolution_scale=1.2),
            lambda r: r["samples"][0].update(ancestry_schema_version=1),
            lambda r: r["samples"][0].update(photon_dressing_dr=.3),
            lambda r: r["samples"][0]["seeds"].update(production=1),
            lambda r: r["sampling"].update(continuum_m4l_gev=None),
        ):
            base, replacement = ssc_fixtures()
            mutation(replacement)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                assembly.assemble(base, replacement, {"base": "old", "replacement": "new"})

    def test_perfect_ssc_comparison_requires_shared_hard_events_and_showers(self):
        ssc, _ = ssc_fixtures()
        perfect, _ = ssc_fixtures("perfect")
        assembly.audit_profile_pair(ssc, perfect)
        for field in ("source_lhe_sha256", "input_lhe_sha256", "seeds", "nevents_requested", "shower_completion"):
            modified = copy.deepcopy(perfect)
            modified["samples"][0][field] = None
            with self.subTest(field=field), self.assertRaises(ValueError):
                assembly.audit_profile_pair(ssc, modified)


class AssemblyExecutionTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.output = self.root / "combined"
        self.inputs = {}
        for profile in ("ssc", "perfect"):
            for folder, manifest in zip(("base", "new"), ssc_fixtures(profile)):
                path = self.root / folder / f"manifest_{profile}.json"
                path.parent.mkdir(exist_ok=True)
                path.write_text(json.dumps(manifest))
                self.inputs[path] = manifest
        self.before = {p: p.read_bytes() for p in self.inputs}
        self.calls = []

    def fake_command(self, command, log):
        command = [str(x) for x in command]
        self.calls.append(command)
        work = Path(command[command.index("--output-dir")+1])
        self.assertFalse(any(work.iterdir()))
        (work / "index.html").write_text("synthetic result")
        if "make_fourlepton_report.py" in command[1]:
            self.assertEqual(command[2:4], [str(self.output / f"manifest_{p}.json") for p in ("perfect", "ssc")])
            result, filename = {"audit": {"ok": True, "errors": []}}, "report.json"
        else:
            self.assertTrue(command[1].endswith("analyze_fourlepton.py"))
            self.assertEqual(command[command.index("--manifest")+1], str(self.output / "manifest_ssc.json"))
            card = json.loads(Path(command[command.index("--config")+1]).read_text())
            self.assertEqual(card["detector_response"], "ssc")
            self.assertEqual(card["selection_profile"], "ssc")
            minimum = 400. if command[2] == "cuts" else 100.
            self.assertEqual(card["min_background_effective_events"], minimum)
            result = {"cuts_full_sample": {"totals": {"statistically_supported": True,
                       "background": {"effective_count": 700.}}}, "threshold": None, "xgboost_test": None}
            filename = "summary.json"
        (work / filename).write_text(json.dumps(result))
        return 0

    def execute(self):
        with mock.patch.object(assembly.revision, "audit_completed_manifest",
                               side_effect=lambda p: (self.inputs[p], {})), \
             mock.patch.object(assembly.revision, "run_logged", side_effect=self.fake_command):
            assembly.execute(self.root / "base/manifest_ssc.json", self.root / "new/manifest_ssc.json",
                             self.output, 100., paired_perfect=True)

    def test_paired_ssc_runs_report_cuts_ml_preserves_inputs_and_resumes(self):
        self.execute()
        self.assertEqual(len(self.calls), 3)
        state = json.loads((self.output / "state.json").read_text())
        self.assertEqual(state["validation_status"], "not_applicable_ssc")
        self.assertEqual(state["statistics_status"], "more_statistics_needed")
        self.assertFalse((self.output / "lhc_validation").exists())
        self.execute()
        self.assertEqual(len(self.calls), 3)
        self.assertEqual({p: p.read_bytes() for p in self.inputs}, self.before)

    def test_modified_finished_product_cannot_be_silently_overwritten(self):
        self.execute()
        (self.output / "xgboost/summary.json").write_text("changed")
        with self.assertRaisesRegex(ValueError, "changed artifact"):
            self.execute()

    def test_pair_mismatch_stops_before_writing_or_launching(self):
        self.inputs[self.root / "new/manifest_perfect.json"]["samples"][0]["input_lhe_sha256"] = "different"
        with self.assertRaisesRegex(ValueError, "perfect/SSC"):
            self.execute()
        self.assertFalse(self.calls)
        self.assertFalse(self.output.exists())

    def test_cms_still_runs_acceptance_check_and_uses_perfect_cards(self):
        base, replacement = fixtures()
        paths = [self.root / folder / "manifest_perfect.json" for folder in ("base", "new")]
        for path, manifest in zip(paths, (base, replacement)):
            path.write_text(json.dumps(manifest))
            self.inputs[path] = manifest

        def run(command, log):
            command = [str(x) for x in command]
            self.calls.append(command)
            work = Path(command[command.index("--output-dir")+1])
            (work / "index.html").write_text("synthetic result")
            result = {}
            if command[1].endswith("make_fourlepton_report.py"):
                filename = "report.json"
            elif command[1].endswith("validate_lhc.py"):
                filename, result = "comparison.json", {"status": "not_validated"}
            else:
                card = json.loads(Path(command[command.index("--config")+1]).read_text())
                self.assertEqual(card["detector_response"], "perfect")
                self.assertEqual(card["selection_profile"], "cms")
                filename = "summary.json"
                result = {"cuts_full_sample": {"totals": {"statistically_supported": False,
                           "background": {"effective_count": 30.}}}, "threshold": None}
            (work / filename).write_text(json.dumps(result))
            return 2 if filename == "comparison.json" else 0

        with mock.patch.object(assembly.revision, "audit_completed_manifest",
                               side_effect=lambda p: (self.inputs[p], {})), \
             mock.patch.object(assembly.revision, "run_logged", side_effect=run):
            assembly.execute(*paths, self.output, 100.)
        self.assertEqual(len(self.calls), 4)
        state = json.loads((self.output / "state.json").read_text())
        self.assertEqual(state["validation_status"], "not_validated")
        self.assertTrue((self.output / "manifest_perfect.json").is_file())
        self.assertFalse((self.output / "manifest_ssc.json").exists())


if __name__ == "__main__":
    unittest.main()
