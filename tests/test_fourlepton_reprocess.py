import copy
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

from hfourlepton import reprocess as revision
import run_pipeline as cli


class ReprocessTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.parent = self.root / "manifest_perfect.json"
        self.output = self.root / "revisions/cms_v2"
        self.summaries = {}
        self.raw = dict(order="ho", selection_profile="cms", detector_profile="perfect",
            collider={"sqrt_s_tev": 13.}, run_tag="original", configuration_fingerprint="generation-proof",
            software={}, tracked_inputs={}, sampling={"preset": "minimal"},
            cut_stages=[{"bit": 5, "label": "Topology"}], samples=[])
        for name in ("first", "second"):
            source = self.root / (name + "_shower.root")
            source.write_bytes(b"original shower")
            inputs = self.root / (name + ".input")
            inputs.write_text(str(source) + "\n")
            original = self.root / (name + "_analysis.root")
            original.write_bytes(b"original analysis")
            weights = dict(entries=3, sumw=1., sumw2=3., sumabsw=3.)
            before = dict(input=str(inputs), output_root=str(original), events_read=3, seed=123,
                sample=name, category="signal", requested_channel="4e", weight_scale=2.,
                muon_resolution_scale=1., selection_profile="cms", response_profile="perfect",
                optional_weight_names_file="", cut_mask_definition="old cuts",
                generated_raw=weights, generated_scaled=weights, selected_scaled=weights)
            summary = original.with_suffix(".summary.json")
            revision.save(summary, before)
            s = dict(name=name, category="signal", channel="4e", status="complete", seed=123,
                output_root=str(original), root_file=str(original), output_root_sha256=revision.digest(original),
                analysis_summary=str(summary), analysis_summary_sha256=revision.digest(summary),
                normalization={"denominator_sumw": 1., "target_cross_section_pb": 2., "weight_scale": 2.},
                optional_weight_names_file=None)
            self.raw["samples"].append(s)
            self.summaries[name] = before
        revision.save(self.parent, self.raw)
        self.original_hashes = {p: revision.digest(p) for p in self.root.iterdir() if p.is_file()}
        self.launched = []
        self.downstream = []

    def fake_reconstruct(self, sample, before, raw, output, executable):
        self.launched.append(sample["name"])
        folder = revision.attempt_directory(output / "analysis/perfect" / sample["name"])
        after = copy.deepcopy(before)
        after.update(output_root=str(folder / (sample["name"] + ".root")),
                     selection_version=revision.SELECTION, photon_dressing_dr=.3,
                     alternative_pair_veto=True, cut_mask_definition="corrected cuts")
        Path(after["output_root"]).write_bytes(b"new analysis")
        path = Path(after["output_root"]).with_suffix(".summary.json")
        revision.save(path, after)
        return revision.updated_sample(sample, before, path, [str(executable)])

    def fake_command(self, command, log):
        command = [str(x) for x in command]
        if command[0] == "make":
            (Path(command[command.index("-C") + 1]) / revision.TARGET).write_bytes(b"test executable")
            return 0
        work = Path(command[command.index("--output-dir") + 1])
        self.assertFalse(any(work.iterdir()), "downstream CLIs require empty output directories")
        self.downstream.append(Path(command[1]).name)
        (work / "index.html").write_text("test result")
        if "validate_lhc.py" in command[1]:
            revision.save(work / "comparison.json", {"status": "not_validated"})
            return 2
        revision.save(work / ("report.json" if "make_fourlepton_report.py" in command[1] else "summary.json"), {})
        return 0

    def execute(self, reconstruct=None, check=False):
        with mock.patch.object(revision, "audit_parent", return_value=(self.raw, self.summaries)), \
             mock.patch.object(revision.importlib, "import_module"), \
             mock.patch.object(revision, "run_logged", side_effect=self.fake_command), \
             mock.patch.object(revision, "reconstruct", side_effect=reconstruct or self.fake_reconstruct):
            revision.execute(self.parent, self.output, 2, 100., check=check)

    def test_revision_preserves_originals_and_resumes_without_execution(self):
        self.execute()
        state = json.loads((self.output / "state.json").read_text())
        self.assertEqual(state["state"], "complete")
        self.assertEqual(state["validation_status"], "not_validated")
        self.assertIn("analyze_fourlepton.py", self.downstream)
        revised = json.loads((self.output / "manifest_perfect.json").read_text())
        self.assertEqual(revised["analysis_revision"]["parent_configuration_fingerprint"], "generation-proof")
        for old, new in zip(self.raw["samples"], revised["samples"]):
            self.assertEqual(old["normalization"], new["normalization"])
            self.assertEqual(new["selection_version"], revision.SELECTION)
            self.assertNotEqual(old["output_root"], new["output_root"])
        for name, minimum in (("cuts", 400.), ("xgboost", 100.)):
            card = json.loads((self.output / "analysis_cards" / (name + ".json")).read_text())
            self.assertEqual(card["min_background_effective_events"], minimum)
        launched, downstream = list(self.launched), list(self.downstream)
        self.execute()
        self.assertEqual(self.launched, launched)
        self.assertEqual(self.downstream, downstream)
        for path, checksum in self.original_hashes.items():
            self.assertEqual(revision.digest(path), checksum, str(path))

    def test_partial_resume_only_reconstructs_missing_sample(self):
        def fail_one(s, *args):
            if s["name"] == "second":
                raise RuntimeError("synthetic worker failure")
            return self.fake_reconstruct(s, *args)
        with self.assertRaisesRegex(RuntimeError, "synthetic worker failure"):
            self.execute(fail_one)
        self.assertFalse((self.output / "manifest_perfect.json").exists())
        self.execute()
        self.assertEqual(self.launched.count("first"), 1)
        self.assertEqual(self.launched.count("second"), 1)

    def test_modified_completed_outputs_and_parent_are_rejected(self):
        self.execute()
        parent_text = self.parent.read_text()
        self.parent.write_text(parent_text + "\n")
        with self.assertRaisesRegex(ValueError, "inputs or software changed"):
            self.execute()
        self.parent.write_text(parent_text)
        path = self.output / "cuts/summary.json"
        path.write_text("changed")
        with self.assertRaisesRegex(ValueError, "changed artifact"):
            self.execute()

    def test_check_is_read_only_and_unmanaged_directory_is_rejected(self):
        self.execute(check=True)
        self.assertFalse(self.output.exists())
        self.output.mkdir(parents=True)
        (self.output / "user-file").write_text("preserve me")
        with self.assertRaisesRegex(ValueError, "unmanaged nonempty"):
            self.execute()

    def test_changed_source_weights_cannot_be_published(self):
        s = self.raw["samples"][0]
        entry = self.fake_reconstruct(s, self.summaries[s["name"]], self.raw, self.output, Path("analyzer"))
        path = Path(entry["analysis_summary"])
        after = json.loads(path.read_text())
        after["generated_raw"]["sumw"] += 1
        revision.save(path, after)
        with self.assertRaisesRegex(ValueError, "changed source weights"):
            revision.updated_sample(s, self.summaries[s["name"]], path, [])

    def test_validation_crash_is_not_a_physics_review_result(self):
        with mock.patch.object(revision, "run_logged", return_value=1):
            with self.assertRaisesRegex(ValueError, "valid comparison"):
                revision.stage(self.output, "lhc_validation", "id", lambda work: [], validation=True)
        self.assertFalse((self.output / "lhc_validation").exists())

    def test_parent_audit_rejects_changed_input_inventory(self):
        samples = [SimpleNamespace(name=s["name"], raw=s, output_root=Path(s["output_root"])) for s in self.raw["samples"]]
        bundle = SimpleNamespace(samples=samples, audit=SimpleNamespace(errors=[]))
        def completed(sample, **_):
            return {"saved_events": 3, "artifacts": {"roots": [{"path": str(self.root / (sample.name + "_shower.root"))}]}}
        with mock.patch.object(revision.report, "load_manifest", return_value=bundle), \
             mock.patch.object(revision.report, "_expected_runtime_metadata", return_value={}), \
             mock.patch.object(revision.report, "_validated_shower_completion", side_effect=completed):
            revision.audit_parent(self.parent)
            (self.root / "first.input").write_text(str(self.root / "second_shower.root") + "\n")
            with self.assertRaisesRegex(ValueError, "inventory differs"):
                revision.audit_parent(self.parent)

    def test_managed_command_never_updates_original_campaign_state(self):
        args = cli.parse_args(["reprocess", "--workflow", "fourlepton-ho", "--energy-tev", "13",
                              "--run-tag", "original", "--revision", "cms_v2"])
        pipeline = SimpleNamespace(args=args, python=Path("python"))
        with mock.patch.object(cli, "physics_artifacts"), mock.patch.object(cli.shared.Pipeline, "call") as call, \
             mock.patch.object(cli, "progress") as progress:
            cli.reprocess(pipeline)
        self.assertIn("hfourlepton/reprocess.py", str(call.call_args.args))
        progress.assert_not_called()
        with self.assertRaises(SystemExit):
            cli.parse_args(["reprocess", "--workflow", "fourlepton-ho", "--run-tag", "original",
                            "--revision", "cms_v2"])

    def test_analysis_rejects_mixed_selection_versions_before_reading_events(self):
        from hfourlepton import analyze_fourlepton as analysis
        samples = [SimpleNamespace(name="first", raw={"selection_version": "cms_fiducial_v2"}),
                   SimpleNamespace(name="second", raw={})]
        bundle = SimpleNamespace(samples=samples, detector_profile="perfect", audit=SimpleNamespace(errors=[]))
        card = dict(detector_response="perfect", selection_profile="cms", required_samples=["first", "second"])
        with mock.patch.object(analysis.report, "load_manifest", return_value=bundle), \
             mock.patch.object(analysis, "audit_disjoint_slices"), \
             mock.patch.object(analysis.report, "load_root_records") as read:
            with self.assertRaisesRegex(ValueError, "different selection versions"):
                analysis.load_dataset(self.parent, card, 100.)
            read.assert_not_called()


if __name__ == "__main__":
    unittest.main()
