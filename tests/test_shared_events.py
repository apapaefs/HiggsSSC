"""Portable data publication and analysis-only student access."""
import copy
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

from hgammagamma import ho_shower_completion as completion
from hfourlepton import shared_events as share
from hfourlepton.reprocess import audit_completed_manifest
import student_events as student
from test_ho_shower_completion import ShowerCompletionTests, SAMPLE


class SharedEventsTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.private = self.root / "private"
        fixture = ShowerCompletionTests().fixture(self.private)
        self.shower, _ = ShowerCompletionTests().build(fixture)
        self.publisher = share.Publisher(self.root / "shared", {"test": True})
        self.addCleanup(self.publisher.lock.close)

    def test_relocated_shower_passes_original_audit_without_private_files(self):
        before = {p: share.digest(p) for p in self.private.rglob("*") if p.is_file()}
        record = self.publisher.shower(self.shower, "nominal")
        for path, checksum in before.items():
            self.assertEqual(share.digest(path), checksum)
        self.assertEqual(record["root_summary"], self.shower["root_summary"])
        for key in completion.MOMENTS:
            self.assertEqual(record["source_weights"][key], self.shower["source_weights"][key])
        self.private.rename(self.root / "unavailable")
        manifest = dict(sample=SAMPLE, nevents_requested=3, lhe_file=record["reader_lhe_path"], lhe={"events": 3})
        completion.validate_completion(record, manifest, directory=Path(record["herwig_dir"]))
        root = Path(record["artifacts"]["roots"][0]["path"])
        root.write_bytes(b"changed events")
        with self.assertRaisesRegex(ValueError, "artifact changed"):
            completion.validate_completion(record, manifest, directory=Path(record["herwig_dir"]))

    def test_manifest_is_usable_with_private_home_removed(self):
        analysis = self.private / "analysis.root"
        analysis.write_bytes(b"analysis ROOT unchanged")
        inputs = self.private / "showers.input"
        inputs.write_text(self.shower["artifacts"]["roots"][0]["path"] + "\n")
        summary = analysis.with_suffix(".summary.json")
        share.save(summary, dict(input=str(inputs), output_root=str(analysis), events_read=3,
            sample=SAMPLE, category="signal", requested_channel="4e", seed=7,
            weight_scale=1., selection_profile="ssc", response_profile="perfect"))
        source = self.shower["source_weights"]["artifact"]
        sample = dict(name=SAMPLE, category="signal", channel="4e", status="complete", order="ho", seed=7,
            output_root=str(analysis), root_file=str(analysis), output_root_sha256=share.digest(analysis),
            analysis_summary=str(summary), analysis_summary_sha256=share.digest(summary),
            input_lhe=source["path"], source_lhe=source["path"], source_lhe_sha256=source["sha256"],
            input_lhe_sha256=source["sha256"], nevents_requested=3, shower_completion=self.shower,
            normalization=dict(native_cross_section_pb=10., target_cross_section_pb=10., weight_scale=1.,
                               denominator_sumw=1.6), selection_profile="ssc", optional_weight_names_file=None)
        manifest = self.private / "manifest_perfect.json"
        share.save(manifest, dict(schema_version=1, status="complete", order="ho", detector_profile="perfect",
            selection_profile="ssc", samples=[sample], collider={"sqrt_s_tev": 40.}))
        published = self.publisher.manifest(manifest, "manifest_perfect.json")
        self.private.rename(self.root / "unavailable")
        raw, summaries = audit_completed_manifest(published)
        self.assertEqual(raw["samples"][0]["output_root_sha256"], sample["output_root_sha256"])
        self.assertEqual(raw["samples"][0]["normalization"], sample["normalization"])
        self.assertEqual(summaries[SAMPLE]["events_read"], 3)

    def simple_dataset(self):
        root = self.publisher.root
        event = self.publisher.write("nominal/analysis/ssc/event.root", "event-data")
        self.publisher.json("manifest_ssc.json", {"samples": [SAMPLE]})
        self.publisher.json("manifest_perfect.json", {"samples": [SAMPLE]})
        cards = {}
        for mode in ("cuts", "xgboost"):
            for profile in ("ssc", "perfect"):
                key = f"{mode}_{profile}"
                cards[key] = f"cards/{key}.json"
                self.publisher.json(cards[key], {"detector_response": profile})
        share.save(root / "dataset.json", dict(schema_version=1, kind="fourlepton_shared_events", status="complete",
            samples=[SAMPLE], manifests={p: f"manifest_{p}.json" for p in ("perfect", "ssc")},
            cards=cards, files=self.publisher.files))
        return root, event

    def test_attach_copies_only_editable_cards_and_is_idempotent(self):
        dataset, _ = self.simple_dataset()
        workspace = self.root / "student"
        with mock.patch.object(student.subprocess, "run", side_effect=AssertionError("must not run a generator")):
            student.attach(dataset, workspace)
            card = workspace / "cards/cuts_ssc.json"
            card.write_text('{"student_cut": 42}\n')
            student.attach(dataset, workspace)
        self.assertEqual(json.loads(card.read_text()), {"student_cut": 42})
        self.assertTrue((workspace / "events").is_symlink())
        self.assertFalse(any(p.suffix == ".root" for p in workspace.rglob("*") if not p.is_symlink()))
        self.assertNotEqual(share.digest(card), share.digest(dataset / "cards/cuts_ssc.json"))
        student.bound(workspace)

    def test_dataset_tampering_and_escaping_paths_fail(self):
        dataset, event = self.simple_dataset()
        with self.assertRaisesRegex(ValueError, "separate directories"):
            student.attach(dataset, dataset / "student")
        event.write_text("bad---data")  # same byte count; --full must detect it
        with self.assertRaisesRegex(ValueError, "checksum changed"):
            share.load_dataset(dataset, full=True)
        with self.assertRaisesRegex(ValueError, "escapes"):
            share.safe_path(dataset, "../private/event.root")
        with self.assertRaisesRegex(ValueError, "invalid name"):
            share.name("../output")

    def test_student_upgrade_preserves_cut_edits_and_adds_new_default_samples(self):
        dataset, _ = self.simple_dataset()
        old_workspace = self.root / "old-workspace"
        relative = "cards/cuts_ssc.json"
        old_default = {"required_samples": [SAMPLE], "cuts": {"pt": 10}}
        def set_card(root, card):
            share.save(root / relative, card)
            index = json.loads((root / "dataset.json").read_text())
            index["files"][relative] = {"size": (root / relative).stat().st_size,
                                       "sha256": share.digest(root / relative)}
            share.save(root / "dataset.json", index)
        set_card(dataset, old_default)
        student.attach(dataset, old_workspace)
        edited = {**old_default, "cuts": {"pt": 22}}
        share.save(old_workspace / relative, edited)
        new_release = self.root / "updated"
        shutil.copytree(dataset, new_release)
        set_card(new_release, {**old_default, "required_samples": [SAMPLE, "bkg_ttbar"]})
        new_workspace = self.root / "new-workspace"
        student.attach(new_release, new_workspace, previous_workspace=old_workspace, channel="fourlepton")
        self.assertEqual(json.loads((new_workspace / relative).read_text()),
                         {**edited, "required_samples": [SAMPLE, "bkg_ttbar"]})
        self.assertEqual(json.loads((old_workspace / relative).read_text()), edited)
        with self.assertRaisesRegex(ValueError, "does not match"):
            student.bound(new_workspace, channel="diphoton")

    def test_diphoton_commands_read_shared_nominal_inputs_and_private_cards(self):
        root = self.root / "diphoton"
        workspace = self.root / "student"
        index = {"kind": "diphoton_shared_events", "run_tag": "nominal",
                 "cards": {"cuts_ssc": "cards/cuts_ssc.yaml", "xgboost_ssc": "cards/xgboost_ssc.yaml"}}
        for action in ("report", "cuts", "xgboost"):
            args = SimpleNamespace(action=action, profile="ssc", luminosity_fb=42., config=None)
            cmd = student.command(args, root, index, workspace, "python", workspace / "results")
            self.assertEqual(cmd[cmd.index("--analysis-root") + 1], str(root))
            self.assertEqual(cmd[cmd.index("--run-tag") + 1], "nominal")
            if action != "report":
                self.assertEqual(cmd[cmd.index("--config") + 1], str(workspace / index["cards"][action + "_ssc"]))
            else:
                self.assertIn("event_xsec", cmd)

    def test_student_commands_only_invoke_analysis_scripts(self):
        dataset, _ = self.simple_dataset()
        index = share.load_dataset(dataset)
        workspace = self.root / "student"
        student.attach(dataset, workspace)
        for mode in ("cuts", "xgboost", "report"):
            args = SimpleNamespace(action=mode, profile="ssc", config=None, luminosity_fb=100.)
            cmd = student.command(args, dataset, index, workspace, Path("python"), workspace / "results")
            self.assertIn(Path(cmd[1]).name, ("analyze_fourlepton.py", "make_fourlepton_report.py"))
            self.assertEqual(Path(cmd[cmd.index("--output-dir") + 1]), workspace / "results")
            self.assertNotIn("run_pipeline.py", " ".join(cmd))
            self.assertNotIn("Herwig", " ".join(cmd))

    def test_setup_repairs_an_interrupted_pip_bootstrap_without_generators(self):
        dataset, _ = self.simple_dataset()
        workspace = self.root / "student"
        student.attach(dataset, workspace)
        python = workspace / ".venv/bin/python"
        python.parent.mkdir(parents=True)
        python.write_bytes(b"existing Python, interrupted pip bootstrap")
        commands = []

        def run(cmd, **kwargs):
            commands.append(cmd)
            return SimpleNamespace(returncode=int(cmd[-2:] == ["pip", "--version"]))

        with mock.patch.object(student.subprocess, "run", side_effect=run):
            student.setup(workspace)
        self.assertIn([str(python), "-m", "ensurepip", "--upgrade"], commands)
        self.assertTrue(all(cmd[0] == str(python) for cmd in commands))
        self.assertFalse(any("run_pipeline.py" in arg for cmd in commands for arg in cmd))

    def test_publication_resume_refuses_changed_sources_and_never_hardlinks(self):
        source = self.private / "data"
        source.write_bytes(b"original")
        checksum = share.digest(source)
        target = self.publisher.copy_file(source, "data", checksum)
        self.assertNotEqual(source.stat().st_ino, target.stat().st_ino)
        with self.assertRaisesRegex(ValueError, "another publisher"):
            share.Publisher(self.publisher.root, {"test": True})
        self.publisher.lock.close()
        resumed = share.Publisher(self.publisher.root, {"test": True})
        self.assertEqual(resumed.copy_file(source, "data", checksum), target)
        resumed.lock.close()
        with self.assertRaisesRegex(ValueError, "different inputs"):
            share.Publisher(self.publisher.root, {"test": False})
        target.write_bytes(b"tampered")
        with share.Publisher(self.publisher.root, {"test": True}) as again:
            with self.assertRaisesRegex(ValueError, "differs"):
                again.copy_file(source, "data", checksum)


if __name__ == "__main__":
    unittest.main()
