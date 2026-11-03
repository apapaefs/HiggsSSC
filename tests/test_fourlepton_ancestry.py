import json
from pathlib import Path
import tempfile
import unittest

from hfourlepton import ancestry_shower, reprocess


class AncestryShowerTest(unittest.TestCase):
    def test_validation_description_requires_both_version_and_schema(self):
        from hfourlepton.validate_lhc import ancestry_aware
        truth = dict(selection_version="cms_fiducial_v3", ancestry_schema_version=1)
        self.assertTrue(ancestry_aware([truth]))
        self.assertFalse(ancestry_aware([]))
        self.assertFalse(ancestry_aware([truth, dict(selection_version="cms_fiducial_v2")]))
        self.assertFalse(ancestry_aware([dict(selection_version="cms_fiducial_v3")]))

    def test_new_campaign_uses_private_binary_and_ancestry(self):
        from hfourlepton import run_four_lepton_campaign as four, higher_order
        cfg = four.parse_config(["--order", "ho", "--selection-profile", "cms", "--run-tag", "truth_test",
                                 "--run-samples", "bkg_qq", "--ancestry-plugin-dir", "/private/truth", "--dry-run"])
        sample = four.selected_samples(cfg)[0]
        self.assertEqual(cfg.analysis_exe.parent, cfg.run_dir / "code")
        card = higher_order.herwig_card(sample, Path("/source.lhe"), 123, 30, cfg)
        self.assertIn("library /private/truth/HwSimAncestry.so", card)
        self.assertIn("SaveAncestry Yes", card)
        self.assertIn("--cms-ancestry", four.analysis_command(sample, Path("/roots.input"), 1, 1., cfg))

    def test_card_changes_only_plugin_and_output_schema(self):
        original = ("library HwSim.so\nset /Herwig/Analysis/HwSim:SaveObjects Yes\n"
                    "set EventGenerator:RandomNumberGenerator:Seed 123\n"
                    "set LesHouchesReader:AllowedToReOpen No\n")
        changed = ancestry_shower.ancestry_card(original, Path("/private/HwSimAncestry.so"))
        restored = changed.replace("library /private/HwSimAncestry.so", "library HwSim.so").replace(
            "set /Herwig/Analysis/HwSim:SaveAncestry Yes\n", "")
        self.assertEqual(restored, original)
        for bad in (original + "library HwSim.so\n", changed, original.replace("SaveObjects Yes", "SaveObjects No")):
            with self.assertRaises(ValueError):
                ancestry_shower.ancestry_card(bad, Path("/private/plugin.so"))

    def test_plugin_and_source_checksums_are_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            library = folder / "HwSimAncestry.so"
            source = folder / "HwSim.cc"
            library.write_bytes(b"plugin")
            source.write_text("source")
            proof = dict(schema_version=1, plugin_sha256=reprocess.digest(library),
                         source_sha256={source.name: reprocess.digest(source)})
            (folder / "build.json").write_text(json.dumps(proof))
            self.assertEqual(ancestry_shower.plugin_artifacts(folder)["library"]["sha256"], proof["plugin_sha256"])
            source.write_text("changed")
            with self.assertRaisesRegex(ValueError, "changed artifact"):
                ancestry_shower.plugin_artifacts(folder)

    def test_ancestry_revision_rejects_missing_schema_and_changed_weights(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            root = folder / "result.root"
            root.write_bytes(b"root")
            sample = dict(name="s", normalization={"weight_scale": 1})
            before = dict(events_read=2, input="old.input", generated_raw={"sumw": 2.},
                          generated_scaled={"sumw": 2.})
            after = {**before, "input": "new.input", "selection_version": "cms_fiducial_v3",
                     "photon_dressing_dr": .3, "alternative_pair_veto": True,
                     "output_root": str(root), "selected_scaled": {"sumw": 1.}, "cut_mask_definition": "test"}
            path = folder / "summary.json"
            reprocess.save(path, after)
            with self.assertRaisesRegex(ValueError, "ancestry schema"):
                reprocess.updated_sample(sample, before, path, [], selection="cms_fiducial_v3", input_path="new.input")
            after["ancestry_schema_version"] = 1
            reprocess.save(path, after)
            entry = reprocess.updated_sample(sample, before, path, [], selection="cms_fiducial_v3", input_path="new.input")
            self.assertEqual(entry["normalization"], sample["normalization"])
            after["generated_raw"] = {"sumw": 3.}
            reprocess.save(path, after)
            with self.assertRaisesRegex(ValueError, "source weights"):
                reprocess.updated_sample(sample, before, path, [], selection="cms_fiducial_v3", input_path="new.input")


if __name__ == "__main__":
    unittest.main()
