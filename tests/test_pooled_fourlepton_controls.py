import json
import tempfile
import unittest
from pathlib import Path

from hfourlepton import make_fourlepton_report as report
from hfourlepton import pooled_controls as controls
from hfourlepton import shared_events
from hfourlepton import student_backgrounds


class PooledControlTest(unittest.TestCase):
    stages = [report.CutStage("objects", "Objects", 0, 1), report.CutStage("selected", "Selected", 1, 3)]

    def shard(self, i):
        candidate = dict(source_index=2*i, cut_mask=3, channel=1, event_weight=0.4, m4l=121.+4*i)
        selected = [candidate, dict(candidate, event_weight=0.2)]
        selected_stats = report.grouped_weight_stats(selected).as_dict()
        upstream = report.WeightStats()
        upstream.add(1.)
        upstream.add(1.)
        return dict(name="bkg_ttbar", label="ttbar", profile="ssc", order="LO", source_kind="LHE",
            native_pb=4., target_pb=4., category="reducible", source_range=[2*i,2*i+2],
            selected=selected, selected_weight=selected_stats, cutflows={"all":[upstream.as_dict(),selected_stats]},
            cuts=[dict(name=s.name,bit=s.bit,required_mask=s.required_mask) for s in self.stages],
            partition=dict(index=i,start=2*i,stop=2*i+2,events=2,source_sha256="parent",sha256=str(i)),
            manifest=f"shard-{i}", manifest_sha256=f"manifest-{i}", root_sha256=f"root-{i}",
            shower_fingerprint=f"shower-{i}", generation_seeds=dict(herwig=i,analysis=i,production=3),
            normalization=dict(denominator_sumw=2.,denominator_sumw2=2.,denominator_sumabsw=2.,denominator_event_count=2,
                denominator_npositive=2, denominator_nnegative=0, denominator_nzero=0),
            generated=upstream.as_dict(),closure=dict(ok=True,source_events=2,evaluated=2,missing=0,failed=0,
                direct_delta_disagreements=0, max_abs_delta=0., mean_abs_delta=0., tolerance=1e-9),
            audit=dict(ok=True,warnings=[],errors=[]),optional=dict(variation_names=[],variations=[],
                source_events=2,events_with_optional_branches=2,events_with_variations=0,synthesized_name_events=0))

    def test_pool_normalizes_parent_once_and_groups_correlated_hypotheses(self):
        sample = controls.pool([self.shard(0),self.shard(1)],4)
        self.assertEqual(sample["normalization"]["denominator_sumw"],4.)
        result = controls.sample_report(sample,Path("compact.json"),self.stages,2.,(124.,126.),1e-9,(120.,130.),2.,False)
        self.assertAlmostEqual(result["selected"]["events"],2400.)
        self.assertAlmostEqual(result["selected"]["raw"]["sumw2"],2*(.4+.2)**2)
        self.assertAlmostEqual(result["signal_region"]["events"],1200.)
        self.assertEqual(result["cutflow"][0]["entries"],4)
        self.assertEqual(result["cutflow"][0]["yield_events"],8000.)
        self.assertEqual(result["comparison_identity"]["seeds"]["herwig"],[0,1])
        self.assertEqual(result["generated"]["entries"],4)

    def test_incomplete_duplicate_or_overlapping_shards_fail(self):
        first, second = self.shard(0),self.shard(1)
        with self.assertRaises(ValueError):
            controls.pool([first],4)
        with self.assertRaises(ValueError):
            controls.pool([first,first],4)
        second["source_range"]=[1,3]
        with self.assertRaises(ValueError):
            controls.pool([first,second],4)

    def test_declared_default_requires_all_controls_and_pinned_products(self):
        sample = controls.pool([self.shard(0)],2)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            refs = {}
            for name in controls.PROCESSES:
                path = root / f"{name}.json"
                path.write_text(json.dumps(dict(sample,name=name)))
                refs[name] = dict(path=path.name,sha256=shared_events.digest(path))
            manifest = root / "manifest.json"
            manifest.write_text(json.dumps(dict(detector_profile="ssc",pooled_backgrounds=refs)))
            self.assertEqual(len(controls.load(manifest,self.stages)),3)
            (root / "bkg_ttbar.json").write_text("{}")
            with self.assertRaisesRegex(ValueError,"changed"):
                controls.load(manifest,self.stages)
            del refs["bkg_ttbar"]
            manifest.write_text(json.dumps(dict(detector_profile="ssc",pooled_backgrounds=refs)))
            with self.assertRaisesRegex(ValueError,"incomplete"):
                controls.load(manifest,self.stages)

    def test_packed_controls_match_public_bytes_and_detector_profile(self):
        import gzip
        import hashlib
        index = json.loads((student_backgrounds.DIRECTORY / "manifest.json").read_bytes())
        self.assertEqual(index["release"], "ssc40-20261118")
        for profile, response in index["profiles"].items():
            raw = dict(order="ho", collider={"sqrt_s_tev": 40}, selection_profile="ssc", detector_profile=profile,
                       detector=response["detector"], samples=[{"selection_version": "ssc_v1"}])
            refs = student_backgrounds.controls_for(raw)
            self.assertEqual(set(refs), set(controls.PROCESSES))
            for name, info in refs.items():
                data = (student_backgrounds.DIRECTORY / info["path"]).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), info["sha256"])
                public = gzip.decompress(data)
                self.assertEqual(hashlib.sha256(public).hexdigest(), info["public_sha256"])
                sample = json.loads(public)
                self.assertEqual((sample["name"], sample["detector_profile"]), (name, profile))
                self.assertEqual(sample["normalization"]["denominator_event_count"], 1_000_000)
                self.assertTrue(sample["audit"]["ok"])
            with self.assertRaisesRegex(ValueError, "detector settings"):
                student_backgrounds.controls_for(dict(raw, detector={"profile": profile}))
            with self.assertRaisesRegex(ValueError, "40 TeV"):
                student_backgrounds.controls_for(dict(raw, collider={"sqrt_s_tev": 13}))


if __name__ == "__main__":
    unittest.main()
