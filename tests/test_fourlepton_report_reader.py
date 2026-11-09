"""Compare chunked production reading with the original full-row calculation."""
from array import array
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np
from hfourlepton import make_fourlepton_report as report
from hfourlepton.report_reader import load_report_records


@unittest.skipUnless(all(importlib.util.find_spec(name) for name in ("ROOT", "uproot")),
                     "ROOT and uproot are required for the real-file reader comparison")
class ChunkedReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "fixture.root"
        self.sources = [dict(source_index=i, generator_weight=w, response_failure_probability=0.,
            response_closure_delta=0., accepted_lepton_count=4,
            optional_weight_names=["up", "down"], optional_weights=[w*1.2, w*.8],
            optional_selected_event_weights=[w*.2*1.2, w*.2*.8]) for i, w in ((7, 1.), (2, -.25))]
        self.rows = []
        for i in range(202):
            for source in self.sources:
                selected = i >= 200
                probability = .1 if selected else .004
                weight = source["generator_weight"] * probability
                self.rows.append(dict(source_index=source["source_index"],
                    channel=1 if selected else i % 3,
                    cut_mask=1023 if selected else (0, 1, 3, 15, 127, 255, 511)[i % 7],
                    hypothesis_probability=probability, event_weight=weight*.98,
                    pretrigger_event_weight=weight, m4l=125. + (i-200) if selected else 0.,
                    mZ1=91., mZ2=25., pt4l=30., y4l=.4, min_delta_r=.2,
                    lepton_pt=[30., 25., 20., 10.], lepton_eta=[.2, .3, -.5, .1]))

    def write(self):
        import ROOT
        output = ROOT.TFile(str(self.path), "RECREATE")
        for name, rows in (("SourceEvents", self.sources), ("FourLepton", self.rows)):
            tree = ROOT.TTree(name, name)
            buffers = {}
            for key, value in rows[0].items():
                if key.startswith("optional_"):
                    buffers[key] = ROOT.std.vector("string" if key.endswith("names") else "double")()
                    tree.Branch(key, buffers[key])
                else:
                    code = "q" if isinstance(value, int) else "d"
                    buffers[key] = array(code, [0] * (len(value) if isinstance(value, list) else 1))
                    shape = f"[{len(value)}]" if isinstance(value, list) else ""
                    tree.Branch(key, buffers[key], f"{key}{shape}/{'L' if code == 'q' else 'D'}")
            for row in rows:
                for key, value in row.items():
                    target = buffers[key]
                    if key.startswith("optional_"):
                        target.clear()
                        for item in value: target.push_back(item)
                    else:
                        for i, item in enumerate(value if isinstance(value, list) else [value]): target[i] = item
                tree.Fill()
            tree.Write()
        output.Close()
        self.path.with_suffix(".summary.json").write_text(json.dumps(dict(events_read=len(self.sources),
            four_lepton_rows=len(self.rows), output_root=str(self.path), diagnostics={})))
        manifest = self.path.with_suffix(".json")
        manifest.write_text(json.dumps(dict(schema_version=1, detector_response="ssc", samples=[dict(
            name="fixture", category="signal", channel="4e", detector_response="ssc", output_root=str(self.path),
            normalization=dict(native_cross_section_pb=1., target_cross_section_pb=1., weight_scale=1.,
                               denominator_sumw=.75))])))
        return report.load_manifest(manifest).samples[0]

    def compare(self, expected, actual, path=""):
        if isinstance(expected, dict):
            self.assertEqual(set(expected), set(actual), path)
            for key in expected: self.compare(expected[key], actual[key], path+"/"+key)
        elif isinstance(expected, list):
            self.assertEqual(len(expected), len(actual), path)
            for i, (left, right) in enumerate(zip(expected, actual)): self.compare(left, right, path+f"/{i}")
        elif isinstance(expected, float):
            self.assertTrue(math.isclose(expected, actual, rel_tol=1.e-11, abs_tol=1.e-12),
                            f"{path}: {expected} != {actual}")
        else:
            self.assertEqual(expected, actual, path)

    def test_full_report_matches_across_chunk_boundaries_and_custom_selection(self):
        sample = self.write()
        original = report._add_summary_sidecar(self.path, report._read_with_pyroot(self.path))
        for stages in (report.DEFAULT_CUT_STAGES, report.DEFAULT_CUT_STAGES[:8]):
            reference = report.analyze_sample(sample, original, stages, 100., (120., 130.), 1.e-9)
            for size in (7, 113, 10000):
                with self.subTest(size=size, stages=len(stages)):
                    compact = load_report_records(self.path, step_size=size,
                                                  selected_mask=report._final_stage(stages).required_mask)
                    self.assertEqual(compact.original_four_lepton_entries, len(self.rows))
                    self.assertLess(len(compact.four_lepton), len(self.rows)//2)
                    result = report.analyze_sample(sample, compact, stages, 100., (120., 130.), 1.e-9)
                    self.compare(reference, result)
                    expected = [r for r in original.four_lepton if report.passes_stage(r, stages[-1])]
                    actual = [r for r in compact.four_lepton if report.passes_stage(r, stages[-1])]
                    self.compare(expected, actual)

    def test_optional_weight_mismatch_remains_visible(self):
        self.sources[1]["optional_weight_names"].reverse()
        sample = self.write()
        compact = load_report_records(self.path, step_size=1)
        result = report.analyze_sample(sample, compact, report.DEFAULT_CUT_STAGES, 100., (120.,130.), 1.e-9)
        self.assertTrue(any("name/order differs" in error for error in result["audit"]["errors"]))

    def test_omitting_rejected_ml_rows_preserves_source_splits_and_signed_yields(self):
        from ho_xgboost_core import grouped_sample_split, signed_yield_summary
        self.write()
        compact = load_report_records(self.path, step_size=13)
        reference = self.rows
        selected = [r for r in compact.four_lepton if r["cut_mask"] == 1023]
        answers = []
        for hypotheses in (reference, selected):
            rows = [dict(source_index=i, event_weight=0., cut_mask=0) for i in range(9)] + hypotheses
            source = np.array([r["source_index"] for r in rows])
            group = np.array(["fixture"] * len(rows))
            weights = np.array([r["event_weight"] for r in rows])
            eligible = np.array([r["cut_mask"] == 1023 for r in rows])
            split = grouped_sample_split(group, source, fractions=(.5, .25, .25), seed=123)
            answers.append(dict(groups=split.sample_group_counts, fractions=split.sample_fractions,
                assignments=split.partitions[eligible].tolist(),
                yields=signed_yield_summary(weights, group, source, selected=eligible)))
        self.compare(answers[0], answers[1])

    def test_bad_rejected_kinematics_and_original_row_count_are_not_hidden(self):
        self.rows[0]["pt4l"] = float("nan")
        sample = self.write()
        compact = load_report_records(self.path, step_size=13)
        compact.metadata["summary"]["four_lepton_rows"] -= 1
        result = report.analyze_sample(sample, compact, report.DEFAULT_CUT_STAGES, 100., (120.,130.), 1.e-9)
        self.assertTrue(any("nonfinite pt4l" in error for error in result["audit"]["errors"]))
        self.assertTrue(any("four_lepton_rows disagrees" in error for error in result["audit"]["errors"]))


if __name__ == "__main__":
    unittest.main()
