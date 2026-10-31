"""Execute the actual MG5 bookkeeping methods on synthetic files, without MEs."""
import ast
from collections import defaultdict
import json
import logging
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

from hfourlepton import mg5_loop_recovery as fix
from hfourlepton import repair_loop_production as repair

UPSTREAM = Path(__file__).resolve().parents[1] / "MG5_aMC_v3_5_15" / fix.SOURCE


class Statistics(defaultdict):
    def __init__(self):
        super().__init__(int)

    def aggregate_statistics(self, other):
        for key, value in other.items():
            self[key] += value

    def nice_output(self, name):
        return ""


class Results(list):
    def __init__(self):
        super().__init__()
        self.run_statistics = Statistics()

    def compute_values(self, **kwargs):
        pass


class Grid:
    def __init__(self, nhel):
        self.results = Results()
        self.target_evt = 100
        self.grids = []

    def add_results_information(self, stream):
        result = SimpleNamespace(**json.load(stream))
        self.results.append(result)
        return result

    def add_one_grid_information(self, stream):
        self.grids.append(stream.read())

    def get_cross_section(self):
        xsec = sum(r.axsec for r in self.results)
        return xsec, xsec, xsec / 10

    def get_max_wgt(self, *args):
        return max(r.th_maxwgt for r in self.results)

    def get_nunwgt(self, *args):
        weights = [r.th_maxwgt for r in self.results]
        return 1 if max(weights) / min(weights) > 1e4 else 100

    def write_grid_for_submission(self, *args, **kwargs):
        self.grid_submission = args


def load_methods(source):
    """Use the real upstream functions; replace only integration/file parsers."""
    tree = ast.parse(source)
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    namespace = dict(os=os, math=math, pjoin=os.path.join, logger=logging.getLogger(__name__),
        combine_grid=SimpleNamespace(grid_information=Grid),
        misc=SimpleNamespace(mult_try_open=open, call=subprocess.call),
        files=SimpleNamespace(cp=shutil.copy2))
    # The synthetic outlier jobs omit events.lhe so MG5's optional event dump
    # is skipped. Event collection itself is exercised separately below.
    for name in ("gensym", "gen_ximprove_share"):
        methods = [node for node in classes[name].body if isinstance(node, ast.FunctionDef)
                   and node.name in ("combine_grid", "combine_iteration")]
        node = ast.ClassDef(name=name, bases=[] if name == "gensym" else [ast.Name(id="gensym", ctx=ast.Load())],
                            keywords=[], body=methods, decorator_list=[])
        module = ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[]))
        exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


def integration_fixture(root, source, weights):
    namespace = load_methods(source)
    obj = namespace["gen_ximprove_share"]()
    obj.run_card = {"nhel": 1}
    obj.splitted_for_dir = lambda p, g: len(weights)
    obj.cross, obj.abscross, obj.sigma, obj.chi2 = (defaultdict(float) for _ in range(4))
    obj.abscross[(str(root), "1")] = 1.
    obj.sigma[(str(root), "1")] = 1.
    obj.run_statistics = {}
    obj.warnings_from_statistics = lambda *args: None
    for i, weight in enumerate(weights):
        path = root / f"G1_{i+1}"
        path.mkdir()
        (path / "results.dat").write_text(json.dumps(dict(axsec=i+1., th_maxwgt=weight, nevents=100)))
        (path / "grid_information").write_text(str(i))
    return obj, namespace


class LoopRecoveryTests(unittest.TestCase):
    def test_reproduces_deleted_result_and_fixes_recursive_combination(self):
        original = UPSTREAM.read_bytes()
        for patched in (False, True):
            with tempfile.TemporaryDirectory() as tmp, self.subTest(patched=patched):
                root = Path(tmp)
                source = fix.corrected_source(original) if patched else original
                obj, _ = integration_fixture(root, source, [2., 1e8, 3.])
                if not patched:
                    with self.assertRaises(FileNotFoundError):
                        obj.combine_grid(str(root), "1", 2)
                else:
                    grid, cross, error = obj.combine_grid(str(root), "1", 2)
                    self.assertEqual(grid.higgsssc_job_indices, [0, 2])
                    self.assertEqual(grid.grids, ["0", "2"])
                    self.assertEqual([r.axsec for r in grid.results], [1., 3.])
                    self.assertEqual(grid.results.run_statistics["skipped_subchannel"], 1)
                    self.assertEqual(obj.run_statistics[(root.name, "1")]["skipped_subchannel"], 1)
                    self.assertTrue(math.isfinite(cross) and math.isfinite(error))
                    self.assertFalse(list(root.glob("G*/results.dat")))
                    self.assertFalse(list(root.glob("G*/grid_information")))

    def test_repeated_exclusions_use_original_job_indices(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = integration_fixture(root, fix.corrected_source(UPSTREAM.read_bytes()), [1e16, 1e10, 2.])
            grid, _, _ = obj.combine_grid(str(root), "1", 2)
            self.assertEqual(grid.higgsssc_job_indices, [2])
            self.assertEqual(grid.grids, ["2"])
            self.assertEqual(grid.results.run_statistics["skipped_subchannel"], 2)

    def test_normal_combination_is_unchanged(self):
        answers = []
        for source in (UPSTREAM.read_bytes(), fix.corrected_source(UPSTREAM.read_bytes())):
            with tempfile.TemporaryDirectory() as tmp:
                obj, _ = integration_fixture(Path(tmp), source, [2., 3., 4.])
                grid, cross, error = obj.combine_grid(tmp, "1", 2)
                answers.append((cross, error, [r.axsec for r in grid.results], grid.grids,
                                dict(grid.results.run_statistics)))
        self.assertEqual(answers[0][:4], answers[1][:4])
        self.assertEqual(answers[1][4].get("skipped_subchannel"), 0)

    def test_failed_read_preserves_earlier_results_for_investigation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = integration_fixture(root, fix.corrected_source(UPSTREAM.read_bytes()), [2., 3.])
            (root / "G1_2/grid_information").unlink()
            with self.assertRaises(FileNotFoundError):
                obj.combine_grid(tmp, "1", 2)
            self.assertTrue((root / "G1_1/results.dat").is_file())
            self.assertTrue((root / "G1_1/grid_information").is_file())

    def test_refine_merges_only_the_jobs_used_in_the_integral(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj, _ = integration_fixture(root, fix.corrected_source(UPSTREAM.read_bytes()), [2., 1e8, 3.])
            grid, _, _ = obj.combine_grid(tmp, "1", 2)
            for i, text in enumerate(("good first\n", "bad outlier\n", "good last\n")):
                (root / f"G1_{i+1}/events.lhe").write_text(text)
            (root / "G1").mkdir()
            obj.combine_grid = lambda *args: (grid, 1., .1)
            obj.goal_lum, obj.err_goal, obj.min_iter = 10, .1, 3
            obj.generated_events = {}
            obj.write_results = lambda *args: None
            obj.combine_iteration(tmp, "1", 2)
            self.assertEqual((root / "G1/events.lhe").read_text(), "good first\ngood last\n")

    def test_survey_merges_consistent_events_and_resubmits_all_jobs(self):
        for step in (1, 6):
            with self.subTest(step=step), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                obj, namespace = integration_fixture(root, fix.corrected_source(UPSTREAM.read_bytes()), [2., 1e8, 3.])
                grid, _, _ = obj.combine_grid(tmp, "1", 2)
                for i, text in enumerate(("first\n", "outlier\n", "last\n")):
                    (root / f"G1_{i+1}/events.lhe").write_text(text)
                (root / "G1_1/ftn25").write_text("grid")
                (root / "G1_1/log.txt").write_text("log")
                obj.combine_grid = lambda *args: (grid, 1., .1)
                obj.cmd = SimpleNamespace(opts={"iterations": 5, "accuracy": .05})
                obj.min_iterations, obj.mode, obj.me_dir = 3, "survey", tmp
                calls = []
                obj.resubmit_survey = lambda *args: calls.append(args)
                obj.write_results = lambda *args: None
                namespace["gensym"].combine_iteration(obj, tmp, "1", step)
                if step == 1:
                    self.assertEqual(calls[0][2], [str(root / f"G1_{i}") for i in (1, 2, 3)])
                else:
                    self.assertEqual((root / "G1/events.lhe").read_text(), "first\nlast\n")

    def test_install_is_version_checked_idempotent_and_does_not_edit_exports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "mg5" / fix.SOURCE
            source.parent.mkdir(parents=True)
            original = UPSTREAM.read_bytes()
            source.write_bytes(original)
            campaign = root / "campaign"
            exported = campaign / "production/done/mg5_process/bin/internal/gen_ximprove.py"
            exported.parent.mkdir(parents=True)
            exported.write_bytes(original)
            fix.prepare_campaign(root / "mg5", campaign)
            self.assertEqual(exported.read_bytes(), original)
            self.assertEqual(source.read_bytes(), fix.corrected_source(original))
            record = json.loads((campaign / "mg5_loop_recovery.json").read_text())
            self.assertEqual(record["exports_already_present"]["done"]["sha256"], fix.UPSTREAM_SHA256)
            before = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in root.rglob("*") if p.is_file()}
            fix.prepare_campaign(root / "mg5", campaign)
            self.assertEqual(before, {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in before})
            source.write_bytes(source.read_bytes() + b"\n# unexpected change\n")
            with self.assertRaisesRegex(ValueError, "source changed"):
                fix.install(root / "mg5")
            with self.assertRaisesRegex(ValueError, "unrecognized MG5"):
                fix.corrected_source(b"unknown version")

    def test_repair_recognizes_only_the_specific_failure_and_never_replaces_final_lhe(self):
        with tempfile.TemporaryDirectory() as tmp:
            group = Path(tmp) / "failed"
            proc = group / "mg5_process"
            source = proc / "bin/internal/gen_ximprove.py"
            source.parent.mkdir(parents=True)
            source.write_bytes(UPSTREAM.read_bytes())
            debug = proc / "test_tag_1_debug.log"
            debug.write_text("gen_ximprove_share.combine_iteration FileNotFoundError G1_1/results.dat")
            discarded = proc / "SubProcesses/P0/DiscardedUnstableEvents/discarded_G1.dat"
            discarded.parent.mkdir(parents=True)
            with self.assertRaisesRegex(ValueError, "no unstable-event"):
                repair.known_failure(group, "test")
            discarded.write_text("outlier evidence")
            report = repair.known_failure(group, "test")
            self.assertEqual(report["mg5_source_sha256"], fix.UPSTREAM_SHA256)
            final = proc / "Events/test/unweighted_events.lhe.gz"
            final.parent.mkdir(parents=True)
            final.write_bytes(b"existing event artifact")
            with self.assertRaisesRegex(ValueError, "final production artifact"):
                repair.known_failure(group, "test")
            self.assertEqual(final.read_bytes(), b"existing event artifact")

    def test_archive_moves_only_failed_group_and_preserves_every_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            group = root / "production/failed"
            group.mkdir(parents=True)
            (group / "partial-grid.dat").write_text("recoverable integration grid")
            sibling = root / "production/completed.lhe"
            sibling.write_text("completed events")
            metadata = root / "manifest.json"
            metadata.write_text("original manifest")
            before = {p.name: p.read_bytes() for p in group.iterdir()}
            destination = root / "recovery/failed"
            repair.archive(group, destination, {"generators_launched": False}, [metadata])
            self.assertFalse(group.exists())
            self.assertEqual(before, {p.name: p.read_bytes() for p in (destination / "failed_production").iterdir()})
            self.assertEqual(sibling.read_text(), "completed events")
            self.assertEqual(metadata.read_text(), "original manifest")
            with self.assertRaises(FileExistsError):
                repair.archive(group, destination, {}, [metadata])


if __name__ == "__main__":
    unittest.main()
