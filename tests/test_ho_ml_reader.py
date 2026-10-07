"""Strict HO source-event metadata reads without requiring a PyROOT install."""

import importlib.util
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


READER_SOURCE = Path(__file__).resolve().parents[1] / "read_root_varfiles.py"


class FakeLeaf:
    def __init__(self, length, type_name="Double_t", counted=False):
        self.length, self.type_name, self.counted = length, type_name, counted

    def GetLenStatic(self):
        return self.length

    def GetTypeName(self):
        return self.type_name

    def GetLeafCount(self):
        return object() if self.counted else None


class FakeBranch:
    def __init__(self, tree, name):
        self.tree, self.name = tree, name

    def GetEntry(self, index):
        return self.tree.branch_bytes.get(self.name, {}).get(index, 40)


class FakeTree:
    def __init__(self, rows):
        self.rows = rows
        self.branches = {"variables", "eventweight", "sourceevent"}
        self.leaves = {
            "variables": FakeLeaf(10), "eventweight": FakeLeaf(1),
            "sourceevent": FakeLeaf(1, "Long64_t"),
        }
        self.branch_bytes, self.entry_bytes = {}, {}
        self.on_entry = None

    def InheritsFrom(self, name):
        return name == "TTree"

    def GetBranch(self, name):
        return FakeBranch(self, name) if name in self.branches else None

    def GetLeaf(self, name):
        return self.leaves.get(name)

    def GetEntries(self):
        return len(self.rows)

    def GetEntry(self, index):
        self.variables, self.eventweight, self.sourceevent = self.rows[index]
        if self.on_entry:
            self.on_entry(index)
        return self.entry_bytes.get(index, 80)


class FakeFile:
    def __init__(self, tree):
        self.tree = tree
        self.recovered = self.zombie = self.closed = False

    def IsZombie(self):
        return self.zombie

    def TestBit(self, bit):
        return self.recovered

    def Get(self, name):
        return self.tree if name == "Data2" else None

    def Close(self):
        self.closed = True


def features(selected=2.0):
    return [125., 60., -.5, 50., .5, 2.1, 2.9, 25., -.4, selected]


class HOMetadataReaderTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / "response.root"
        self.path.write_bytes(b"synthetic ROOT fixture")
        self.handle = FakeFile(FakeTree([(features(), [1.], [3])]))
        root = SimpleNamespace(
            gROOT=SimpleNamespace(SetBatch=mock.Mock()),
            TFile=SimpleNamespace(kRecovered=1, Open=lambda path: self.handle),
        )
        spec = importlib.util.spec_from_file_location("ho_ml_reader_fixture", READER_SOURCE)
        self.reader = importlib.util.module_from_spec(spec)
        with mock.patch.dict(sys.modules, {"ROOT": root}):
            spec.loader.exec_module(self.reader)

    def read(self):
        return self.reader.read_ho_ROOT_varfile(self.path)

    def test_full_population_and_exact_large_source_ids_are_preserved(self):
        large_id = 2**53 + 1
        sentinel = [-999.] * 9 + [0.]
        self.handle.tree.rows = [
            (sentinel, [-.25], [large_id]),
            (features(), [-.75], [large_id]),
            (features(), 2., 2**63 - 1),
            (features(), [0.], [0]),
        ]
        rows, weights, source_events, entries = self.read()
        self.assertEqual(rows[0], dict(zip(self.reader.FEATURE_NAMES, sentinel)))
        self.assertEqual(weights, [-.25, -.75, 2., 0.])
        self.assertEqual(source_events, [large_id, large_id, 2**63 - 1, 0])
        self.assertEqual(entries, [0, 1, 2, 3])
        self.assertTrue(all(isinstance(value, int) for value in source_events))
        self.assertEqual(set(rows[1]), set(self.reader.FEATURE_NAMES))
        self.assertTrue(self.handle.closed)

    def test_missing_source_branch_is_rejected_and_closed(self):
        self.handle.tree.branches.remove("sourceevent")
        with self.assertRaisesRegex(KeyError, "sourceevent"):
            self.read()
        self.assertTrue(self.handle.closed)

    def test_malformed_source_leaf_is_rejected(self):
        for leaf in (None, FakeLeaf(0, "Long64_t"), FakeLeaf(2, "Long64_t"),
                     FakeLeaf(1, "Double_t"), FakeLeaf(1, "ULong64_t"),
                     FakeLeaf(1, "Long64_t", counted=True)):
            with self.subTest(leaf=leaf):
                self.handle.closed = False
                self.handle.tree.leaves["sourceevent"] = leaf
                with self.assertRaisesRegex(ValueError, "sourceevent\\[1\\]/L"):
                    self.read()
                self.assertTrue(self.handle.closed)

    def test_invalid_source_values_are_not_coerced(self):
        for value in (-1, 2**63, 1.5, 1., float("nan"), float("inf"), "1", True):
            with self.subTest(value=value):
                self.handle.closed = False
                self.handle.tree.rows = [(features(), [1.], [value])]
                with self.assertRaisesRegex(ValueError, "sourceevent at entry 0"):
                    self.read()
                self.assertTrue(self.handle.closed)

    def test_each_required_branch_must_be_readable(self):
        for name in ("variables", "eventweight", "sourceevent"):
            for entry_bytes in (0, -1):
                with self.subTest(branch=name, bytes=entry_bytes):
                    self.handle.closed = False
                    self.handle.tree.branch_bytes = {name: {0: entry_bytes}}
                    with self.assertRaisesRegex(OSError, "Unreadable HO " + name):
                        self.read()
                    self.assertTrue(self.handle.closed)

    def test_unreadable_tree_entry_is_rejected(self):
        self.handle.tree.entry_bytes[0] = 0
        with self.assertRaisesRegex(OSError, "Unreadable HO response tree entry 0"):
            self.read()
        self.assertTrue(self.handle.closed)

    def test_unselected_nonfinite_rows_are_rejected(self):
        for bad_weight in (False, True):
            with self.subTest(weight=bad_weight):
                values, weight = features(0.), 1.
                if bad_weight:
                    weight = float("nan")
                else:
                    values[0] = float("nan")
                self.handle.tree.rows = [(values, [weight], [0])]
                with self.assertRaisesRegex(ValueError, "Non-finite HO response tree"):
                    self.read()
                self.assertTrue(self.handle.closed)

    def test_recovered_or_zombie_file_is_rejected_and_closed(self):
        for problem in ("recovered", "zombie"):
            with self.subTest(problem=problem):
                self.handle = FakeFile(FakeTree([(features(), [1.], [0])]))
                setattr(self.handle, problem, True)
                with self.assertRaises(OSError):
                    self.read()
                self.assertTrue(self.handle.closed)

    def test_change_during_read_is_rejected_and_closed(self):
        self.handle.tree.on_entry = lambda entry: self.path.write_bytes(b"modified fixture")
        with self.assertRaisesRegex(OSError, "changed during inspection"):
            self.read()
        self.assertTrue(self.handle.closed)

    def test_empty_tree_returns_empty_aligned_lists(self):
        self.handle.tree.rows = []
        self.assertEqual(self.read(), ([], [], [], []))
        self.assertTrue(self.handle.closed)

    def test_legacy_reader_does_not_require_source_metadata(self):
        self.handle.tree.branches.remove("sourceevent")
        self.handle.tree.leaves.pop("sourceevent")
        rows, labels, weights = self.reader.read_ROOT_varfile(self.path, 7, strict=True)
        self.assertEqual((rows, labels, weights), ([features()], [7], [1.]))

    def test_real_root_long64_source_ids_round_trip_when_available(self):
        try:
            import ROOT
        except ImportError:
            self.skipTest("PyROOT is not installed")
        from array import array

        source_ids = [0, 2**53 + 1, 2**63 - 1]
        output = ROOT.TFile(str(self.path), "RECREATE")
        tree = ROOT.TTree("Data2", "HO metadata fixture")
        values, weight, source = array("d", features()), array("d", [1.]), array("q", [0])
        tree.Branch("variables", values, "variables[10]/D")
        tree.Branch("eventweight", weight, "eventweight[1]/D")
        tree.Branch("sourceevent", source, "sourceevent[1]/L")
        for index, source_id in enumerate(source_ids):
            source[0], weight[0] = source_id, (-1.)**index
            tree.Fill()
        tree.Write()
        output.Close()

        spec = importlib.util.spec_from_file_location("real_ho_ml_reader_fixture", READER_SOURCE)
        reader = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reader)
        rows, weights, sources, entries = reader.read_ho_ROOT_varfile(self.path)
        self.assertEqual(sources, source_ids)
        self.assertEqual(weights, [1., -1., 1.])
        self.assertEqual(entries, [0, 1, 2])
        self.assertEqual(rows, [dict(zip(reader.FEATURE_NAMES, features()))] * 3)


if __name__ == "__main__":
    unittest.main()
