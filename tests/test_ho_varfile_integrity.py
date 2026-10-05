"""HO response-tree integrity checks with fake ROOT objects; no PyROOT required."""

import importlib.util
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock


READER_SOURCE = Path(__file__).resolve().parents[1] / "read_root_varfiles.py"


class FakeLeaf:
    def __init__(self, length, counted=False):
        self.length = length
        self.counted = counted

    def GetLenStatic(self):
        return self.length

    def GetLeafCount(self):
        return object() if self.counted else None


class FakeBranch:
    def __init__(self, tree, name):
        self.tree = tree
        self.name = name

    def GetEntry(self, index):
        entry_bytes = self.tree.branch_bytes.get(self.name)
        return 40 if entry_bytes is None else entry_bytes[index]


class FakeTree:
    def __init__(self, rows, *, is_tree=True, entry_bytes=None, on_entry=None):
        self.rows = rows
        self.is_tree = is_tree
        self.entry_bytes = entry_bytes
        self.on_entry = on_entry
        self.branches = {"variables", "eventweight"}
        self.branch_bytes = {}
        self.leaves = {"variables": FakeLeaf(10), "eventweight": FakeLeaf(1)}

    def InheritsFrom(self, name):
        return self.is_tree and name == "TTree"

    def GetBranch(self, name):
        return FakeBranch(self, name) if name in self.branches else None

    def GetLeaf(self, name):
        return self.leaves.get(name)

    def GetEntries(self):
        return len(self.rows)

    def GetEntry(self, index):
        self.variables, self.eventweight = self.rows[index]
        if self.on_entry is not None:
            self.on_entry(index)
        return 80 if self.entry_bytes is None else self.entry_bytes[index]


class FakeFile:
    def __init__(self, tree, *, recovered=False, zombie=False):
        self.tree = tree
        self.recovered = recovered
        self.zombie = zombie
        self.closed = False

    def IsZombie(self):
        return self.zombie

    def TestBit(self, bit):
        return self.recovered

    def Get(self, name):
        return self.tree if name == "Data2" else None

    def Close(self):
        self.closed = True


def finite_features(selected=2.0):
    return [125.0, 60.0, -1.2, 50.0, 0.8, 2.1, 2.9, 25.0, -0.4, selected]


class HOResponseTreeIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "response.root"
        self.path.write_bytes(b"synthetic ROOT fixture")
        self.handle = FakeFile(FakeTree([(finite_features(), [1.0])]))
        fake_root = SimpleNamespace(
            gROOT=SimpleNamespace(SetBatch=mock.Mock()),
            TFile=SimpleNamespace(kRecovered=1, Open=lambda path: self.handle),
        )
        # Load a separate module with a local fake ROOT dependency, leaving any
        # real ROOT or previously imported reader module untouched.
        spec = importlib.util.spec_from_file_location("ho_varfile_reader_fixture", READER_SOURCE)
        self.reader = importlib.util.module_from_spec(spec)
        with mock.patch.dict(sys.modules, {"ROOT": fake_root}):
            spec.loader.exec_module(self.reader)

    def read(self, **kwargs):
        return self.reader.read_ROOT_varfile(self.path, sample_id=7, **kwargs)

    def test_strict_preserves_signed_weights_and_finite_unselected_sentinels(self):
        sentinel = [-999.0] * 9 + [0.0]
        selected = finite_features()
        self.handle.tree = FakeTree([(sentinel, [-0.25]), (selected, 1.5)])
        rows, labels, weights = self.read(strict=True, xsec=10.0)
        self.assertEqual(rows, [sentinel, selected])
        self.assertEqual(labels, [7, 7])
        self.assertEqual(weights, [-2.5, 15.0])
        self.assertTrue(self.handle.closed)

    def test_named_reader_forwards_strict_and_preserves_feature_names(self):
        rows, weights = self.reader.read_named_ROOT_varfile(self.path, strict=True)
        self.assertEqual(rows, [dict(zip(self.reader.FEATURE_NAMES, finite_features()))])
        self.assertEqual(weights, [1.0])
        self.handle = FakeFile(FakeTree([(finite_features(), [1.0])]), recovered=True)
        with self.assertRaisesRegex(OSError, "Recovered ROOT"):
            self.reader.read_named_ROOT_varfile(self.path, strict=True)

    def test_recovered_file_is_rejected_and_closed(self):
        self.handle.recovered = True
        with self.assertRaisesRegex(OSError, "Recovered ROOT"):
            self.read(strict=True)
        self.assertTrue(self.handle.closed)

    def test_object_with_tree_like_methods_must_be_a_ttree(self):
        self.handle.tree.is_tree = False
        with self.assertRaisesRegex(TypeError, "not a TTree"):
            self.read(strict=True)
        self.assertTrue(self.handle.closed)

    def test_missing_required_branch_is_rejected(self):
        for name in ("variables", "eventweight"):
            with self.subTest(branch=name):
                self.handle = FakeFile(FakeTree([(finite_features(), [1.0])]))
                self.handle.tree.branches.remove(name)
                with self.assertRaisesRegex(KeyError, name):
                    self.read(strict=True)
                self.assertTrue(self.handle.closed)

    def test_missing_wrong_length_or_counted_leaves_are_rejected(self):
        for name, length in (("variables", 10), ("eventweight", 1)):
            for leaf in (None, FakeLeaf(length - 1), FakeLeaf(length + 1),
                         FakeLeaf(length, counted=True)):
                with self.subTest(branch=name, leaf=leaf):
                    self.handle = FakeFile(FakeTree([(finite_features(), [1.0])]))
                    self.handle.tree.leaves[name] = leaf
                    with self.assertRaisesRegex(ValueError, "invalid " + name + " branch length"):
                        self.read(strict=True)
                    self.assertTrue(self.handle.closed)

    def test_unreadable_entries_are_rejected_and_file_is_closed(self):
        for entry_bytes in (0, -1):
            with self.subTest(entry_bytes=entry_bytes):
                self.handle = FakeFile(FakeTree([(finite_features(), [1.0])],
                                                entry_bytes=[entry_bytes]))
                with self.assertRaisesRegex(OSError, "Unreadable HO response tree entry 0"):
                    self.read(strict=True)
                self.assertTrue(self.handle.closed)

    def test_nonfinite_features_and_weights_are_rejected(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            for kind in ("feature", "weight"):
                with self.subTest(kind=kind, value=value):
                    features = finite_features()
                    weight = [1.0]
                    if kind == "feature":
                        features[0] = value
                    else:
                        weight[0] = value
                    self.handle = FakeFile(FakeTree([(features, weight)]))
                    with self.assertRaisesRegex(ValueError, "Non-finite HO response tree entry 0"):
                        self.read(strict=True)
                    self.assertTrue(self.handle.closed)

    def test_unreadable_required_branch_is_rejected_when_total_entry_is_readable(self):
        for name in ("variables", "eventweight"):
            for entry_bytes in (0, -1):
                with self.subTest(branch=name, entry_bytes=entry_bytes):
                    self.handle = FakeFile(FakeTree([(finite_features(), [1.0])]))
                    self.handle.tree.branch_bytes[name] = [entry_bytes]
                    with self.assertRaisesRegex(OSError, "Unreadable.*" + name):
                        self.read(strict=True)
                    self.assertTrue(self.handle.closed)

    def test_nonfinite_unselected_row_cannot_hide_behind_selection_filter(self):
        features = finite_features(selected=0.0)
        features[0] = float("nan")
        self.handle.tree = FakeTree([(features, [1.0])])
        with self.assertRaisesRegex(ValueError, "Non-finite HO response tree entry"):
            self.read(strict=True, selected_only=True)

    def test_changed_file_is_rejected_even_if_all_rows_are_finite(self):
        def change_file(index):
            with self.path.open("ab") as handle:
                handle.write(b"changed during inspection")

        self.handle.tree = FakeTree([(finite_features(), [1.0])], on_entry=change_file)
        with self.assertRaisesRegex(OSError, "changed during inspection"):
            self.read(strict=True)
        self.assertTrue(self.handle.closed)

    def test_legacy_default_retains_permissive_schema_and_entry_handling(self):
        self.handle.recovered = True
        self.handle.tree = FakeTree([(finite_features(), [-2.0])],
                                    is_tree=False, entry_bytes=[0])
        self.handle.tree.leaves.clear()
        rows, labels, weights = self.read()
        self.assertEqual(rows, [finite_features()])
        self.assertEqual(labels, [7])
        self.assertEqual(weights, [-2.0])

    def test_legacy_default_skips_nonfinite_rows_and_keeps_finite_signed_rows(self):
        bad_features = finite_features()
        bad_features[0] = float("nan")
        self.handle.tree = FakeTree([
            (bad_features, [1.0]),
            (finite_features(), [float("inf")]),
            (finite_features(), [-0.5]),
        ])
        rows, labels, weights = self.read()
        self.assertEqual(rows, [finite_features()])
        self.assertEqual(labels, [7])
        self.assertEqual(weights, [-0.5])


if __name__ == "__main__":
    unittest.main()
