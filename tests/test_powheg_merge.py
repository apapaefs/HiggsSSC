"""POWHEG merging preserves completed products and interrupted attempts."""
import contextlib
import io
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

from hgammagamma import run_powheg_hjminnlo as powheg


class PowhegMergeTests(unittest.TestCase):
    def seed(self, directory, name, weights, header="first"):
        events = "".join(
            f"<event>\n1 1 {weight} 125.0 0.0073 0.118\n"
            "25 1 0 0 0 0 0 0 0 125 125 0 9\n</event>\n"
            for weight in weights
        )
        preamble = (f'<LesHouchesEvents version="1.0">\n<header>{header}</header>\n'
                    "<init>\n2212 2212 20000 20000 0 0 336100 336100 -4 1\n1 0 1 1\n</init>\n")
        path = directory / name
        path.write_text(preamble + events + "</LesHouchesEvents>\n")
        return path, preamble, events

    def signature(self, path):
        stat = path.stat()
        return (path.read_bytes(), stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)

    def test_first_merge_preserves_header_signed_events_and_order(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            first, header, events_a = self.seed(directory, "pwgevents-0001.lhe", [2.0, -1.0])
            second, _, events_b = self.seed(directory, "pwgevents-0002.lhe", [3.0], header="second")
            output = directory / "merged.lhe"
            self.assertEqual(powheg.merge_lhe_files([first, second], output, False), 3)
            self.assertEqual(output.read_text(), header + events_a + events_b + "</LesHouchesEvents>\n")
            self.assertEqual(list(directory.glob(".merged.lhe.merge-*")), [])

    def test_identical_resume_preserves_inode_timestamps_and_contents(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            directory = Path(tmp)
            first, _, _ = self.seed(directory, "pwgevents-0001.lhe", [2.0, -1.0])
            output = directory / "merged.lhe"
            powheg.merge_lhe_files([first], output, False)
            os.utime(output, ns=(1_500_000_000_000_000_000, 1_500_000_000_000_000_000))
            before = self.signature(output)
            self.assertEqual(powheg.merge_lhe_files([first], output, False), 2)
            self.assertEqual(self.signature(output), before)
            self.assertEqual(list(directory.glob(".merged.lhe.merge-*")), [])

    def test_changed_seed_preserves_existing_merge_and_retains_candidate(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            first, _, _ = self.seed(directory, "pwgevents-0001.lhe", [2.0])
            output = directory / "merged.lhe"
            powheg.merge_lhe_files([first], output, False)
            before = self.signature(output)
            self.seed(directory, first.name, [-1.0])
            with self.assertRaisesRegex(SystemExit, "existing merged LHE differs.*preserved"):
                powheg.merge_lhe_files([first], output, False)
            self.assertEqual(self.signature(output), before)
            attempts = list(directory.glob(".merged.lhe.merge-*"))
            self.assertEqual(len(attempts), 1)
            self.assertEqual(attempts[0].read_bytes(), first.read_bytes())

    def test_incomplete_source_never_publishes_or_overwrites_a_merge(self):
        for existing in (False, True):
            for broken in ("<LesHouchesEvents>\n", "<LesHouchesEvents>\n<event>\n",
                           "<LesHouchesEvents>\n<event>\n</event>\n"):
                with self.subTest(existing=existing, broken=broken), TemporaryDirectory() as tmp, \
                        contextlib.redirect_stdout(io.StringIO()):
                    directory = Path(tmp)
                    first, _, _ = self.seed(directory, "pwgevents-0001.lhe", [2.0])
                    second = directory / "pwgevents-0002.lhe"
                    second.write_text(broken)
                    output = directory / "merged.lhe"
                    if existing:
                        output.write_text("saved partial merge from an earlier failed attempt")
                        before = self.signature(output)
                    with self.assertRaises(SystemExit):
                        powheg.merge_lhe_files([first, second], output, False)
                    if existing:
                        self.assertEqual(self.signature(output), before)
                    else:
                        self.assertFalse(output.exists())
                    attempts = list(directory.glob(".merged.lhe.merge-*"))
                    self.assertEqual(len(attempts), 1)
                    self.assertTrue(attempts[0].stat().st_size)

    def test_late_read_failure_preserves_existing_output_and_partial_attempt(self):
        with TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            directory = Path(tmp)
            first, _, _ = self.seed(directory, "pwgevents-0001.lhe", [2.0])
            output = directory / "merged.lhe"
            output.write_text("existing saved product")
            before = self.signature(output)
            with self.assertRaises(FileNotFoundError):
                powheg.merge_lhe_files([first, directory / "missing.lhe"], output, False)
            self.assertEqual(self.signature(output), before)
            self.assertEqual(len(list(directory.glob(".merged.lhe.merge-*"))), 1)

    def test_atomic_publication_does_not_replace_a_concurrent_destination(self):
        with TemporaryDirectory() as tmp:
            directory = Path(tmp)
            first, _, _ = self.seed(directory, "pwgevents-0001.lhe", [2.0])
            output = directory / "merged.lhe"
            link = os.link
            def competing_publication(source, destination):
                destination.write_text("another completed writer")
                return link(source, destination)
            with mock.patch.object(powheg.os, "link", side_effect=competing_publication):
                with self.assertRaisesRegex(SystemExit, "existing merged LHE differs"):
                    powheg.merge_lhe_files([first], output, False)
            self.assertEqual(output.read_text(), "another completed writer")
            self.assertEqual(len(list(directory.glob(".merged.lhe.merge-*"))), 1)


if __name__ == "__main__":
    unittest.main()
