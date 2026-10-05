"""Mocked ihixs parser-repair checks; no executable or integration is run."""

import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest import mock

from hgammagamma import ho_signal_normalization as norm
from hgammagamma import run_ihixs_normalization as wrapper


class ParserSourceRepairTests(unittest.TestCase):
    def source(self):
        name = "src/tools/user_interface.cpp"
        path = wrapper.REPO_ROOT / "external/ihixs" / name
        if not path.is_file():
            self.skipTest("initialize the pinned external/ihixs submodule for source-anchor checks")
        return name, path.read_text()

    def test_repaired_parser_allocates_the_terminator_and_bounds_the_sentinel(self):
        name, original = self.source()
        repaired = wrapper.patched_source(name, original, "upstream", parser_repair=True)
        self.assertRegex(repaired, r"new char\[strlen\(options\[i\]\.name\.c_str\(\)\)\s*\+\s*1\]")
        self.assertNotIn("long_options[N+1].", repaired)
        for member in ("name", "has_arg", "flag", "val"):
            self.assertIn(f"long_options[N].{member}", repaired)
        self.assertEqual(wrapper.patched_source(name, original, "lhapdf", parser_repair=True), repaired)

    def test_legacy_source_and_adapter_remain_unchanged_without_repair(self):
        name, original = self.source()
        self.assertEqual(wrapper.patched_source(name, original, "upstream"), original)
        self.assertEqual(wrapper.patched_source(name, original, "lhapdf"), original)
        adapter_name, (before, after) = next(iter(wrapper.ADAPTER_REPLACEMENTS.items()))
        self.assertEqual(wrapper.patched_source(adapter_name, before, "upstream"), before)
        self.assertEqual(wrapper.patched_source(adapter_name, before, "lhapdf"), after)

    def test_each_parser_anchor_is_checked_before_replacing_source(self):
        name, original = self.source()
        for before, _ in wrapper.PARSER_REPAIRS[name]:
            with self.subTest(anchor=before):
                changed = original.replace(before, "synthetic changed upstream anchor", 1)
                with self.assertRaisesRegex(ValueError, "anchor"):
                    wrapper.patched_source(name, changed, "lhapdf", parser_repair=True)

    def test_repair_flag_requires_resume_and_an_execution_stage(self):
        for stage in ("build", "benchmark", "calculate"):
            args = wrapper.parse_args(["--stage", stage, "--resume", "--repair-parser"])
            self.assertTrue(args.repair_parser)
            self.assertTrue(args.resume)
        for argv in (["--repair-parser"],
                     ["--stage", "prepare", "--resume", "--repair-parser"],
                     ["--stage", "build", "--repair-parser"],
                     ["--stage", "benchmark", "--repair-parser"],
                     ["--stage", "calculate", "--repair-parser"]):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    wrapper.parse_args(argv)
                self.assertEqual(raised.exception.code, 2)


class ParserRepairBuildTests(unittest.TestCase):
    def fixture(self, tmp):
        args = SimpleNamespace(work_dir=Path(tmp) / "work", resume=True, repair_parser=True)
        repair_args = wrapper.parser_repair_args(args)
        library = Path(tmp) / "libLHAPDF.synthetic"
        library.write_text("synthetic linked library")
        source_files = {"ihixs.cpp": "a" * 64, "src/tools/user_interface.cpp": "b" * 64}
        identity = {"ihixs_commit": norm.IHIXS_COMMIT,
                    "source_sha256": norm.canonical_sha256(source_files),
                    "adapter_sha256": norm.canonical_sha256({
                        "replacements": wrapper.ADAPTER_REPLACEMENTS, "counts": wrapper.ADAPTER_COUNTS}),
                    "cc": "cc", "cxx": "c++", "lhapdf_dir": str(Path(tmp)),
                    "cuba_dir": str(Path(tmp) / "cuba"), "boost_dir": "None",
                    "probe_source_sha256": hashlib.sha256(wrapper.PDF_PROBE.encode()).hexdigest()}
        builds = []
        for current, suffix in ((args, "original"), (repair_args, "repaired")):
            hashes = {}
            for variant in ("upstream", "lhapdf"):
                binary = current.work_dir / f"build-{variant}" / "ihixs"
                binary.parent.mkdir(parents=True)
                binary.write_text(f"synthetic {suffix} {variant} executable, never executed")
                hashes[variant] = wrapper.file_sha256(binary)
            probe = current.work_dir / "lhapdf_probe"
            probe.write_text(f"synthetic {suffix} probe, never executed")
            build_identity = copy.deepcopy(identity)
            if suffix == "repaired":
                build_identity["parser_patch_sha256"] = wrapper.parser_repair_sha256()
            data = {"schema_version": 1, "identity": build_identity, "source_files": source_files,
                    "executable_sha256": hashes, "probe_sha256": wrapper.file_sha256(probe),
                    "lhapdf_library_path": str(library),
                    "lhapdf_library_sha256": wrapper.file_sha256(library)}
            (current.work_dir / "build-manifest.json").write_text(json.dumps(data))
            builds.append(data)
        return args, repair_args, builds[0], builds[1]

    def test_repair_workspace_is_separate_and_original_args_are_preserved(self):
        with TemporaryDirectory() as tmp:
            args, repair_args, _, _ = self.fixture(tmp)
            self.assertIsNot(args, repair_args)
            self.assertEqual(args.work_dir, Path(tmp) / "work")
            self.assertTrue(args.repair_parser)
            self.assertEqual(repair_args.work_dir, args.work_dir / "parser-repair")
            self.assertFalse(repair_args.repair_parser)

    def test_verified_repair_matches_original_source_adapter_and_lhapdf_library(self):
        with TemporaryDirectory() as tmp:
            args, repair_args, base, repair = self.fixture(tmp)
            with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")):
                selected_args, selected_build = wrapper.load_parser_repair(args, base)
            self.assertEqual(selected_args.work_dir, repair_args.work_dir)
            self.assertEqual(selected_build, repair)
            self.assertEqual(selected_build["source_files"], base["source_files"])

    def test_mismatched_source_adapter_patch_and_library_are_rejected(self):
        mutations = (
            lambda data: data["identity"].update(source_sha256="f" * 64),
            lambda data: data["identity"].update(adapter_sha256="f" * 64),
            lambda data: data["identity"].update(parser_patch_sha256="f" * 64),
            lambda data: data["source_files"].update({"ihixs.cpp": "f" * 64}),
            lambda data: data.update(lhapdf_library_sha256="f" * 64),
        )
        for mutate in mutations:
            with self.subTest(mutation=mutate), TemporaryDirectory() as tmp:
                args, repair_args, base, repair = self.fixture(tmp)
                bad = copy.deepcopy(repair)
                mutate(bad)
                (repair_args.work_dir / "build-manifest.json").write_text(json.dumps(bad))
                with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")):
                    with self.assertRaises(ValueError):
                        wrapper.load_parser_repair(args, base)

    def test_repaired_binary_hash_cannot_change_after_build(self):
        with TemporaryDirectory() as tmp:
            args, repair_args, base, _ = self.fixture(tmp)
            (repair_args.work_dir / "build-lhapdf/ihixs").write_text("modified synthetic executable")
            with self.assertRaisesRegex(ValueError, "executable changed"):
                wrapper.load_parser_repair(args, base)

    def benchmark_fixture(self, args, build):
        directory = args.work_dir / "runs/benchmark"
        directory.mkdir(parents=True)
        (directory / "input.card").write_text("synthetic published benchmark card")
        (directory / "ihixs.out").write_text("synthetic published benchmark output")
        data = {"passed": True, "ihixs_commit": norm.IHIXS_COMMIT,
                "input_sha256": wrapper.file_sha256(directory / "input.card"),
                "output_sha256": wrapper.file_sha256(directory / "ihixs.out"),
                "executable_sha256": build["executable_sha256"]["upstream"],
                "production_executable_sha256": build["executable_sha256"]["lhapdf"]}
        if "parser_patch_sha256" in build["identity"]:
            data["parser_patch_sha256"] = build["identity"]["parser_patch_sha256"]
        (args.work_dir / "benchmark.json").write_text(json.dumps(data))
        return data

    def test_original_and_repaired_benchmarks_bind_to_their_own_artifacts_and_builds(self):
        with TemporaryDirectory() as tmp:
            args, repair_args, base, repair = self.fixture(tmp)
            for current, build in ((args, base), (repair_args, repair)):
                with self.subTest(work_dir=current.work_dir):
                    expected = self.benchmark_fixture(current, build)
                    self.assertEqual(wrapper.load_benchmark(current, build), expected)

    def test_changed_benchmark_output_binary_or_parser_patch_cannot_be_reused(self):
        for key in ("output_sha256", "executable_sha256", "production_executable_sha256",
                    "parser_patch_sha256"):
            with self.subTest(key=key), TemporaryDirectory() as tmp:
                _, repair_args, _, repair = self.fixture(tmp)
                data = self.benchmark_fixture(repair_args, repair)
                data[key] = "f" * 64
                (repair_args.work_dir / "benchmark.json").write_text(json.dumps(data))
                with self.assertRaises(ValueError):
                    wrapper.load_benchmark(repair_args, repair)


class ParserRepairResumeTests(unittest.TestCase):
    def fixture(self, tmp, *, refine_failed=False):
        work = Path(tmp) / "work"
        args = SimpleNamespace(work_dir=work, resume=True, refine_failed=refine_failed,
                               repair_parser=True)
        repair_args = wrapper.parser_repair_args(args)
        original = work / "build-lhapdf/ihixs"
        repaired = repair_args.work_dir / "build-lhapdf/ihixs"
        for executable, text in ((original, "original"), (repaired, "parser repaired")):
            executable.parent.mkdir(parents=True)
            executable.write_text(f"synthetic {text} executable, never executed")
        repair_build = {"executable_sha256": {"lhapdf": wrapper.file_sha256(repaired)},
                        "probe_sha256": "e" * 64,
                        "identity": {"parser_patch_sha256": wrapper.parser_repair_sha256()}}
        settings = json.loads(wrapper.DEFAULT_SETTINGS.read_text())
        settings["integration"].update(epsrel=5e-5, epsabs=0., mineval=200000,
                                       maxeval=200000000, nstart=40000, nincrease=4000)
        point = {"label": "nnlo_native", "mur_gev": 62.5, "muf_gev": 62.5,
                 "pdf_member": 0, "pdf_set": norm.NATIVE_PDF_SET, "qcd_order": "NNLO"}
        identity = {"settings_sha256": norm.canonical_sha256(settings), "pdf_info_sha256": "b" * 64,
                    "pdf_member_sha256": "c" * 64, "probe_sha256": "d" * 64, "variant": "lhapdf"}
        pdf = {"alpha_s_mur": .125, "alpha_s_91_1876": .118}
        return args, settings, point, original, repaired, identity, pdf, (repair_args, repair_build)

    @staticmethod
    def output(error=.01):
        return (f"eftnnlo = 1.00000000e+02 [{error}]\n"
                "as_at_mz = 1.18000000e-01\nas_at_mur = 1.25000000e-01\n")

    def complete(self, args, label, card, executable, identity, *, error=.01):
        directory = args.work_dir / "runs" / label
        directory.mkdir(parents=True)
        (directory / "input.card").write_text(card)
        (directory / "ihixs.out").write_text(self.output(error))
        (directory / "ihixs.log").write_text("synthetic successful exit")
        payload = {"identity": {**identity, "input_sha256": hashlib.sha256(card.encode()).hexdigest(),
                                "executable_sha256": wrapper.file_sha256(executable)},
                   "output_sha256": wrapper.file_sha256(directory / "ihixs.out")}
        (directory / "complete.json").write_text(json.dumps(payload))
        return directory, payload

    @staticmethod
    def snapshot(root):
        return {str(path.relative_to(root)): path.read_bytes()
                for path in root.rglob("*") if path.is_file()}

    def test_legacy_completed_point_is_reused_with_its_original_binary(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, original, _, identity, pdf, repair = self.fixture(tmp)
            directory, saved = self.complete(args, point["label"], wrapper.production_card(settings, point),
                                             original, identity)
            before = self.snapshot(args.work_dir)
            with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")), \
                    mock.patch.object(wrapper, "run_card", wraps=wrapper.run_card) as runner:
                selected, completed, _ = wrapper.run_precision_checked(
                    args, settings, point, original, identity, pdf, repair_context=repair)
            self.assertEqual(selected, directory)
            self.assertEqual(completed, saved)
            self.assertEqual(runner.call_args.args[3], original)
            self.assertEqual(runner.call_args.args[4], identity)
            self.assertEqual(self.snapshot(args.work_dir), before)

    def test_mixed_legacy_base_and_repaired_precision_cache_choose_their_own_binaries(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, original, repaired, identity, pdf, repair = self.fixture(tmp, refine_failed=True)
            self.complete(args, point["label"], wrapper.production_card(settings, point),
                          original, identity, error=.1)
            refined = copy.deepcopy(settings)
            refined["integration"]["epsrel"] = 1e-6
            refined["integration"]["epsabs"] *= 1e-6 / settings["integration"]["epsrel"]
            fixed_identity = {**identity, "settings_sha256": norm.canonical_sha256(refined),
                              "variant": "lhapdf-parser-repair"}
            retry_label = point["label"] + "__precision_1"
            expected, saved = self.complete(args, retry_label, wrapper.production_card(refined, point),
                                            repaired, fixed_identity)
            before = self.snapshot(args.work_dir)
            with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")), \
                    mock.patch.object(wrapper, "run_card", wraps=wrapper.run_card) as runner, \
                    contextlib.redirect_stdout(io.StringIO()):
                selected, completed, result = wrapper.run_precision_checked(
                    args, settings, point, original, identity, pdf, repair_context=repair)
            self.assertEqual(selected, expected)
            self.assertEqual(completed, saved)
            self.assertEqual([call.args[3] for call in runner.call_args_list], [original, repaired])
            self.assertEqual(result["precision_refinements"], 1)
            self.assertEqual(self.snapshot(args.work_dir), before)

    def test_aborted_output_is_archived_and_requires_a_successful_repaired_run(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, original, repaired, identity, pdf, repair = self.fixture(tmp)
            partial = args.work_dir / "runs" / point["label"]
            partial.mkdir(parents=True)
            (partial / "input.card").write_text(wrapper.production_card(settings, point))
            (partial / "ihixs.out").write_text(self.output())
            (partial / "ihixs.log").write_text("corrupted size vs. prev_size")
            aborted = {path.name: path.read_bytes() for path in partial.iterdir()}

            def successful_command(current, argv, *, cwd=None, log=None):
                self.assertEqual(current.work_dir, args.work_dir)
                self.assertEqual(argv[0], repaired)
                self.assertFalse((cwd / "complete.json").exists())
                self.assertFalse((cwd / "ihixs.out").exists())
                (cwd / "ihixs.out").write_text(self.output(.005))
                log.write_text("synthetic repaired successful exit")
                return ""

            with mock.patch.object(wrapper, "command", side_effect=successful_command) as command:
                directory, completed, _ = wrapper.run_precision_checked(
                    args, settings, point, original, identity, pdf, repair_context=repair)
            command.assert_called_once()
            backups = list((args.work_dir / "resume-backups").iterdir())
            self.assertEqual(len(backups), 1)
            self.assertEqual({path.name: path.read_bytes() for path in backups[0].iterdir()}, aborted)
            self.assertEqual(directory, partial)
            self.assertEqual(completed["identity"]["variant"], "lhapdf-parser-repair")
            self.assertEqual(completed["identity"]["probe_sha256"], identity["probe_sha256"])
            self.assertEqual(completed["identity"]["executable_sha256"], wrapper.file_sha256(repaired))

    def test_repaired_abort_cannot_publish_a_completion_or_trigger_precision_retry(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, original, repaired, identity, pdf, repair = self.fixture(tmp, refine_failed=True)

            def aborted_command(current, argv, *, cwd=None, log=None):
                self.assertEqual(argv[0], repaired)
                (cwd / "ihixs.out").write_text(self.output())
                log.write_text("synthetic SIGABRT after writing output")
                raise RuntimeError("command failed (-6)")

            with mock.patch.object(wrapper, "command", side_effect=aborted_command) as command:
                with self.assertRaisesRegex(RuntimeError, r"command failed \(-6\)"):
                    wrapper.run_precision_checked(args, settings, point, original, identity, pdf,
                                                  repair_context=repair)
            command.assert_called_once()
            self.assertTrue((args.work_dir / "runs/nnlo_native/ihixs.out").is_file())
            self.assertFalse(list((args.work_dir / "runs").rglob("complete.json")))
            self.assertFalse((args.work_dir / "runs/nnlo_native__precision_1").exists())

    def test_unknown_variant_and_changed_cached_binary_fail_closed(self):
        for mutation in ({"variant": "unknown-build"}, {"executable_sha256": "f" * 64}):
            with self.subTest(mutation=mutation), TemporaryDirectory() as tmp:
                args, settings, point, original, _, identity, pdf, repair = self.fixture(tmp)
                directory, saved = self.complete(args, point["label"], wrapper.production_card(settings, point),
                                                 original, identity)
                saved["identity"].update(mutation)
                (directory / "complete.json").write_text(json.dumps(saved))
                before = self.snapshot(args.work_dir)
                with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")):
                    with self.assertRaises(ValueError):
                        wrapper.run_precision_checked(args, settings, point, original, identity, pdf,
                                                      repair_context=repair)
                self.assertEqual(self.snapshot(args.work_dir), before)


if __name__ == "__main__":
    unittest.main()
