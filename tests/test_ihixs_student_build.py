"""Fresh-clone workflow fixtures; every compiler and ihixs call is mocked."""

import contextlib
import copy
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

from hgammagamma import ho_signal_normalization as norm
from hgammagamma import run_ihixs_normalization as wrapper


class StudentBuildTests(unittest.TestCase):
    def fixture(self, tmp):
        root = Path(tmp)
        source = root / "external/ihixs"
        originals = {
            "CMakeLists.txt": "# Synthetic fixture: never passed to a compiler.\n",
            "runcard/default.card": "# Synthetic published benchmark input.\noutput_filename = old.out\n",
            "src/tools/user_interface.cpp": (
                "char *name = new char[strlen(options[i].name.c_str())];\n"
                "long_options[N+1].name = 0;\n"
                "long_options[N+1].has_arg = 0;\n"
                "long_options[N+1].flag = 0;\n"
                "long_options[N+1].val = 0;\n"),
        }
        for name, (before, _) in wrapper.ADAPTER_REPLACEMENTS.items():
            originals[name] = (before + "\n") * wrapper.ADAPTER_COUNTS[name]
        for name, text in originals.items():
            path = source / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        hashes = {name: wrapper.file_sha256(source / name) for name in originals}
        library = root / "lhapdf/lib/libLHAPDF.synthetic"
        library.parent.mkdir(parents=True)
        library.write_text("Synthetic linked library: never loaded.\n")
        include = root / "lhapdf/include"
        include.mkdir(parents=True)
        args = wrapper.parse_args([
            "--work-dir", str(root / "normalization/student-fresh"),
            "--record", str(root / "normalization/student-rate.json"),
            "--ihixs-source", str(source), "--lhapdf-dir", str(root / "lhapdf"),
        ])
        settings = json.loads(wrapper.DEFAULT_SETTINGS.read_text())
        calls = []

        def fake_command(current, argv, *, cwd=None, log=None):
            """Write boundary artifacts without spawning any subprocess."""
            argv = [str(item) for item in argv]
            calls.append(argv)
            if argv[0] == current.cmake and "-S" in argv:
                directory = Path(argv[argv.index("-B") + 1])
                (directory / "CMakeCache.txt").write_text(
                    f"LHAPDF_INC_DIR:PATH={include}\nLHAPDF_LIB_NAMES:FILEPATH={library}\n")
            elif argv[0] == current.cmake and "--build" in argv:
                directory = Path(argv[argv.index("--build") + 1])
                (directory / "ihixs").write_text(
                    f"Synthetic {directory.name} executable: never executed.\n")
            elif argv[0] == current.cxx and "-o" in argv:
                Path(argv[argv.index("-o") + 1]).write_text("Synthetic probe: never executed.\n")
            elif len(argv) == 3 and argv[1:] == ["-i", "input.card"]:
                directory = Path(cwd)
                card = (directory / "input.card").read_text()
                if directory.name == "benchmark":
                    raw, key = 45.1816, "eftn3lo"
                else:
                    options = dict(line.split(" = ", 1) for line in card.splitlines()
                                   if " = " in line and not line.startswith("#"))
                    raw = 100. + .001 * int(options["pdf_member"])
                    key = "eftnnlo" if options["qcd_perturbative_order"] == "NNLO" else "eftn3lo"
                (directory / "ihixs.out").write_text(
                    f"{key} = {raw:.12g} [0.001]\nas_at_mz = 0.118\nas_at_mur = 0.125\n")
            else:
                raise AssertionError(f"unexpected external command: {argv}")
            if log is not None:
                Path(log).write_text("Synthetic command succeeded; no program was run.\n")
            return ""

        def fake_pdf_probe(current, point, expected_order=None):
            order = 3 if point["pdf_set"] == norm.PDF_SET else 2
            if expected_order is not None:
                self.assertEqual(expected_order, order)
            return {"version": "synthetic", "info_sha256": "b" * 64,
                    "member_sha256": "c" * 64, "OrderQCD": order, "AlphaS_OrderQCD": order,
                    "AlphaS_MZ": ".118", "alpha_s_mur": .125, "alpha_s_91_1876": .118,
                    "photon": True, "NumFlavors": "5", "NumMembers": "101",
                    "ErrorType": "replicas"}

        return args, settings, originals, hashes, calls, fake_command, fake_pdf_probe

    def test_fresh_build_repairs_both_parsers_and_only_adapts_production_physics(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings)
            parser_name = "src/tools/user_interface.cpp"
            for variant in ("upstream", "lhapdf"):
                parser = (args.work_dir / f"source-{variant}" / parser_name).read_text()
                self.assertIn("new char[strlen(options[i].name.c_str()) + 1]", parser)
                self.assertNotIn("long_options[N+1].", parser)
                for member in ("name", "has_arg", "flag", "val"):
                    self.assertIn(f"long_options[N].{member}", parser)
                for name, (before, after) in wrapper.ADAPTER_REPLACEMENTS.items():
                    text = (args.work_dir / f"source-{variant}" / name).read_text()
                    expected = (originals[name] if variant == "upstream"
                                else originals[name].replace(before, after))
                    self.assertEqual(text, expected)
            for name, original in originals.items():
                self.assertEqual((args.ihixs_source / name).read_text(), original)
                self.assertEqual(wrapper.file_sha256(args.ihixs_source / name), hashes[name])
            build = wrapper.load_build(args)
            self.assertEqual(build["source_files"], hashes)
            self.assertEqual(build["identity"]["source_sha256"], norm.canonical_sha256(hashes))
            self.assertEqual(build["identity"]["parser_patch_sha256"], wrapper.parser_repair_sha256())

    def test_verified_fixed_build_reuses_without_any_compiler_call(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, calls, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings)
            manifest = (args.work_dir / "build-manifest.json").read_bytes()
            self.assertTrue(calls)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=AssertionError("unexpected compiler call")), \
                    contextlib.redirect_stdout(io.StringIO()):
                wrapper.build(args, settings)
            self.assertEqual((args.work_dir / "build-manifest.json").read_bytes(), manifest)

    def test_legacy_build_is_preserved_instead_of_rewritten_by_a_normal_build(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings, parser_repair=False)
            before = {str(path.relative_to(args.work_dir)): path.read_bytes()
                      for path in args.work_dir.rglob("*") if path.is_file()}
            capture = io.StringIO()
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=AssertionError("legacy build was rebuilt")), \
                    contextlib.redirect_stdout(capture):
                wrapper.build(args, settings)
            after = {str(path.relative_to(args.work_dir)): path.read_bytes()
                     for path in args.work_dir.rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            self.assertIn("legacy", capture.getvalue().lower())
            self.assertNotIn("parser_patch_sha256", wrapper.load_build(args)["identity"])

    def test_changed_legacy_build_inputs_still_fail_closed(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings, parser_repair=False)
            changed = copy.copy(args)
            changed.cc = "different-compiler"
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=AssertionError("legacy build was rebuilt")):
                with self.assertRaisesRegex(ValueError, "build inputs changed"):
                    wrapper.build(changed, settings)

    def test_unfinished_legacy_calculation_refuses_before_loading_a_benchmark_or_pdf(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings, parser_repair=False)
            manifest = (args.work_dir / "build-manifest.json").read_bytes()
            with mock.patch.object(wrapper, "load_benchmark",
                                   side_effect=AssertionError("legacy calculation loaded a benchmark")), \
                    mock.patch.object(wrapper, "pdf_probe",
                                      side_effect=AssertionError("legacy calculation inspected a PDF")), \
                    mock.patch.object(wrapper, "command", side_effect=AssertionError("legacy executable was run")):
                with self.assertRaisesRegex(ValueError, "legacy ihixs build.*--repair-parser --resume"):
                    wrapper.calculate(args, settings)
            self.assertFalse((args.work_dir / "runs").exists())
            self.assertFalse(args.record.exists())
            self.assertEqual((args.work_dir / "build-manifest.json").read_bytes(), manifest)

    def test_legacy_benchmark_refuses_before_reading_a_card_or_running_a_probe(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, _ = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command):
                wrapper.build(args, settings, parser_repair=False)
            manifest = (args.work_dir / "build-manifest.json").read_bytes()
            # A missing card would raise FileNotFoundError if the safety guard
            # failed to reject the old executable before opening its inputs.
            (args.work_dir / "source-upstream/runcard/default.card").unlink()
            before = {str(path.relative_to(args.work_dir)): path.read_bytes()
                      for path in args.work_dir.rglob("*") if path.is_file()}
            with mock.patch.object(wrapper, "pdf_probe",
                                   side_effect=AssertionError("legacy benchmark inspected a PDF")), \
                    mock.patch.object(wrapper, "command", side_effect=AssertionError("legacy executable was run")):
                with self.assertRaisesRegex(ValueError, "legacy ihixs benchmark build.*--repair-parser --resume"):
                    wrapper.benchmark(args, settings)
            after = {str(path.relative_to(args.work_dir)): path.read_bytes()
                     for path in args.work_dir.rglob("*") if path.is_file()}
            self.assertEqual(before, after)
            self.assertEqual((args.work_dir / "build-manifest.json").read_bytes(), manifest)
            self.assertFalse((args.work_dir / "benchmark.json").exists())
            self.assertFalse((args.work_dir / "runs").exists())

    def test_fresh_fixed_workflow_installs_a_complete_record_without_legacy_recovery(self):
        with TemporaryDirectory() as tmp:
            args, settings, originals, hashes, _, command, probe = self.fixture(tmp)
            with mock.patch.object(wrapper, "source_inventory", return_value=(list(originals), hashes)), \
                    mock.patch.object(wrapper, "command", side_effect=command), \
                    mock.patch.object(wrapper, "pdf_probe", side_effect=probe), \
                    mock.patch.object(wrapper, "check_repaired_central",
                                      side_effect=AssertionError("fresh build required a legacy central result")), \
                    contextlib.redirect_stdout(io.StringIO()):
                wrapper.prepare(args, settings)
                wrapper.build(args, settings)
                wrapper.benchmark(args, settings)
                self.assertFalse((args.work_dir / "runs/central/complete.json").exists())
                wrapper.calculate(args, settings)
            record = norm.load_record(args.record)
            self.assertEqual(len(record["runs"]), 109)
            self.assertEqual({run["build_variant"] for run in record["runs"]}, {"lhapdf"})
            self.assertEqual(record["provenance"]["parser_patch_sha256"], wrapper.parser_repair_sha256())
            self.assertEqual(record["provenance"]["parser_patches"],
                             json.loads(json.dumps({"replacements": wrapper.PARSER_REPAIRS,
                                                   "counts": wrapper.PARSER_REPAIR_COUNTS})))
            self.assertEqual(record["validations"]["benchmark"]["parser_patch_sha256"],
                             record["provenance"]["parser_patch_sha256"])
            self.assertNotIn("parser_repair", record["provenance"])
            self.assertNotIn("parser_repair_central", record["validations"])
            self.assertFalse((args.work_dir / "parser-repair").exists())
            self.assertEqual({run["label"] for run in record["runs"][-2:]},
                             {"nnlo_native", "n3lo_nnlo_pdf"})


if __name__ == "__main__":
    unittest.main()
