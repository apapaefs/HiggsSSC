"""HO rate integration tests; generator and detector commands are mocked."""

import contextlib
import io
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hgammagamma import run_gammagamma_ho_campaign as ho


def signal_lhe_text():
    events = "".join(
        f"<event>\n1 1 {weight} 100 .007 .118\n"
        "25 1 1 2 0 0 0 0 0 125 125 0 9\n</event>\n"
        for weight in (10., -2., 8.)
    )
    return (
        '<LesHouchesEvents version="3.0">\n<header>\n'
        "minnlo 1\nhmass 125\nlhans1 336100\nlhans2 336100\n"
        "alphas_from_pdf 1\nuse_NNLOPS_pdfs 1\n</header>\n"
        "<init>\n2212 2212 20000 20000 -1 -1 336100 336100 -4 1\n"
        "5.333333333333333 .2 10 1\n</init>\n"
        f"{events}</LesHouchesEvents>\n"
    )


RECORD = {
    "higgs_mass_gev": 125., "native_parameters": {"hmass": 125.},
    "fingerprint": "a" * 64,
}


def sidecar_fixture(path, manifest, native_lhe, weight_scale):
    # Full record validation is covered by the shared normalization module's
    # tests. This fixture isolates the runner's orchestration and migration.
    return {
        "sample": manifest["sample"], "run_tag": manifest["run_tag"],
        "normalization_kind": "ihixs_n3lo", "cross_section_pb": 400.,
        "cross_section_error_pb": .05, "weight_scale": weight_scale,
        "native_cross_section_pb": native_lhe["cross_section_pb"],
        "native_cross_section_error_pb": native_lhe["cross_section_error_pb"],
        "native_lhe": native_lhe, "ihixs_record_sha256": RECORD["fingerprint"],
        "ihixs": RECORD,
    }


class HoNormalizationCampaignTests(unittest.TestCase):
    def args(self, directory, *extra):
        return ho.parse_args([
            "--output-dir", str(directory), "--no-herwig-module", *extra,
        ])

    def completed_campaign(self, root):
        args = self.args(root, "--run-tag", "stored_run", "--nevents", "3", "--higgs-br", ".003")
        sample = ho.SAMPLES[0]
        directory = ho.sample_dir(args, sample)
        directory.mkdir(parents=True)
        args.signal_lhe = Path(root) / "signal.lhe"
        args.signal_lhe.write_text(signal_lhe_text())
        event_path = directory / "herwig/events/input.root"
        event_path.parent.mkdir(parents=True)
        event_path.write_bytes(b"existing HwSim events")
        input_list = directory / f"{sample.name}_hwsim_roots.input"
        input_list.write_text(f"{event_path}\n")
        prefix = directory / f"{sample.name}_hwsim_roots-{args.run_tag}"
        summary = {
            "input": str(input_list),
            "events_read": "3", "weight_scale": ".003", "sum_weight": ".048",
            "sum_abs_weight": ".06", "sum_weight_squared": ".001512",
            "sum_tree_weight": ".048", "sum_diphoton_weight": ".027",
        }
        Path(str(prefix) + ".dat").write_text("".join(f"{key} {value}\n" for key, value in summary.items()))
        Path(str(prefix) + ".top").write_text("existing signed histograms\n")
        Path(str(prefix) + "_var.root").write_bytes(b"existing analysis tree")
        manifest = {
            **ho.signal_campaign_identity(args), "fingerprint": "unchanged-generation-fingerprint",
            "completed": ["generate", "shower", "analyze"], "nevents_requested": 3,
            "weight_scale": .003, "lhe_file": str(args.signal_lhe),
            "lhe": ho.inspect_lhe(args.signal_lhe, sample, args),
            "root_files": ho.root_inventory(directory), "analysis": summary,
        }
        (directory / "campaign.json").write_text(json.dumps(manifest))
        (directory / "normalization-stored_run.json").write_text(json.dumps({
            "sample": sample.name, "run_tag": args.run_tag,
            "cross_section_pb": manifest["lhe"]["cross_section_pb"], "weight_scale": .003,
        }))
        return args, directory, manifest

    def test_invalid_record_blocks_all_before_building_or_generating(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp, "--stage", "all", "--signal-lhe", str(Path(tmp) / "signal.lhe"))
            with patch.object(ho.signal_normalization, "load_signal_record", side_effect=ValueError("not calculated")), \
                    patch.object(ho, "command") as command, patch.object(ho, "prepare") as prepare:
                with self.assertRaisesRegex(ValueError, "validated --signal-normalization"):
                    ho.run_campaign(args)
                command.assert_not_called()
                prepare.assert_not_called()

    def test_dryrun_prints_pending_normalization_without_writes(self):
        with TemporaryDirectory() as tmp:
            output = Path(tmp) / "not-created"
            args = self.args(output, "--stage", "analyze", "--run-samples", "signal_gg_h_aa",
                             "--signal-lhe", str(Path(tmp) / "signal.lhe"), "--dry-run")
            buffer = io.StringIO()
            with patch.object(ho.signal_normalization, "load_signal_record", side_effect=ValueError("not calculated")), \
                    contextlib.redirect_stdout(buffer):
                self.assertEqual(ho.run_campaign(args), 0)
            self.assertIn("pending HO signal normalization", buffer.getvalue())
            self.assertFalse(output.exists())

    def test_normalize_infers_metadata_and_preserves_all_event_products(self):
        with TemporaryDirectory() as tmp:
            original_args, directory, old_manifest = self.completed_campaign(tmp)
            products = [path for path in directory.rglob("*") if path.is_file()
                        and path.name != "campaign.json" and not path.name.startswith("normalization-")]
            before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in products}
            args = self.args(tmp, "--stage", "normalize", "--ebeam", "6500", "--higgs-br", ".004")
            with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                    patch.object(ho.signal_normalization, "signal_sidecar", side_effect=sidecar_fixture), \
                    patch.object(ho, "command") as command, patch.object(ho, "prepare") as prepare:
                self.assertEqual(ho.run_campaign(args), 0)
                command.assert_not_called()
                prepare.assert_not_called()
            sidecar = json.loads((directory / "normalization-stored_run.json").read_text())
            self.assertEqual(sidecar["cross_section_pb"], 400.)
            self.assertEqual(sidecar["weight_scale"], .003)
            self.assertEqual(sidecar["native_cross_section_pb"], old_manifest["lhe"]["cross_section_pb"])
            updated = json.loads((directory / "campaign.json").read_text())
            self.assertEqual(updated["fingerprint"], old_manifest["fingerprint"])
            self.assertEqual(updated["completed"], old_manifest["completed"])
            self.assertEqual(updated["analysis"], old_manifest["analysis"])
            self.assertEqual(updated["rate_normalization"]["record_fingerprint"], RECORD["fingerprint"])
            for path, expected in before.items():
                self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), expected)
            self.assertFalse((directory / "normalization-ho_run_01.json").exists())
            self.assertEqual(original_args.run_tag, "stored_run")

    def test_normalize_rejects_changed_lhe_root_and_branching_weights(self):
        for changed in ("lhe", "root", "branching", "closure"):
            with self.subTest(changed=changed), TemporaryDirectory() as tmp:
                original_args, directory, _ = self.completed_campaign(tmp)
                if changed == "lhe":
                    original_args.signal_lhe.write_text(signal_lhe_text().replace("1 1 10.0", "1 1 11.0"))
                elif changed == "root":
                    (directory / "herwig/events/input.root").write_bytes(b"modified HwSim events")
                else:
                    dat = directory / "signal_gg_h_aa_hwsim_roots-stored_run.dat"
                    before, after = ("weight_scale .003", "weight_scale .006") if changed == "branching" else (
                        "sum_tree_weight .048", "sum_tree_weight .040")
                    dat.write_text(dat.read_text().replace(before, after))
                sidecar_before = (directory / "normalization-stored_run.json").read_bytes()
                with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                        patch.object(ho.signal_normalization, "signal_sidecar") as sidecar:
                    with self.assertRaises(ValueError):
                        ho.run_campaign(self.args(tmp, "--stage", "normalize"))
                    sidecar.assert_not_called()
                self.assertEqual((directory / "normalization-stored_run.json").read_bytes(), sidecar_before)

    def test_normalize_skips_backgrounds_without_touching_them(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp, "--stage", "normalize", "--run-samples", "backgrounds")
            with patch.object(ho, "normalize_existing") as normalize, patch.object(ho, "command") as command:
                self.assertEqual(ho.run_campaign(args), 0)
                normalize.assert_not_called()
                command.assert_not_called()

    def test_normalize_rejects_an_explicit_wrong_run_tag(self):
        with TemporaryDirectory() as tmp:
            self.completed_campaign(tmp)
            with patch.object(ho.signal_normalization, "load_signal_record") as load:
                with self.assertRaisesRegex(ValueError, "--run-tag differs"):
                    ho.run_campaign(self.args(tmp, "--stage", "normalize", "--run-tag", "wrong_run"))
                load.assert_not_called()

    def test_legacy_missing_inventory_is_adopted_with_explicit_snapshot_provenance(self):
        with TemporaryDirectory() as tmp:
            _, directory, manifest = self.completed_campaign(tmp)
            del manifest["root_files"]
            (directory / "campaign.json").write_text(json.dumps(manifest))
            with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                    patch.object(ho.signal_normalization, "signal_sidecar", side_effect=sidecar_fixture), \
                    patch.object(ho, "command") as command:
                ho.run_campaign(self.args(tmp, "--stage", "normalize"))
                command.assert_not_called()
            migrated = json.loads((directory / "campaign.json").read_text())
            self.assertEqual(migrated["root_files"], ho.root_inventory(directory))
            provenance = migrated["root_inventory_provenance"]
            self.assertEqual(provenance["kind"], "legacy_snapshot_at_normalize")
            self.assertFalse(provenance["historical_inventory_available"])
            self.assertEqual(len(provenance["root_sha256"]), 1)
            # Once adopted, content changes are detected even with unchanged
            # size and timestamp metadata.
            root = directory / "herwig/events/input.root"
            old_stat = root.stat()
            root.write_bytes(b"X" + root.read_bytes()[1:])
            os.utime(root, ns=(old_stat.st_atime_ns, old_stat.st_mtime_ns))
            with self.assertRaisesRegex(ValueError, "ROOT contents changed"):
                ho.validate_root_inventory(directory, migrated)

    def test_legacy_missing_inventory_rejects_newer_roots_or_mismatched_input_list(self):
        for changed in ("newer_root", "input_list"):
            with self.subTest(changed=changed), TemporaryDirectory() as tmp:
                _, directory, manifest = self.completed_campaign(tmp)
                del manifest["root_files"]
                if changed == "newer_root":
                    root = directory / "herwig/events/input.root"
                    newer = (directory / "signal_gg_h_aa_hwsim_roots-stored_run.dat").stat().st_mtime_ns + 1000000
                    os.utime(root, ns=(newer, newer))
                else:
                    (directory / "signal_gg_h_aa_hwsim_roots.input").write_text("different-root.root\n")
                with self.assertRaises(ValueError):
                    ho.validate_root_inventory(directory, manifest, allow_legacy_snapshot=True)
                self.assertNotIn("root_files", manifest)

    def test_signal_analyze_passes_br_without_rate_factor_to_detector(self):
        with TemporaryDirectory() as tmp:
            args, directory, manifest = self.completed_campaign(tmp)
            with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                    patch.object(ho.signal_normalization, "signal_sidecar", side_effect=sidecar_fixture), \
                    patch.object(ho, "command") as command:
                ho.analyze(args, ho.SAMPLES[0], manifest)
            argv = command.call_args.args[1]
            self.assertEqual(argv[argv.index("-w") + 1], .003)
            sidecar = json.loads((directory / "normalization-stored_run.json").read_text())
            self.assertEqual(sidecar["cross_section_pb"], 400.)
            self.assertAlmostEqual(1000. * 100. * sidecar["cross_section_pb"] * .003 * (.027 / .048), 67500.)

    def test_resume_requires_matching_record_fingerprint(self):
        with TemporaryDirectory() as tmp:
            args, directory, manifest = self.completed_campaign(tmp)
            sidecar = sidecar_fixture(args.signal_normalization, manifest, manifest["lhe"], .003)
            sidecar["ihixs_record_sha256"] = "b" * 64
            (directory / "normalization-stored_run.json").write_text(json.dumps(sidecar))
            with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                    patch.object(ho.signal_normalization, "validate_sidecar", return_value=sidecar):
                with self.assertRaisesRegex(ValueError, "--stage normalize"):
                    ho.verify_completed(args, ho.SAMPLES[0], "analyze", manifest)

    def test_old_native_sidecar_blocks_resume_before_build_commands(self):
        with TemporaryDirectory() as tmp:
            args, _, _ = self.completed_campaign(tmp)
            args.stage = "analyze"
            args.resume = True
            with patch.object(ho.signal_normalization, "load_signal_record", return_value=RECORD), \
                    patch.object(ho.signal_normalization, "validate_sidecar", side_effect=ValueError("old native sidecar")), \
                    patch.object(ho, "command") as command:
                with self.assertRaisesRegex(ValueError, "--stage normalize"):
                    ho.run_campaign(args)
                command.assert_not_called()

    def test_native_higgs_and_explicit_hard_parameters_must_match(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp, "--nevents", "3")
            args.expected_native_parameters = {"hmass": 125., "tmass": 173.2}
            path = Path(tmp) / "signal.lhe"
            for old, new, message in (("hmass 125", "hmass 126", "hmass"),
                                      ("minnlo 1", "minnlo 1\ntmass 170", "tmass"),
                                      ("alphas_from_pdf 1", "alphas_from_pdf 0", "alphas_from_pdf")):
                with self.subTest(message=message):
                    path.write_text(signal_lhe_text().replace(old, new))
                    with self.assertRaisesRegex(ValueError, message):
                        ho.inspect_lhe(path, ho.SAMPLES[0], args)

    def test_native_checks_ignore_event_counts_and_allow_finite_width_mass(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp, "--nevents", "3")
            args.expected_native_parameters = {"hmass": 125., "numevts": 20000, "iseed": 10}
            path = Path(tmp) / "signal.lhe"
            path.write_text(signal_lhe_text().replace("hmass 125", "hmass 125\nnumevts 100000\niseed 20")
                            .replace("125 125 0 9", "125.001 125.001 0 9"))
            self.assertEqual(ho.inspect_lhe(path, ho.SAMPLES[0], args)["events"], 3)

    def test_native_pure_heft_gate_checks_actual_topmass_and_rate_rescaling(self):
        with TemporaryDirectory() as tmp:
            args = self.args(tmp, "--nevents", "3")
            path = Path(tmp) / "signal.lhe"
            for setting in ("topmass 170", "quarkmasseffects 1", "quarkmasseffects 2",
                            "quarkmasseffects 3", "quarkmasseffects 4", "nnlo 1", "nnloint 1",
                            "nnlopsreweight 1", "lhapdf_in_hoppet 1"):
                with self.subTest(setting=setting):
                    path.write_text(signal_lhe_text().replace("minnlo 1", f"minnlo 1\n{setting}"))
                    with self.assertRaisesRegex(ValueError, setting.split()[0]):
                        ho.inspect_lhe(path, ho.SAMPLES[0], args)
            for settings in ("topmass 173.2\nquarkmasseffects 0\nnnlo 0\nnnloint 0",
                             "topmass -1000000\nquarkmasseffects -1000000\nnnlo -1000000\nnnloint -1000000"):
                with self.subTest(settings=settings):
                    path.write_text(signal_lhe_text().replace("minnlo 1", f"minnlo 1\n{settings}"))
                    self.assertEqual(ho.inspect_lhe(path, ho.SAMPLES[0], args)["events"], 3)


if __name__ == "__main__":
    unittest.main()
