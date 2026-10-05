"""Synthetic fixtures only: these tests never call ihixs or calculate a rate."""

import contextlib
import copy
import hashlib
import io
import json
import math
from pathlib import Path
import statistics
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest import mock

from hgammagamma import ho_signal_normalization as norm
from hgammagamma import run_ihixs_normalization as wrapper


def synthetic_record():
    """Build a deliberately synthetic complete record for schema/closure tests."""
    digest = "a" * 64
    settings = {"synthetic_fixture": True}
    multiplier = norm.GF * math.pi / (math.sqrt(2.) * 288.) * 389379660. / 35.0309
    runs = []
    for index, (r, f) in enumerate(norm.SCALE_POINTS):
        runs.append({"label": "central" if index == 0 else f"scale_{index}",
                     "mur_gev": 62.5 * r, "muf_gev": 62.5 * f, "pdf_member": 0,
                     "pdf_set": norm.PDF_SET, "qcd_order": "N3LO", "cross_section_pb": 100. + index})
    for member in range(1, 101):
        runs.append({"label": f"replica_{member:03d}", "mur_gev": 62.5, "muf_gev": 62.5,
                     "pdf_member": member, "pdf_set": norm.PDF_SET, "qcd_order": "N3LO",
                     "cross_section_pb": 99. + member / 100.})
    for label, order in (("nnlo_native", "NNLO"), ("n3lo_nnlo_pdf", "N3LO")):
        runs.append({"label": label, "mur_gev": 62.5, "muf_gev": 62.5, "pdf_member": 0,
                     "pdf_set": norm.NATIVE_PDF_SET, "qcd_order": order, "cross_section_pb": 98.})
    for run in runs:
        run.update(cross_section_error_pb=.01, input_sha256=digest, output_sha256=digest,
                   pdf_info_sha256=digest, pdf_member_sha256=digest,
                   alpha_s_mur=.125, pdf_alpha_s_mur=.125,
                   alpha_s_at_91_1876=.118,
                   result_key="eftn3lo" if run["qcd_order"] == "N3LO" else "eftnnlo",
                   raw_cross_section_pb=run["cross_section_pb"] / multiplier,
                   raw_cross_section_error_pb=.01 / multiplier)
    pdfs = {}
    for run in runs:
        key = f"{run['pdf_set']}/{run['pdf_member']}"
        if key not in pdfs:
            pdfs[key] = {"info_sha256": digest, "member_sha256": digest,
                         "OrderQCD": 3 if run["pdf_set"] == norm.PDF_SET else 2,
                         "AlphaS_OrderQCD": 3 if run["pdf_set"] == norm.PDF_SET else 2,
                         "AlphaS_MZ": ".118", "alpha_s_91_1876": .118,
                         "photon": True, "NumFlavors": "5", "NumMembers": "101",
                         "ErrorType": "replicas", "alpha_s_by_q_gev": {}}
        pdfs[key]["alpha_s_by_q_gev"][str(run["mur_gev"])] = .125
    record = {"schema_version": 1, "profile": norm.PROFILE, "normalization_kind": "ihixs_n3lo",
              "process": "ggF", "hard_model": "pure_heft", "qcd_order": "N3LO", "top_scheme": "on-shell",
              "pdf_set": norm.PDF_SET, "pdf_member": 0, "electroweak_corrections": False,
              "finite_mass_corrections": False, "resummation": False, "branching_fraction_included": False,
              "alpha_s_source": "LHAPDF::PDF::alphasQ", "sqrt_s_gev": 40000., "higgs_mass_gev": 125.,
              "top_mass_gev": 173.2, "gf_gev_minus2": norm.GF, "mur_gev": 62.5, "muf_gev": 62.5,
              "alpha_s_mz": .118, "cross_section_pb": 100., "cross_section_error_pb": .01,
              "numerical_relative_tolerance": .0005, "runs": runs, "pdfs": pdfs,
              "native_parameters": {"hmass": 125., "hwidth": .004152, "tmass": 173.2, "topmass": 173.2,
                  "ebeam1": 20000., "ebeam2": 20000., "lhans1": 336100., "lhans2": 336100.,
                  "ih1": 1., "ih2": 1., "minlo": 1., "minnlo": 1., "alphas_from_pdf": 1.,
                  "use_NNLOPS_pdfs": 1., "renscfact": 1., "facscfact": 1., "gf_gev_minus2": norm.GF,
                  "hard_model": "pure_heft", "quarkmasseffects": 0., "nnloint": 0., "nnlo": 0.,
                  "nnlopsreweight": 0., "lhapdf_in_hoppet": 0., "alpha_s_evolution": "HOPPET three-loop"},
              "prefactor_correction": {"upstream_prefactor_pb": 35.0309,
                  "gev_minus2_to_pb": 389379660., "multiplier": multiplier},
              "uncertainties": {"scale_up_pb": 6., "scale_down_pb": 0.,
                  "pdf_stddev_pb": statistics.stdev(run["cross_section_pb"] for run in runs[7:107])},
              "provenance": {"ihixs_commit": norm.IHIXS_COMMIT, "settings": settings,
                  "settings_sha256": norm.canonical_sha256(settings), "source_sha256": digest,
                  "adapter_sha256": digest, "executable_sha256": digest, "build_manifest_sha256": digest,
                  "powheg_input_sha256": digest, "artifact_root": "not-copied-from-original-machine",
                  "artifacts": [{"path": "runs/central/ihixs.out", "sha256": digest}]},
              "validations": {"benchmark": {"passed": True, "ihixs_commit": norm.IHIXS_COMMIT,
                  "executable_sha256": digest, "production_executable_sha256": digest,
                  "input_sha256": digest, "output_sha256": digest, "result_key": "eftn3lo",
                  "pdf_set": "PDF4LHC15_nnlo_100", "expected_pb": 45.1816, "relative_tolerance": .005,
                  "raw_cross_section_pb": 45.1816, "raw_cross_section_error_pb": .001},
                  "alpha_s": {"passed": True, "checked_runs": 109}}}
    return norm.seal_record(record)


def synthetic_manifest():
    return {"run_tag": "synthetic-ho", "sample": "signal_gg_h_aa", "matching": "POWHEG",
            "hard_accuracy": "NNLO+PS (HJMiNNLO)", "process": "HJMiNNLO; h -> gamma gamma",
            "pdf": norm.NATIVE_PDF_SET, "lhaid": norm.NATIVE_PDF_ID, "ebeam_gev": 20000.,
            "weight_scale": .00227}


class NormalizationRecordTests(unittest.TestCase):
    def test_complete_record_and_canonical_hash_are_portable(self):
        record = synthetic_record()
        self.assertIs(norm.validate_record(record), record)
        reformatted = json.loads(json.dumps(record, indent=8))
        self.assertEqual(norm.record_sha256(reformatted), record["fingerprint"])
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "migrated.json"
            path.write_text(json.dumps(record))
            self.assertEqual(norm.load_record(path), record)

    def test_mutated_contents_and_missing_record_fail_closed(self):
        record = synthetic_record()
        record["cross_section_pb"] = 200.
        with self.assertRaisesRegex(ValueError, "fingerprint"):
            norm.validate_record(record)
        with TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "calculation first"):
                norm.load_record(Path(tmp) / "not-computed.json")

    def test_no_partial_or_failed_calculation_can_be_published(self):
        for mutation in (lambda r: r["runs"].pop(),
                         lambda r: r["validations"]["benchmark"].update(passed=False),
                         lambda r: r["validations"]["alpha_s"].update(passed=False),
                         lambda r: r["runs"][0].update(cross_section_error_pb=1., raw_cross_section_error_pb=1.)):
            record = synthetic_record()
            mutation(record)
            with self.assertRaises(ValueError):
                norm.validate_record(norm.seal_record(record))

    def test_native_hard_model_mass_pdf_and_scale_must_match(self):
        for key, wrong in (("topmass", 172.5), ("hmass", 126.), ("lhans1", 331100.),
                           ("ebeam1", 6500.), ("minnlo", 0.), ("quarkmasseffects", 1.),
                           ("nnloint", 1.), ("lhapdf_in_hoppet", 1.)):
            with self.subTest(key=key):
                record = synthetic_record()
                record["native_parameters"][key] = wrong
                with self.assertRaises(ValueError):
                    norm.validate_record(norm.seal_record(record))

    def test_scale_and_pdf_uncertainties_are_recomputed_from_runs(self):
        for key in ("scale_up_pb", "scale_down_pb", "pdf_stddev_pb"):
            record = synthetic_record()
            record["uncertainties"][key] += 1.
            with self.assertRaisesRegex(ValueError, key):
                norm.validate_record(norm.seal_record(record))

    def test_replica_comparison_order_and_pdf_coupling_are_checked(self):
        mutations = ((7, "pdf_member", 0), (7, "mur_gev", 125.),
                     (108, "muf_gev", 125.), (108, "pdf_set", norm.PDF_SET),
                     (10, "alpha_s_mur", .118), (10, "pdf_member_sha256", "not-a-hash"))
        for index, key, value in mutations:
            record = synthetic_record()
            record["runs"][index][key] = value
            with self.assertRaises(ValueError):
                norm.validate_record(norm.seal_record(record))

    def test_separate_numerical_and_theory_uncertainty(self):
        record = synthetic_record()
        self.assertLess(record["cross_section_error_pb"], record["uncertainties"]["pdf_stddev_pb"])
        self.assertEqual(norm.validate_record(record)["cross_section_error_pb"], .01)

    def test_portable_sidecar_keeps_native_signed_diagnostics_and_br_once(self):
        manifest = synthetic_manifest()
        native = {"cross_section_pb": 80., "cross_section_error_pb": 2.,
                  "sum_weight": 16., "sum_abs_weight": 20., "negative_events": 1}
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "record.json"
            path.write_text(json.dumps(synthetic_record()))
            sidecar = norm.signal_sidecar(path, manifest, native, .00227)
        self.assertEqual(norm.validate_sidecar(sidecar, manifest, .00227)["cross_section_pb"], 100.)
        self.assertEqual(sidecar["native_lhe"]["sum_weight"], 16.)
        self.assertEqual(sidecar["native_lhe"]["sum_abs_weight"], 20.)
        self.assertEqual(sidecar["native_cross_section_pb"], 80.)
        self.assertEqual(sidecar["weight_scale"], .00227)

    def test_sidecar_identity_rate_and_branching_fraction_cannot_drift(self):
        manifest = synthetic_manifest()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "record.json"
            path.write_text(json.dumps(synthetic_record()))
            base = norm.signal_sidecar(path, manifest, {"cross_section_pb": 80., "cross_section_error_pb": 1.}, .00227)
        for key, value in (("cross_section_pb", 80.), ("run_tag", "other"),
                           ("weight_scale", .00227 * .00227), ("ihixs_record_sha256", "0" * 64)):
            sidecar = copy.deepcopy(base)
            sidecar[key] = value
            with self.assertRaises(ValueError):
                norm.validate_sidecar(sidecar, manifest, .00227)

    def test_only_actual_ho_signal_is_eligible(self):
        manifest = synthetic_manifest()
        self.assertTrue(norm.is_ho_signal(manifest))
        for key, value in (("matching", "LO"), ("hard_accuracy", "LO"), ("sample", "bkg_prompt_aa")):
            self.assertFalse(norm.is_ho_signal({**manifest, key: value}))
        for key, value in (("ebeam_gev", 6500.), ("pdf", "wrong"), ("lhaid", 1)):
            bad = {**manifest, key: value}
            self.assertTrue(norm.is_ho_signal(bad))
            with self.assertRaises(ValueError):
                norm.validate_manifest(synthetic_record(), bad)


class IhixsWrapperTests(unittest.TestCase):
    def test_raw_eft_parser_rejects_rescaled_only_output(self):
        text = "Higgs XS = 108.0 [0.01]\nR_LO*eftn3lo = 108.0 [0.01]\n"
        with self.assertRaisesRegex(ValueError, "raw eftn3lo"):
            wrapper.parse_output(text)
        text += "eftn3lo = 1.00000000e+02 [1.00000000e-02]\nas_at_mz = 1.18000000e-01\nas_at_mur = 1.25000000e-01\n"
        result = wrapper.parse_output(text)
        self.assertEqual(result["raw_cross_section_pb"], 100.)
        self.assertEqual(result["raw_cross_section_error_pb"], .01)
        self.assertEqual(result["alpha_s_mur"], .125)
        with self.assertRaises(ValueError):
            wrapper.parse_output(text + "eftn3lo = 100 [0.01]\n")

    def test_campaign_has_exact_seven_point_envelope_all_replicas_and_comparisons(self):
        settings = wrapper.load_settings(wrapper.DEFAULT_SETTINGS)
        points = wrapper.calculation_points(settings)
        self.assertEqual(len(points), 109)
        self.assertNotIn((31.25, 125.), {(p["mur_gev"], p["muf_gev"]) for p in points[:7]})
        self.assertNotIn((125., 31.25), {(p["mur_gev"], p["muf_gev"]) for p in points[:7]})
        self.assertEqual({p["pdf_member"] for p in points[7:107]}, set(range(1, 101)))
        self.assertEqual(points[-2]["qcd_order"], "NNLO")
        self.assertEqual(points[-1]["qcd_order"], "N3LO")
        self.assertEqual(points[-1]["pdf_set"], norm.NATIVE_PDF_SET)

    def test_cards_enforce_pure_eft_and_keep_uncertainties_external(self):
        settings = wrapper.load_settings(wrapper.DEFAULT_SETTINGS)
        card = wrapper.production_card(settings, wrapper.calculation_points(settings)[0])
        for option in ("top_scheme = on-shell", "mt_on_shell = 173.2", "y_bot = 0.0", "y_charm = 0.0",
                       "with_exact_qcd_corrections = false", "with_ew_corrections = false",
                       "with_resummation = false", "with_scet = false", "with_scale_variation = false",
                       "with_pdf_error = false", "qcd_order_evol = 3"):
            self.assertIn(option, card)
        self.assertEqual(settings["gev_minus2_to_pb"], 389379660.)
        self.assertNotEqual(wrapper.prefactor_multiplier(settings), norm.GF / 1.16637e-5)

    def test_native_mass_reweight_and_scale_conflicts_are_rejected(self):
        original = (wrapper.SCRIPT_DIR / "HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input").read_text()
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "powheg.input"
            path.write_text(original + "\ntopmass -1\nquarkmasseffects -1\nnnloint -1\n")
            native = wrapper.native_parameters(path)
            self.assertEqual(native["topmass"], 173.2)
            self.assertEqual(native["quarkmasseffects"], 0.)
            self.assertEqual(native["use_NNLOPS_pdfs"], 1.)
            for addition in ("topmass 172.5", "quarkmasseffects 1", "nnloint 1", "renscfact 2",
                             "lhapdf_in_hoppet 1"):
                path.write_text(original + "\n" + addition + "\n")
                with self.assertRaises(ValueError):
                    wrapper.native_parameters(path)

    def test_dry_run_never_calls_external_commands_or_writes_outputs(self):
        with TemporaryDirectory() as tmp:
            work = Path(tmp) / "never-created"
            with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")):
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(wrapper.main(["--stage", "calculate", "--dry-run", "--work-dir", str(work)]), 0)
            self.assertFalse(work.exists())

    def test_prepare_writes_no_fabricated_result_record(self):
        with TemporaryDirectory() as tmp:
            work = Path(tmp) / "work"
            record = Path(tmp) / "result.json"
            with mock.patch.object(wrapper, "command", side_effect=AssertionError("external command called")):
                with contextlib.redirect_stdout(io.StringIO()):
                    wrapper.main(["--work-dir", str(work), "--record", str(record)])
            self.assertTrue((work / "plan.json").is_file())
            self.assertEqual(len(list((work / "cards").glob("*.card"))), 109)
            self.assertFalse(record.exists())

    def test_modified_settings_and_loose_precision_fail_closed(self):
        settings = json.loads(wrapper.DEFAULT_SETTINGS.read_text())
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "settings.json"
            for key, value in (("top_mass_gev", 172.5), ("numerical_relative_tolerance", .001),
                               ("electroweak_corrections", True)):
                path.write_text(json.dumps({**settings, key: value}))
                with self.assertRaises(ValueError):
                    wrapper.load_settings(path)


class PrecisionRefinementTests(unittest.TestCase):
    """Mock integrations: refinement changes numerical inputs, never physics."""

    def fixture(self, tmp, *, refine_failed=False):
        work = Path(tmp) / "work"
        args = SimpleNamespace(work_dir=work, resume=True, refine_failed=refine_failed)
        settings = json.loads(wrapper.DEFAULT_SETTINGS.read_text())
        settings["integration"].update(epsrel=5e-5, epsabs=.025,
                                       mineval=200000, maxeval=200000000,
                                       nstart=40000, nincrease=4000)
        point = {"label": "scale_2", "mur_gev": 31.25, "muf_gev": 62.5,
                 "pdf_member": 0, "pdf_set": norm.PDF_SET, "qcd_order": "N3LO"}
        identity = {"settings_sha256": norm.canonical_sha256(settings), "pdf_info_sha256": "b" * 64,
                    "pdf_member_sha256": "c" * 64, "probe_sha256": "d" * 64,
                    "variant": "lhapdf"}
        pdf = {"alpha_s_mur": .125, "alpha_s_91_1876": .118}
        executable = work / "build-lhapdf/ihixs"
        return args, settings, point, executable, identity, pdf

    @staticmethod
    def result(error=.01, **updates):
        return {"result_key": "eftn3lo", "raw_cross_section_pb": 100.,
                "raw_cross_section_error_pb": error, "alpha_s_mur": .125,
                "alpha_s_at_91_1876": .118, **updates}

    def fake_runs(self):
        """A run_card boundary double retaining prior inputs/output verbatim."""
        calls = []
        created = []

        def run(args, label, card, executable, identity):
            calls.append((label, card, executable, copy.deepcopy(identity)))
            directory = args.work_dir / "runs" / label
            completed = {"identity": {**identity,
                         "input_sha256": hashlib.sha256(card.encode()).hexdigest()},
                         "output_sha256": hashlib.sha256(label.encode()).hexdigest()}
            if directory.exists():
                self.assertEqual((directory / "input.card").read_text(), card)
                self.assertEqual(json.loads((directory / "complete.json").read_text()), completed)
            else:
                directory.mkdir(parents=True)
                (directory / "input.card").write_text(card)
                (directory / "ihixs.out").write_text("synthetic output: " + label)
                (directory / "complete.json").write_text(json.dumps(completed))
                created.append(label)
            return directory, completed

        return run, calls, created

    @staticmethod
    def card_options(card):
        return dict(line.split(" = ", 1) for line in card.splitlines()
                    if " = " in line and not line.lstrip().startswith("#"))

    def test_valid_base_point_needs_no_flag_or_refinement(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp)
            # Older callers without the new optional attribute retain strict behavior.
            del args.refine_failed
            run, calls, created = self.fake_runs()
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", return_value=self.result()) as parsed:
                directory, completed, result = wrapper.run_precision_checked(
                    args, settings, point, executable, identity, pdf)
            self.assertEqual(directory, args.work_dir / "runs/scale_2")
            self.assertEqual(created, ["scale_2"])
            self.assertEqual(calls[0][1], wrapper.production_card(settings, point))
            self.assertEqual(calls[0][3], identity)
            self.assertEqual(result["precision_refinements"], 0)
            self.assertEqual(result["run_directory"], "runs/scale_2")
            multiplier = wrapper.prefactor_multiplier(settings)
            self.assertAlmostEqual(result["cross_section_pb"], 100. * multiplier)
            self.assertAlmostEqual(result["cross_section_error_pb"], .01 * multiplier)
            self.assertEqual(completed["identity"]["settings_sha256"], identity["settings_sha256"])
            parsed.assert_called_once_with("synthetic output: scale_2", "N3LO")

    def test_precision_failure_without_flag_stops_after_base(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp)
            run, calls, created = self.fake_runs()
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", return_value=self.result(.1)):
                with self.assertRaisesRegex(ValueError, "numerical uncertainty"):
                    wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
            self.assertEqual([call[0] for call in calls], ["scale_2"])
            self.assertEqual(created, ["scale_2"])

    def test_only_failed_point_is_refined_and_all_original_artifacts_survive(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
            original_settings = copy.deepcopy(settings)
            central = {**point, "label": "central", "mur_gev": 62.5}
            run, calls, created = self.fake_runs()
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", side_effect=[
                        self.result(), self.result(.1), self.result(.01)]):
                central_dir, _, central_result = wrapper.run_precision_checked(
                    args, settings, central, executable, identity, pdf)
                central_contents = {p.name: p.read_bytes() for p in central_dir.iterdir()}
                directory, _, result = wrapper.run_precision_checked(
                    args, settings, point, executable, identity, pdf)
            self.assertEqual(created, ["central", "scale_2", "scale_2__precision_1"])
            self.assertEqual([call[0] for call in calls], created)
            self.assertEqual(central_result["precision_refinements"], 0)
            self.assertEqual(directory, args.work_dir / "runs/scale_2__precision_1")
            self.assertEqual(result["precision_refinements"], 1)
            self.assertEqual(result["run_directory"], "runs/scale_2__precision_1")
            self.assertEqual({p.name: p.read_bytes() for p in central_dir.iterdir()}, central_contents)
            self.assertEqual((args.work_dir / "runs/scale_2/input.card").read_text(),
                             wrapper.production_card(original_settings, point))
            self.assertEqual((args.work_dir / "runs/scale_2/ihixs.out").read_text(),
                             "synthetic output: scale_2")
            self.assertEqual(settings, original_settings)
            base_options = self.card_options(calls[1][1])
            retry_options = self.card_options(calls[2][1])
            for key in base_options:
                if key not in ("epsrel", "epsabs"):
                    self.assertEqual(retry_options[key], base_options[key], key)
            self.assertAlmostEqual(float(retry_options["epsrel"]), 1e-6)
            self.assertAlmostEqual(float(retry_options["epsabs"]), .0005)
            self.assertEqual(result["integration"]["epsrel"], 1e-6)
            self.assertAlmostEqual(result["integration"]["epsabs"], .0005)
            self.assertEqual(calls[2][2], executable)
            expected_settings = copy.deepcopy(original_settings)
            expected_settings["integration"]["epsrel"] = 1e-6
            expected_settings["integration"]["epsabs"] *= 1e-6 / original_settings["integration"]["epsrel"]
            self.assertEqual(calls[2][3], {**identity,
                             "settings_sha256": norm.canonical_sha256(expected_settings)})

    def test_verified_refinement_is_revisited_via_same_run_card_identity(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
            run, calls, created = self.fake_runs()
            outputs = [self.result(.1), self.result(.01)] * 2
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", side_effect=outputs):
                first = wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
                saved = {str(p.relative_to(args.work_dir)): p.read_bytes()
                         for p in (args.work_dir / "runs").rglob("*") if p.is_file()}
                second = wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
            self.assertEqual(first, second)
            self.assertEqual(created, ["scale_2", "scale_2__precision_1"])
            self.assertEqual([call[0] for call in calls], created * 2)
            self.assertEqual(calls[:2], calls[2:])
            self.assertEqual({str(p.relative_to(args.work_dir)): p.read_bytes()
                              for p in (args.work_dir / "runs").rglob("*") if p.is_file()}, saved)

    def test_three_refinements_are_tenfold_tighter_and_exhaustion_fails_closed(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
            run, calls, created = self.fake_runs()
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", return_value=self.result(.1)):
                with self.assertRaisesRegex(ValueError, "numerical uncertainty"):
                    wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
            labels = ["scale_2"] + [f"scale_2__precision_{index}" for index in range(1, 4)]
            self.assertEqual(created, labels)
            self.assertEqual([call[0] for call in calls], labels)
            for index, call in enumerate(calls[1:]):
                options = self.card_options(call[1])
                self.assertAlmostEqual(float(options["epsrel"]), 1e-6 / (10 ** index), places=12)
                self.assertAlmostEqual(float(options["epsabs"]), .0005 / (10 ** index), places=12)

    def test_already_tight_base_gets_tighter_and_second_success_is_selected(self):
        with TemporaryDirectory() as tmp:
            args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
            settings["integration"]["epsrel"] = 5e-7
            identity["settings_sha256"] = norm.canonical_sha256(settings)
            run, calls, created = self.fake_runs()
            with mock.patch.object(wrapper, "run_card", side_effect=run), \
                    mock.patch.object(wrapper, "parse_output", side_effect=[
                        self.result(.1), self.result(.08), self.result(.01)]):
                directory, _, result = wrapper.run_precision_checked(
                    args, settings, point, executable, identity, pdf)
            self.assertEqual(created, ["scale_2", "scale_2__precision_1", "scale_2__precision_2"])
            self.assertEqual(directory, args.work_dir / "runs/scale_2__precision_2")
            self.assertEqual(result["run_directory"], "runs/scale_2__precision_2")
            self.assertEqual(result["precision_refinements"], 2)
            self.assertAlmostEqual(float(self.card_options(calls[1][1])["epsrel"]), 5e-8, places=14)
            self.assertAlmostEqual(result["integration"]["epsrel"], 5e-9, places=14)
            self.assertAlmostEqual(result["integration"]["epsabs"], .00025)

    def test_alpha_s_mismatch_is_never_retried_even_with_large_integration_error(self):
        for wrong in ({"alpha_s_mur": .126}, {"alpha_s_at_91_1876": .119}):
            with self.subTest(wrong=wrong), TemporaryDirectory() as tmp:
                args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
                run, calls, _ = self.fake_runs()
                with mock.patch.object(wrapper, "run_card", side_effect=run), \
                        mock.patch.object(wrapper, "parse_output", return_value=self.result(.1, **wrong)):
                    with self.assertRaisesRegex(ValueError, "alpha_s"):
                        wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
                self.assertEqual([call[0] for call in calls], ["scale_2"])

    def test_execution_and_parser_errors_are_not_precision_refinements(self):
        for error, target in ((RuntimeError("synthetic process failure"), "run_card"),
                              (ValueError("ambiguous raw EFT output"), "parse_output")):
            with self.subTest(target=target), TemporaryDirectory() as tmp:
                args, settings, point, executable, identity, pdf = self.fixture(tmp, refine_failed=True)
                run, calls, _ = self.fake_runs()
                with mock.patch.object(wrapper, "run_card", side_effect=run) as runner, \
                        mock.patch.object(wrapper, "parse_output", return_value=self.result()) as parsed:
                    (runner if target == "run_card" else parsed).side_effect = error
                    with self.assertRaisesRegex(type(error), str(error)):
                        wrapper.run_precision_checked(args, settings, point, executable, identity, pdf)
                self.assertEqual(runner.call_count, 1)
                if target == "parse_output":
                    self.assertEqual([call[0] for call in calls], ["scale_2"])

    def test_refinement_flag_requires_calculate_and_resume(self):
        self.assertFalse(wrapper.parse_args(["--stage", "calculate", "--resume"]).refine_failed)
        args = wrapper.parse_args(["--stage", "calculate", "--resume", "--refine-failed"])
        self.assertTrue(args.refine_failed)
        for argv in (["--stage", "calculate", "--refine-failed"],
                     ["--refine-failed", "--resume"],
                     ["--stage", "build", "--resume", "--refine-failed"],
                     ["--stage", "benchmark", "--resume", "--refine-failed"]):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    wrapper.parse_args(argv)
                self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
