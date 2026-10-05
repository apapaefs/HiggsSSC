"""Validated, portable inclusive normalization for the HJMiNNLO HO signal.

Only the HO diphoton ggF sample consumes these records.  The production rate
excludes the Higgs branching fraction, detector response and generator cuts.
Numerical integration errors and theory uncertainties remain separate.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
import re
import statistics


DEFAULT_RECORD = Path(__file__).resolve().parent / "HOAnalysis/normalization/ggf-ssc40-n3lo.json"
IHIXS_COMMIT = "316d5cdf3e88f2ba69cc43668d66b4f69ff833ce"
PDF_SET = "NNPDF40_an3lo_as_01180_qed_mhou"
NATIVE_PDF_SET = "NNPDF40_nnlo_as_01180_qed"
NATIVE_PDF_ID = 336100
GF = 1.1663899987193257e-5
NUMERICAL_TOLERANCE = 0.0005
PROFILE = "ssc40-ho-pure-heft-v1"
SCALE_POINTS = ((1., 1.), (.5, .5), (.5, 1.), (1., .5), (1., 2.), (2., 1.), (2., 2.))


def canonical_sha256(payload):
    """Hash JSON content, independent of formatting and file location."""
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def record_sha256(record):
    return canonical_sha256({key: value for key, value in record.items() if key != "fingerprint"})


def seal_record(record):
    result = copy.deepcopy(record)
    result["fingerprint"] = record_sha256(result)
    return result


def finite_number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if positive and value <= 0 or nonnegative and value < 0:
        raise ValueError(f"{name} has an invalid sign")
    return float(value)


def _equal_number(value, expected, name):
    if not math.isclose(finite_number(value, name), expected, rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f"{name} differs from the {PROFILE} normalization profile")


def _hash(value, name):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{name} is not a SHA256 digest")


def _validate_benchmark(benchmark, production_executable, name="benchmark"):
    """Bind a published-example check to its production executable."""
    if not isinstance(benchmark, dict) or benchmark.get("passed") is not True:
        raise ValueError(f"{name} validation must pass")
    if benchmark["ihixs_commit"] != IHIXS_COMMIT:
        raise ValueError(f"{name} used a different ihixs source")
    for key in ("executable_sha256", "production_executable_sha256", "input_sha256", "output_sha256"):
        _hash(benchmark[key], f"{name} {key}")
    if benchmark["production_executable_sha256"] != production_executable:
        raise ValueError(f"{name} validation refers to a different production build")
    if benchmark["result_key"] != "eftn3lo" or benchmark["pdf_set"] != "PDF4LHC15_nnlo_100":
        raise ValueError(f"{name} is not the published raw-EFT example")
    _equal_number(benchmark["expected_pb"], 45.1816, f"{name} reference")
    tolerance = finite_number(benchmark["relative_tolerance"], f"{name} tolerance", positive=True)
    sigma = finite_number(benchmark["raw_cross_section_pb"], f"{name} rate", positive=True)
    error = finite_number(benchmark["raw_cross_section_error_pb"], f"{name} error", nonnegative=True)
    if tolerance > .005 or abs(sigma / 45.1816 - 1.) > tolerance or error / sigma > .0005:
        raise ValueError(f"published ihixs {name} tolerance failed")


def _validate_parser_repair_central(check, record, central, repair):
    """Check the repaired binary against the preserved original central rate."""
    if central["build_variant"] != "lhapdf":
        raise ValueError("parser repair central baseline must use the original production build")
    if not isinstance(check, dict) or check.get("passed") is not True:
        raise ValueError("parser repair central validation must pass")
    expected = {"result_key": "eftn3lo", "qcd_order": "N3LO", "pdf_set": PDF_SET,
                "pdf_member": 0, "build_variant": "lhapdf-parser-repair"}
    for key, value in expected.items():
        if check[key] != value:
            raise ValueError(f"parser repair central has incompatible {key}")
    for key in ("mur_gev", "muf_gev"):
        _equal_number(check[key], 62.5, f"parser repair central {key}")
    for key in ("input_sha256", "output_sha256", "pdf_info_sha256", "pdf_member_sha256",
                "executable_sha256"):
        _hash(check[key], f"parser repair central {key}")
    if check["executable_sha256"] != repair["executable_sha256"]:
        raise ValueError("parser repair central used a different repaired executable")
    for key in ("pdf_info_sha256", "pdf_member_sha256"):
        if check[key] != central[key]:
            raise ValueError("parser repair central PDF hashes differ from the original central")
    for key in ("cross_section_pb", "cross_section_error_pb"):
        _equal_number(check[f"original_{key}"], central[key], f"parser repair central original {key}")
    tolerance = record["numerical_relative_tolerance"]
    _equal_number(check["relative_tolerance"], tolerance, "parser repair central tolerance")
    sigma = finite_number(check["cross_section_pb"], "parser repair central rate", positive=True)
    error = finite_number(check["cross_section_error_pb"], "parser repair central error", nonnegative=True)
    difference = abs(sigma / central["cross_section_pb"] - 1.)
    _equal_number(check["relative_difference"], difference, "parser repair central relative difference")
    if difference > tolerance or error / sigma > tolerance:
        raise ValueError("parser repair central numerical agreement failed")
    pdf = record["pdfs"][f"{PDF_SET}/0"]
    for observed, expected_as in ((check["alpha_s_mur"], central["pdf_alpha_s_mur"]),
                                 (check["pdf_alpha_s_mur"], central["pdf_alpha_s_mur"]),
                                 (check["alpha_s_at_91_1876"], pdf["alpha_s_91_1876"])):
        if not math.isclose(finite_number(observed, "parser repair central alpha_s", positive=True),
                            expected_as, rel_tol=1e-7, abs_tol=1e-10):
            raise ValueError("parser repair central alpha_s differs from its PDF")
    multiplier = GF * math.pi / (math.sqrt(2.) * 288.) * 389379660. / 35.0309
    for raw_key, key, nonnegative in (("raw_cross_section_pb", "cross_section_pb", False),
                                     ("raw_cross_section_error_pb", "cross_section_error_pb", True)):
        raw = finite_number(check[raw_key], f"parser repair central {raw_key}",
                            positive=not nonnegative, nonnegative=nonnegative)
        _equal_number(raw * multiplier, check[key], f"parser repair central GF-corrected {key}")
    refinements = check["precision_refinements"]
    if isinstance(refinements, bool) or not isinstance(refinements, int) or not 0 <= refinements <= 3:
        raise ValueError("parser repair central precision refinements are invalid")
    suffix = f"__precision_{refinements}" if refinements else ""
    if check["run_directory"] != f"parser-repair/runs/central_check{suffix}":
        raise ValueError("parser repair central run directory disagrees with its refinement")
    expected_integration = copy.deepcopy(central["integration"] if "integration" in central
                                          else record["provenance"]["settings"]["integration"])
    if not isinstance(expected_integration, dict):
        raise ValueError("parser repair central original integration settings are missing")
    if refinements:
        base_epsrel = finite_number(expected_integration["epsrel"], "original integration epsrel", positive=True)
        finite_number(expected_integration["epsabs"], "original integration epsabs", nonnegative=True)
        epsrel = min(base_epsrel / 10., 1e-6) / 10. ** (refinements - 1)
        expected_integration["epsrel"] = epsrel
        expected_integration["epsabs"] *= epsrel / base_epsrel
    integration = check["integration"]
    if not isinstance(integration, dict) or set(integration) != set(expected_integration):
        raise ValueError("parser repair central integration settings are incomplete")
    for key, value in expected_integration.items():
        _equal_number(integration[key], value, f"parser repair central integration {key}")


def validate_record(record):
    """Check calculation completeness and provenance without requiring old paths."""
    if not isinstance(record, dict):
        raise ValueError("ihixs record must be a JSON object")
    try:
        _hash(record["fingerprint"], "record fingerprint")
        if record["fingerprint"] != record_sha256(record):
            raise ValueError("ihixs record fingerprint does not match its contents")
        expected = {"schema_version": 1, "profile": PROFILE, "normalization_kind": "ihixs_n3lo",
                    "process": "ggF", "hard_model": "pure_heft", "qcd_order": "N3LO",
                    "top_scheme": "on-shell", "pdf_set": PDF_SET, "pdf_member": 0,
                    "electroweak_corrections": False, "resummation": False,
                    "finite_mass_corrections": False, "branching_fraction_included": False,
                    "alpha_s_source": "LHAPDF::PDF::alphasQ"}
        for key, value in expected.items():
            if record[key] != value:
                raise ValueError(f"ihixs record has incompatible {key}")
        for key, value in {"sqrt_s_gev": 40000., "higgs_mass_gev": 125.,
                           "top_mass_gev": 173.2, "gf_gev_minus2": GF,
                           "mur_gev": 62.5, "muf_gev": 62.5, "alpha_s_mz": .118}.items():
            _equal_number(record[key], value, key)
        sigma = finite_number(record["cross_section_pb"], "cross section", positive=True)
        error = finite_number(record["cross_section_error_pb"], "numerical error", nonnegative=True)
        tolerance = finite_number(record["numerical_relative_tolerance"], "numerical tolerance", positive=True)
        if tolerance > NUMERICAL_TOLERANCE or error / sigma > tolerance:
            raise ValueError("ihixs integration error exceeds the required 0.05% tolerance")
        provenance = record["provenance"]
        if provenance["ihixs_commit"] != IHIXS_COMMIT:
            raise ValueError("ihixs source commit differs from the pinned version")
        for key in ("settings_sha256", "source_sha256", "adapter_sha256", "executable_sha256",
                    "build_manifest_sha256", "powheg_input_sha256"):
            _hash(provenance[key], key)
        artifacts = provenance["artifacts"]
        if not isinstance(artifacts, list) or not artifacts:
            raise ValueError("ihixs record is missing its input/output artifact hashes")
        for artifact in artifacts:
            _hash(artifact["sha256"], "artifact SHA256")
            if not isinstance(artifact["path"], str) or not artifact["path"]:
                raise ValueError("ihixs artifact path is missing")
        validations = record["validations"]
        if validations["benchmark"]["passed"] is not True or validations["alpha_s"]["passed"] is not True:
            raise ValueError("ihixs benchmark and PDF alpha_s validation must pass")
        _validate_benchmark(validations["benchmark"], provenance["executable_sha256"])
        executables = {"lhapdf": provenance["executable_sha256"]}
        has_parser_repair = "parser_repair" in provenance
        if has_parser_repair:
            repair = provenance["parser_repair"]
            if not isinstance(repair, dict) or repair["ihixs_commit"] != IHIXS_COMMIT:
                raise ValueError("parser repair used a different ihixs source")
            for key in ("source_sha256", "adapter_sha256", "parser_patch_sha256",
                        "executable_sha256", "build_manifest_sha256"):
                _hash(repair[key], f"parser repair {key}")
            for key in ("source_sha256", "adapter_sha256"):
                if repair[key] != provenance[key]:
                    raise ValueError(f"parser repair {key} differs from the original build")
            _validate_benchmark(validations["parser_repair_benchmark"], repair["executable_sha256"],
                                "parser repair benchmark")
            benchmark_patch = validations["parser_repair_benchmark"]["parser_patch_sha256"]
            _hash(benchmark_patch, "parser repair benchmark parser_patch_sha256")
            if benchmark_patch != repair["parser_patch_sha256"]:
                raise ValueError("parser repair benchmark used a different parser patch")
            executables["lhapdf-parser-repair"] = repair["executable_sha256"]
        elif "parser_repair_benchmark" in validations or "parser_repair_central" in validations:
            raise ValueError("parser repair validation lacks its build provenance")
        if validations["alpha_s"]["checked_runs"] != 109:
            raise ValueError("alpha_s validation did not cover all calculation points")
        # The central, six noncentral scales, all replicas and two order/PDF
        # comparisons must be present. A central-only run is never publishable.
        runs = record["runs"]
        if len(runs) != 109 or len({run["label"] for run in runs}) != 109:
            raise ValueError("ihixs record requires 7 scale points, 100 replicas and 2 comparisons")
        by_label = {run["label"]: run for run in runs}
        for label in ["central", *[f"scale_{index}" for index in range(1, 7)],
                      *[f"replica_{index:03d}" for index in range(1, 101)],
                      "nnlo_native", "n3lo_nnlo_pdf"]:
            if label not in by_label:
                raise ValueError(f"ihixs record is missing {label}")
        for run in runs:
            # Historical records omit both fields. Mixed-build records must
            # identify every run, including those retained from the base build.
            if has_parser_repair or "executable_sha256" in run or "build_variant" in run:
                _hash(run["executable_sha256"], "run executable hash")
                variant = run["build_variant"]
                if not isinstance(variant, str) or variant not in executables:
                    raise ValueError(f"{run['label']} has an unknown build variant")
                if run["executable_sha256"] != executables[variant]:
                    raise ValueError(f"{run['label']} executable differs from its build variant")
            val = finite_number(run["cross_section_pb"], "run cross section", positive=True)
            err = finite_number(run["cross_section_error_pb"], "run numerical error", nonnegative=True)
            if err / val > tolerance:
                raise ValueError(f"{run['label']} exceeds the numerical tolerance")
            _hash(run["input_sha256"], "run input hash")
            _hash(run["output_sha256"], "run output hash")
            _hash(run["pdf_info_sha256"], "run PDF info hash")
            _hash(run["pdf_member_sha256"], "run PDF member hash")
            if run["result_key"] != {"N3LO": "eftn3lo", "NNLO": "eftnnlo"}.get(run["qcd_order"]):
                raise ValueError("run result is not the raw EFT result at its requested order")
            pdf = record["pdfs"][f"{run['pdf_set']}/{run['pdf_member']}"]
            if pdf["info_sha256"] != run["pdf_info_sha256"] or pdf["member_sha256"] != run["pdf_member_sha256"]:
                raise ValueError("run PDF hashes disagree with its PDF provenance")
            order = 3 if run["pdf_set"] == PDF_SET else 2
            if int(pdf["OrderQCD"]) != order or int(pdf["AlphaS_OrderQCD"]) != order:
                raise ValueError("PDF/alpha_s order differs from the normalization reference")
            _equal_number(float(pdf["AlphaS_MZ"]), .118, "PDF alpha_s(MZ)")
            if pdf["photon"] is not True or int(pdf["NumFlavors"]) != 5:
                raise ValueError("PDF provenance lacks the required five-flavour QED convention")
            if run["pdf_set"] == PDF_SET and (int(pdf["NumMembers"]) != 101 or pdf["ErrorType"] != "replicas"):
                raise ValueError("aN3LO PDF replica ensemble is incompatible")
            observed = finite_number(run["alpha_s_mur"], "run alpha_s", positive=True)
            expected_as = finite_number(run["pdf_alpha_s_mur"], "PDF alpha_s", positive=True)
            if not math.isclose(observed, expected_as, rel_tol=1e-7, abs_tol=1e-10):
                raise ValueError(f"{run['label']} alpha_s does not match its PDF")
            _equal_number(pdf["alpha_s_by_q_gev"][str(run["mur_gev"])], expected_as, "PDF alpha_s grid sample")
            if not math.isclose(finite_number(run["alpha_s_at_91_1876"], "input alpha_s(MZ)", positive=True),
                                finite_number(pdf["alpha_s_91_1876"], "PDF alpha_s at 91.1876", positive=True),
                                rel_tol=1e-7, abs_tol=1e-10):
                raise ValueError("run alpha_s(MZ) differs from the PDF reference")
        _equal_number(by_label["central"]["cross_section_pb"], sigma, "central result")
        _equal_number(by_label["central"]["cross_section_error_pb"], error, "central error")
        if has_parser_repair:
            _validate_parser_repair_central(validations["parser_repair_central"], record,
                                            by_label["central"], repair)
        for run, point in zip([by_label["central"], *[by_label[f"scale_{i}"] for i in range(1, 7)]], SCALE_POINTS):
            _equal_number(run["mur_gev"], 62.5 * point[0], "scale muR")
            _equal_number(run["muf_gev"], 62.5 * point[1], "scale muF")
            if run["pdf_set"] != PDF_SET or run["pdf_member"] != 0 or run["qcd_order"] != "N3LO":
                raise ValueError("scale run PDF/order mismatch")
        for index in range(1, 101):
            run = by_label[f"replica_{index:03d}"]
            if run["pdf_member"] != index or run["pdf_set"] != PDF_SET or run["qcd_order"] != "N3LO":
                raise ValueError("replica run PDF/member/order mismatch")
            _equal_number(run["mur_gev"], 62.5, "replica muR")
            _equal_number(run["muf_gev"], 62.5, "replica muF")
        expected_pdf_keys = {f"{PDF_SET}/{index}" for index in range(101)} | {f"{NATIVE_PDF_SET}/0"}
        if set(record["pdfs"]) != expected_pdf_keys:
            raise ValueError("PDF provenance does not contain the full required ensemble")
        for label, order in (("nnlo_native", "NNLO"), ("n3lo_nnlo_pdf", "N3LO")):
            run = by_label[label]
            if run["pdf_set"] != NATIVE_PDF_SET or run["pdf_member"] != 0 or run["qcd_order"] != order:
                raise ValueError("NNLO comparison PDF/order mismatch")
            _equal_number(run["mur_gev"], 62.5, "comparison muR")
            _equal_number(run["muf_gev"], 62.5, "comparison muF")
        for key in ("scale_up_pb", "scale_down_pb", "pdf_stddev_pb"):
            finite_number(record["uncertainties"][key], key, nonnegative=True)
        rates = [by_label["central"]["cross_section_pb"],
                 *[by_label[f"scale_{i}"]["cross_section_pb"] for i in range(1, 7)]]
        replicas = [by_label[f"replica_{i:03d}"]["cross_section_pb"] for i in range(1, 101)]
        for key, expected_uncertainty in {"scale_up_pb": max(rates) - sigma,
                                          "scale_down_pb": sigma - min(rates),
                                          "pdf_stddev_pb": statistics.stdev(replicas)}.items():
            _equal_number(record["uncertainties"][key], expected_uncertainty, key)
        native = record["native_parameters"]
        for key, value in {"hmass": 125., "hwidth": .004152, "tmass": 173.2, "topmass": 173.2,
                           "ebeam1": 20000., "ebeam2": 20000., "lhans1": 336100., "lhans2": 336100.,
                           "ih1": 1., "ih2": 1., "minlo": 1., "minnlo": 1.,
                           "alphas_from_pdf": 1., "use_NNLOPS_pdfs": 1.,
                           "renscfact": 1., "facscfact": 1., "gf_gev_minus2": GF}.items():
            _equal_number(native[key], value, f"native {key}")
        if native["hard_model"] != "pure_heft":
            raise ValueError("native signal hard model differs from pure HEFT")
        if any(finite_number(native[key], f"native {key}") != 0. for key in
               ("quarkmasseffects", "nnloint", "nnlo", "nnlopsreweight", "lhapdf_in_hoppet")):
            raise ValueError("native signal mass effects/reweighting/PDF evolution differs from its reference")
        if not isinstance(native["alpha_s_evolution"], str) or "three-loop" not in native["alpha_s_evolution"]:
            raise ValueError("native alpha_s evolution convention is missing")
        correction = record["prefactor_correction"]
        _equal_number(correction["upstream_prefactor_pb"], 35.0309, "upstream prefactor")
        _equal_number(correction["gev_minus2_to_pb"], 389379660., "POWHEG pb conversion")
        multiplier = GF * math.pi / (math.sqrt(2.) * 288.) * 389379660. / 35.0309
        _equal_number(correction["multiplier"], multiplier, "GF prefactor correction")
        for run in runs:
            _equal_number(run["raw_cross_section_pb"] * multiplier, run["cross_section_pb"], "GF-corrected run")
            _equal_number(run["raw_cross_section_error_pb"] * multiplier, run["cross_section_error_pb"], "GF-corrected error")
        if canonical_sha256(provenance["settings"]) != provenance["settings_sha256"]:
            raise ValueError("calculation settings hash does not match the saved settings")
    except (KeyError, TypeError) as exc:
        raise ValueError(f"incomplete ihixs normalization record: {exc}") from exc
    return record


def load_record(path):
    path = Path(path)
    try:
        record = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read ihixs normalization record {path}; run the ihixs calculation first") from exc
    validate_record(record)
    # Raw files are checked when carried with the record. Embedded sidecars
    # remain self-contained when moved to another host without the work tree.
    artifact_root = path.parent / record["provenance"].get("artifact_root", ".")
    if artifact_root.is_dir():
        for artifact in record["provenance"]["artifacts"]:
            candidate = artifact_root / artifact["path"]
            if candidate.is_file() and hashlib.sha256(candidate.read_bytes()).hexdigest() != artifact["sha256"]:
                raise ValueError(f"ihixs artifact hash changed: {candidate}")
    return record


def is_ho_signal(manifest):
    """Recognize HO campaign identity; incompatible physics is rejected later."""
    return (isinstance(manifest, dict) and manifest.get("sample") == "signal_gg_h_aa"
            and manifest.get("matching") == "POWHEG"
            and "HJMiNNLO" in str(manifest.get("hard_accuracy", ""))
            and "HJMiNNLO" in str(manifest.get("process", "")))


def validate_manifest(record, manifest):
    if not is_ho_signal(manifest):
        raise ValueError("ihixs normalization is only valid for the HO HJMiNNLO diphoton signal")
    _equal_number(manifest.get("ebeam_gev"), record["sqrt_s_gev"] / 2., "campaign beam energy")
    if manifest.get("pdf") != NATIVE_PDF_SET or manifest.get("lhaid") != NATIVE_PDF_ID:
        raise ValueError("HO signal native PDF differs from the normalization reference")
    for key in ("hard_model", "higgs_mass_gev", "top_mass_gev", "gf_gev_minus2"):
        if key in manifest and manifest[key] != record[key]:
            raise ValueError(f"HO signal {key} differs from the normalization reference")
    return record


def load_signal_record(path, manifest):
    return validate_manifest(load_record(path), manifest)


def signal_sidecar(record_path, manifest, native_lhe, weight_scale):
    record = load_signal_record(record_path, manifest)
    finite_number(weight_scale, "signal branching fraction", positive=True)
    finite_number(native_lhe["cross_section_pb"], "native cross section", positive=True)
    finite_number(native_lhe["cross_section_error_pb"], "native error", nonnegative=True)
    result = {"run_tag": manifest["run_tag"], "sample": manifest["sample"],
              "normalization_kind": "ihixs_n3lo", "cross_section_pb": record["cross_section_pb"],
              "cross_section_error_pb": record["cross_section_error_pb"],
              "weight_scale": weight_scale,
              "source": "ihixs pure-HEFT N3LO inclusive ggF; before detector response and BR",
              "native_cross_section_pb": native_lhe["cross_section_pb"],
              "native_cross_section_error_pb": native_lhe["cross_section_error_pb"],
              "native_lhe": copy.deepcopy(native_lhe),
              "ihixs_record_sha256": record["fingerprint"], "ihixs": copy.deepcopy(record)}
    if "shower_completion" in manifest:
        result["shower_completion_sha256"] = manifest["shower_completion"]["fingerprint"]
    return validate_sidecar(result, manifest, weight_scale)


def validate_sidecar(sidecar, manifest, weight_scale=None):
    """Return a verified sidecar, using its embedded record after migration."""
    try:
        if sidecar["normalization_kind"] != "ihixs_n3lo":
            raise ValueError("HO signal sidecar requires ihixs_n3lo normalization")
        record = validate_manifest(validate_record(sidecar["ihixs"]), manifest)
        if sidecar["ihixs_record_sha256"] != record["fingerprint"]:
            raise ValueError("sidecar record digest differs from its embedded record")
        if sidecar["sample"] != manifest["sample"] or sidecar["run_tag"] != manifest["run_tag"]:
            raise ValueError("sidecar campaign identity does not match")
        if "shower_completion" in manifest and sidecar.get("shower_completion_sha256") != manifest["shower_completion"]["fingerprint"]:
            raise ValueError("sidecar shower population differs from the campaign")
        if "shower_completion" not in manifest and "shower_completion_sha256" in sidecar:
            raise ValueError("sidecar shower population lacks its campaign completion record")
        _equal_number(sidecar["cross_section_pb"], record["cross_section_pb"], "sidecar cross section")
        _equal_number(sidecar["cross_section_error_pb"], record["cross_section_error_pb"], "sidecar error")
        br = finite_number(sidecar["weight_scale"], "sidecar BR", positive=True)
        if br > 1:
            raise ValueError("signal branching fraction exceeds one")
        _equal_number(manifest["weight_scale"], br, "manifest BR")
        if weight_scale is not None:
            _equal_number(weight_scale, br, "analysis BR")
        finite_number(sidecar["native_cross_section_pb"], "native cross section", positive=True)
        finite_number(sidecar["native_cross_section_error_pb"], "native numerical error", nonnegative=True)
        if sidecar["native_lhe"]["cross_section_pb"] != sidecar["native_cross_section_pb"]:
            raise ValueError("native LHE diagnostic cross section is inconsistent")
        if sidecar["native_lhe"]["cross_section_error_pb"] != sidecar["native_cross_section_error_pb"]:
            raise ValueError("native LHE diagnostic error is inconsistent")
    except (KeyError, TypeError) as exc:
        raise ValueError(f"incomplete HO signal normalization sidecar: {exc}") from exc
    return sidecar
