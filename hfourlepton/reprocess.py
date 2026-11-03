#!/usr/bin/env python3
"""Reconstruct saved CMS-profile showers in an independent, resumable revision.

By default this reads saved showers. With --ancestry-plugin-dir it repeats the
same seeded Herwig showers in private directories; hard events and decays are
unchanged. It never invokes MG5, MadSpin, or the campaign runner.
The original manifest and all original physics products are read-only inputs.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import fcntl
import hashlib
import html
import importlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from hfourlepton import make_fourlepton_report as report

CODE = REPO / "hfourlepton/LOAnalysis/Code"
TARGET = "HwSimPostAnalysis_fourlepton"
SELECTION = "cms_fiducial_v2"
SELECTION_FIELDS = ("selection_version", "photon_dressing_dr", "alternative_pair_veto")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def save(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def verify(path, expected):
    path = Path(path)
    if not expected or not path.is_file() or digest(path) != expected:
        raise ValueError(f"missing or changed artifact: {path}")


def artifact(path):
    return {"path": str(path), "sha256": digest(path), "exists": True,
            "checksum_status": "available"}


def audit_parent(path):
    raw = json.loads(path.read_text())
    if (raw.get("order") != "ho" or raw.get("selection_profile") != "cms"
            or raw.get("detector_profile") != "perfect"
            or raw.get("collider", {}).get("sqrt_s_tev") != 13.):
        raise ValueError("reprocessing currently supports completed 13 TeV HO CMS/perfect manifests")
    bundle = report.load_manifest(path)
    if bundle.audit.errors:
        raise ValueError("; ".join(bundle.audit.errors))
    summaries = {}
    for sample in bundle.samples:
        s = sample.raw
        if s.get("status") != "complete" or not re.fullmatch(r"[A-Za-z0-9_]+", sample.name):
            raise ValueError(f"incomplete or invalid sample: {sample.name}")
        verify(sample.output_root, s.get("output_root_sha256"))
        verify(s["analysis_summary"], s.get("analysis_summary_sha256"))
        completed = report._validated_shower_completion(sample, verify_files=True)
        summary = json.loads(Path(s["analysis_summary"]).read_text())
        roots = [str(Path(item["path"]).resolve()) for item in completed["artifacts"]["roots"]]
        listed = [str(Path(line.strip()).resolve()) for line in Path(summary["input"]).read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#")]
        if listed != roots or len(set(listed)) != len(listed):
            raise ValueError(f"shower input inventory differs for {sample.name}")
        if summary["events_read"] != completed["saved_events"]:
            raise ValueError(f"analysis/shower event count differs for {sample.name}")
        for key, expected in report._expected_runtime_metadata(sample).items():
            if key in ("sample", "category", "requested_channel", "seed", "weight_scale", "selection_profile", "response_profile"):
                if not report._runtime_values_equal(key, expected, summary.get(key)):
                    raise ValueError(f"parent {key} differs for {sample.name}")
        optional = s.get("optional_weight_names_file")
        if optional:
            verify(optional["path"], optional["sha256"])
            if summary.get("optional_weight_names_file") != optional["path"]:
                raise ValueError(f"optional-weight mapping differs for {sample.name}")
        summaries[sample.name] = summary
    if not summaries:
        raise ValueError("manifest has no completed samples")
    return raw, summaries


def run_logged(command, log):
    print("Running " + log.stem + "; log: " + str(log), flush=True)
    with log.open("a") as stream:
        stream.write("COMMAND: " + json.dumps([str(v) for v in command]) + "\n")
        stream.flush()
        result = subprocess.run([str(v) for v in command], stdout=stream, stderr=subprocess.STDOUT)
    return result.returncode


def attempt_directory(parent):
    parent.mkdir(parents=True, exist_ok=True)
    i = 1
    while (parent / f"attempt_{i:03d}").exists():
        i += 1
    path = parent / f"attempt_{i:03d}"
    path.mkdir()
    return path


def analysis_command(s, summary, raw, output, executable):
    command = [executable, "--input-list", summary["input"], "--response-profile", "perfect",
               "--selection-profile", "cms", "--seed", str(summary["seed"]),
               "--weight-scale", str(summary["weight_scale"]), "--muon-resolution-scale",
               str(summary["muon_resolution_scale"]), "--tag", raw["run_tag"],
               "--output-dir", output, "--sample", s["name"], "--category", s["category"],
               "--channel", s["channel"]]
    if s.get("optional_weight_names_file"):
        command += ["--optional-weight-names-file", s["optional_weight_names_file"]["path"]]
    return [str(v) for v in command]


def updated_sample(s, before, summary_path, command, *, selection=SELECTION, input_path=None):
    after = json.loads(summary_path.read_text())
    if (after.get("selection_version") != selection or after.get("photon_dressing_dr") != .3
            or after.get("alternative_pair_veto") is not True):
        raise ValueError("analyzer did not execute the corrected CMS selection")
    if selection == "cms_fiducial_v3" and after.get("ancestry_schema_version") != 1:
        raise ValueError("CMS v3 analysis did not verify ancestry schema 1")
    for key in ("events_read", "sample", "category", "requested_channel", "seed", "weight_scale",
                "muon_resolution_scale", "response_profile", "selection_profile", "input",
                "optional_weight_names_file"):
        expected = input_path if key == "input" and input_path is not None else before.get(key)
        if after.get(key) != expected:
            raise ValueError(f"reprocessing changed {key} for {s['name']}")
    for key in ("generated_raw", "generated_scaled"):
        for field, value in before[key].items():
            if not math.isclose(after[key][field], value, rel_tol=1.e-10, abs_tol=1.e-12):
                raise ValueError(f"reprocessing changed source weights: {s['name']} {key}.{field}")
    entry = copy.deepcopy(s)
    entry.update(output_root=after["output_root"], root_file=after["output_root"],
                 output_root_sha256=digest(after["output_root"]), output_root_checksum_status="available",
                 analysis_summary=str(summary_path), analysis_summary_sha256=digest(summary_path),
                 analyzer_cli=command, cut_mask_definition=after["cut_mask_definition"],
                 weight_sums={k: after[k] for k in ("generated_raw", "generated_scaled", "selected_scaled")},
                 **{key: after[key] for key in SELECTION_FIELDS})
    return entry


def reconstruct(s, before, raw, output, executable):
    folder = attempt_directory(output / "analysis/perfect" / s["name"])
    command = analysis_command(s, before, raw, folder, executable)
    log = folder / "analyzer.log"
    if run_logged(command, log):
        raise RuntimeError(f"sample {s['name']} failed; inspect {log}; rerun to resume")
    summary = folder / f"{s['name']}_{raw['run_tag']}_perfect.summary.json"
    return updated_sample(s, before, summary, command)


def stage(output, name, identity, command, *, validation=False):
    """Publish complete downstream results atomically; preserve failed attempts."""
    final = output / name
    if final.exists():
        proof = json.loads((final / "completion.json").read_text())
        if proof["identity"] != identity:
            raise ValueError(f"different analysis settings in {final}; use a new revision")
        for relative, checksum in proof["files"].items():
            verify(final / relative, checksum)
        return proof
    work = attempt_directory(output / "attempts" / name)
    log = work.parent / (work.name + ".log")
    code = run_logged(command(work), log)
    result = json.loads((work / "comparison.json").read_text()) if validation and code in (0, 2) else {}
    if validation and result.get("status") not in ("not_validated", "pass_approximate_acceptance_check"):
        raise ValueError(f"validation did not produce a valid comparison: {work}")
    if code and not (validation and code == 2 and result["status"] == "not_validated"):
        raise RuntimeError(f"{name} failed (exit {code}); inspect {log}; rerun to resume")
    required = ("index.html", "comparison.json") if validation else ("index.html", "report.json") if name == "report" else ("index.html", "summary.json")
    if not all((work / p).is_file() for p in required):
        raise ValueError(f"incomplete {name} output: {work}")
    proof = {"identity": identity, "result": result.get("status", "complete"),
             "files": {str(p.relative_to(work)): digest(p) for p in sorted(work.rglob("*")) if p.is_file()}}
    save(work / "completion.json", proof)
    work.rename(final)
    return proof


def execute(parent, output, jobs, luminosity, *, check=False, ancestry_plugin_dir=None):
    raw, summaries = audit_parent(parent)
    if parent.is_relative_to(output):
        raise ValueError("revision directory must not contain the original manifest")
    sources = [CODE / (TARGET + ".cc"), CODE / "Makefile", CODE / "HwSimAncestry.h", Path(__file__).resolve(),
               REPO / "hfourlepton/make_fourlepton_report.py", REPO / "hfourlepton/validate_lhc.py",
               REPO / "hfourlepton/analyze_fourlepton.py", REPO / "ho_xgboost_core.py",
               *[REPO / f"hfourlepton/analysis_cards/lhc13_{name}.json" for name in ("cuts", "xgboost")]]
    plugin = None
    if ancestry_plugin_dir:
        from hfourlepton import ancestry_shower
        plugin = ancestry_shower.plugin_artifacts(ancestry_plugin_dir)
        sources += [Path(ancestry_shower.__file__), REPO / "hgammagamma/ho_shower_completion.py"]
    identity = {"parent": artifact(parent), "selection_version": "cms_fiducial_v3" if plugin else SELECTION,
                "source_hashes": {str(p.relative_to(REPO)): digest(p) for p in sources},
                "luminosity_fb": luminosity}
    if plugin:
        identity["ancestry_plugin"] = plugin
    if check:
        print(f"Verified {len(summaries)} samples and saved showers. Planned revision: {output}")
        return
    for module in ("ROOT", "numpy", "matplotlib", "xgboost", "sklearn"):
        importlib.import_module(module)
    if output.exists() and any(output.iterdir()) and not (output / "state.json").is_file():
        raise ValueError(f"unmanaged nonempty revision directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("this revision is already running") from exc
        state_path = output / "state.json"
        state = json.loads(state_path.read_text()) if state_path.exists() else {"identity": identity, "samples": {}}
        if state["identity"] != identity:
            raise ValueError("revision inputs or software changed; use a new revision name")
        save(state_path, state)
        code = output / "code"
        code.mkdir(exist_ok=True)
        for src in sources[:3]:
            target = code / src.name
            if target.exists():
                verify(target, digest(src))
            else:
                shutil.copyfile(src, target)
        executable = code / TARGET
        if state.get("executable_sha256"):
            verify(executable, state["executable_sha256"])
        else:
            if run_logged(["make", "-B", "-C", code, TARGET], output / "build.log"):
                raise RuntimeError(f"analyzer build failed; inspect {output / 'build.log'}")
            state["executable_sha256"] = digest(executable)
            save(state_path, state)
        for s in state["samples"].values():
            verify(s["output_root"], s["output_root_sha256"])
            verify(s["analysis_summary"], s["analysis_summary_sha256"])
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(ancestry_shower.reconstruct if plugin else reconstruct,
                                   s, summaries[s["name"]], raw, output, executable, *([plugin] if plugin else [])): s["name"]
                       for s in raw["samples"] if s["name"] not in state["samples"]}
            failures = []
            for future in as_completed(futures):
                try:
                    entry = future.result()
                except Exception as exc:
                    failures.append(str(exc))
                    continue
                state["samples"][entry["name"]] = entry
                save(state_path, state)
                print(f"Completed {entry['name']} ({len(state['samples'])}/{len(raw['samples'])})", flush=True)
            if failures:
                raise RuntimeError("\n".join(failures))
        manifest = output / "manifest_perfect.json"
        verify(parent, identity["parent"]["sha256"])
        if state.get("manifest_sha256"):
            verify(manifest, state["manifest_sha256"])
        else:
            revised = copy.deepcopy(raw)
            revised["samples"] = [state["samples"][s["name"]] for s in raw["samples"]]
            revised["analysis_revision"] = {**identity, "name": output.name,
                "parent_configuration_fingerprint": raw.get("configuration_fingerprint")}
            revised["configuration_fingerprint"] = fingerprint(identity)
            revised["software"]["analyzer_executable"] = artifact(executable)
            revised["software"]["analysis_reprocessor"] = artifact(Path(__file__).resolve())
            revised["tracked_inputs"].update(analyzer_source=artifact(code / (TARGET + ".cc")),
                                             analyzer_makefile=artifact(code / "Makefile"),
                                             analyzer_ancestry_header=artifact(code / "HwSimAncestry.h"))
            for cut in revised["cut_stages"]:
                if cut["bit"] == 5:
                    cut["label"] += "; alternative-pair veto"
            save(manifest, revised)
            state["manifest_sha256"] = digest(manifest)
            save(state_path, state)
        stage_id = fingerprint({"identity": identity, "manifest": state["manifest_sha256"]})
        py = sys.executable
        stage(output, "report", stage_id, lambda work: [py, REPO / "hfourlepton/make_fourlepton_report.py",
              manifest, "--output-dir", work, "--luminosity-fb", luminosity, "--strict"])
        validation = stage(output, "lhc_validation", stage_id, lambda work: [py, REPO / "hfourlepton/validate_lhc.py",
                           "--manifest", manifest, "--output-dir", work], validation=True)
        cards = output / "analysis_cards"
        cards.mkdir(exist_ok=True)
        for name in ("cuts", "xgboost"):
            card = json.loads((REPO / f"hfourlepton/analysis_cards/lhc13_{name}.json").read_text())
            if (raw.get("sampling") or {}).get("preset") == "minimal":
                card["min_background_effective_events"] = 400. if name == "cuts" else 100.
            path = cards / (name + ".json")
            if path.exists() and json.loads(path.read_text()) != card:
                raise ValueError(f"analysis card changed: {path}")
            if not path.exists():
                save(path, card)
            stage(output, name, stage_id, lambda work: [py, REPO / "hfourlepton/analyze_fourlepton.py", name,
                  "--manifest", manifest, "--config", path, "--output-dir", work, "--luminosity-fb", luminosity])
        state.update(state="complete", validation_status=validation["result"])
        save(state_path, state)
        links = ''.join(f'<li><a href="{name}/index.html">{name}</a></li>'
                        for name in ("report", "lhc_validation", "cuts", "xgboost"))
        (output / "index.html").write_text('<!doctype html><meta charset="utf-8"><title>CMS analysis revision</title>'
            f'<h1>{html.escape(output.name)}</h1><p>Analysis complete. Acceptance check: '
            f'{html.escape(validation["result"])}</p><ul>{links}</ul>')
        print(f"Reprocessing complete. Acceptance check: {validation['result']}. Results: {output / 'index.html'}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--luminosity-fb", type=float, default=100.)
    p.add_argument("--ancestry-plugin-dir", type=Path, help="repeat seeded showers using a private schema-1 HwSim plugin and apply CMS v3")
    p.add_argument("--check", action="store_true", help="audit saved inputs without writing files or executing analysis")
    args = p.parse_args(argv)
    if args.jobs < 1 or not math.isfinite(args.luminosity_fb) or args.luminosity_fb <= 0:
        p.error("jobs and luminosity must be positive and finite")
    execute(args.manifest.resolve(), args.output_dir.resolve(), args.jobs, args.luminosity_fb, check=args.check,
            ancestry_plugin_dir=args.ancestry_plugin_dir)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        sys.exit(f"Error: {exc}")
