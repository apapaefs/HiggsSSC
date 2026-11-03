"""Audited, private ancestry showers of an existing completed LHE population."""
from pathlib import Path
import json
import subprocess
import sys

from hfourlepton import reprocess as revision
from hgammagamma import ho_shower_completion as completion
from hgammagamma.run_gammagamma_ho_campaign import ROOT_METADATA_SCRIPT


def plugin_artifacts(directory):
    directory = Path(directory).resolve()
    proof = json.loads((directory / "build.json").read_text())
    library = directory / "HwSimAncestry.so"
    revision.verify(library, proof["plugin_sha256"])
    if proof.get("schema_version") != 1:
        raise ValueError("unsupported ancestry plugin schema")
    for name, sha in proof["source_sha256"].items():
        revision.verify(directory / name, sha)
    return {"library": revision.artifact(library), "build": revision.artifact(directory / "build.json")}


def ancestry_card(original, library):
    if original.count("library HwSim.so") != 1 or original.count("set /Herwig/Analysis/HwSim:SaveObjects Yes") != 1:
        raise ValueError("expected one explicit HwSim library and SaveObjects Yes in parent card")
    if "SaveAncestry" in original:
        raise ValueError("parent already has ancestry settings; use a fresh original campaign")
    return original.replace("library HwSim.so", "library " + str(library)).replace(
        "set /Herwig/Analysis/HwSim:SaveObjects Yes",
        "set /Herwig/Analysis/HwSim:SaveObjects Yes\nset /Herwig/Analysis/HwSim:SaveAncestry Yes")


def logged(command, folder, name):
    with (folder / (name + ".log")).open("w") as stream:
        return subprocess.run([str(v) for v in command], cwd=folder,
                              stdout=stream, stderr=subprocess.STDOUT).returncode


def reconstruct(sample, before, raw, output, executable, plugin):
    """Repeat the same seeded shower and require unchanged count/weight moments.

    No hard-event generation, decay, source population, or normalization is
    changed. Failed sample attempts stay available; successful showers can be
    reused if just the analyzer step fails.
    """
    name = sample["name"]
    library = Path(plugin["library"]["path"])
    revision.verify(library, plugin["library"]["sha256"])
    lhe = Path(sample["input_lhe"])
    revision.verify(lhe, sample["input_lhe_sha256"])
    source_card = sample["cards"]["herwig"]
    revision.verify(source_card["path"], source_card["sha256"])
    card = ancestry_card(Path(source_card["path"]).read_text(), library)
    parent = output / "showers" / name
    completed = parent / "completed.json"
    if completed.exists():
        saved = json.loads(completed.read_text())
        folder = Path(saved["herwig_dir"])
        if (folder / (name + ".in")).read_text() != card:
            raise ValueError("ancestry shower card changed")
    else:
        attempt = revision.attempt_directory(parent)
        folder = attempt / "herwig"
        (folder / "events").mkdir(parents=True)
        (folder / (name + ".in")).write_text(card)
        count = sample["shower_completion"]["lhe_events"]
        requested = sample["nevents_requested"]
        plan = completion.finite_shower_plan(name, requested, count)
        revision.save(folder / "shower_execution.json", plan)
        print("Showering with ancestry: " + name, flush=True)
        if logged(["Herwig", "read", name + ".in", "-L", library.parent], folder, "read"):
            raise RuntimeError(f"Herwig read failed: {folder / 'read.log'}")
        code = logged(["Herwig", "run", name + ".run", "-L", library.parent,
                       "-N" + str(plan["execution_events"])], folder, "run")
        roots = sorted((folder / "events").glob("*.root"))
        if logged([sys.executable, "-c", ROOT_METADATA_SCRIPT, *roots], folder, "root-metadata"):
            raise RuntimeError(f"cannot audit shower ROOT files: {folder}")
        metadata = json.loads((folder / "root-metadata.log").read_text().splitlines()[-1])
        saved = completion.build_completion(folder, name, requested, count, metadata,
                                            "source_exhausted" if code else "requested_events")
        saved["source_weights"] = completion.source_weight_summary(lhe, saved["attempted_events"])
        saved["input_lhe_sha256"] = sample["input_lhe_sha256"]
        saved["fingerprint"] = completion.completion_fingerprint(saved)
    completion.validate_completion(saved, {"sample": name, "nevents_requested": sample["nevents_requested"],
        "lhe_file": str(lhe), "lhe": {"events": saved["lhe_events"]}}, directory=folder.parent)
    old = sample["shower_completion"]
    for field in ("saved_events", "attempted_events", "generated_events", "source_weights", "root_summary"):
        if saved[field] != old[field]:
            raise ValueError(f"metadata-only shower changed {field} for {name}; investigate before publishing")
    if not completed.exists():
        revision.save(completed, saved)
    inputs = parent / "roots.input"
    inventory = "".join(item["path"] + "\n" for item in saved["artifacts"]["roots"])
    if inputs.exists() and inputs.read_text() != inventory:
        raise ValueError("ancestry shower inventory changed")
    if not inputs.exists():
        inputs.write_text(inventory)
    analysis = revision.attempt_directory(output / "analysis/perfect" / name)
    command = revision.analysis_command(sample, {**before, "input": str(inputs)}, raw, analysis, executable)
    command.append("--cms-ancestry")
    if revision.run_logged(command, analysis / "analyzer.log"):
        raise RuntimeError(f"ancestry analysis failed: {analysis / 'analyzer.log'}")
    summary = analysis / f"{name}_{raw['run_tag']}_perfect.summary.json"
    entry = revision.updated_sample(sample, before, summary, command,
                                    selection="cms_fiducial_v3", input_path=str(inputs))
    entry["ancestry_schema_version"] = 1
    entry["shower_completion"] = saved
    entry["cards"]["herwig"] = revision.artifact(folder / (name + ".in"))
    entry["ancestry_plugin"] = plugin
    return entry
