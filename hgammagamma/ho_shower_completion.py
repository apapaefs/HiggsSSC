"""Evidence for completed HO showers and their signed source-weight population.

These helpers only read persisted artifacts. An exhausted finite LHE stream is
accepted only when Herwig finalized its outputs and accounts for every source
attempt. Failed shower events have zero simulated response; their signed source
weights remain in the denominator instead of renormalizing away their loss.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import xml.etree.ElementTree as ET


EOF_MESSAGE = "More events requested than available in LesHouchesReader"
KINDS = {"requested_events", "source_exhausted"}
SEVERITIES = ("info", "warning", "eventerror", "runerror", "setuperror",
              "maybeabort", "abortnow", "unknown")
MOMENTS = ("sum_weight", "sum_abs_weight", "sum_weight_squared")
REL_TOL = 1.e-7


def completion_fingerprint(record):
    """Canonical digest; call again if adding source_weights to an existing record."""
    payload = {key: value for key, value in record.items() if key != "fingerprint"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _integer(value, name, *, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _finite(value, name):
    if isinstance(value, bool):
        raise ValueError(f"non-finite {name}")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"non-finite {name}") from error
    if not math.isfinite(result):
        raise ValueError(f"non-finite {name}")
    return result


def _number(text):
    return _finite(str(text).replace("D", "E").replace("d", "e"), "LHE value")


def _lhe_integer(text, name):
    value = _number(text)
    if not value.is_integer():
        raise ValueError(f"noninteger LHE {name}")
    return int(value)


def _close(first, second):
    return math.isclose(first, second, rel_tol=REL_TOL, abs_tol=1.e-12)


def _artifact(path):
    path = Path(path).resolve()
    before = path.stat()
    if not path.is_file() or before.st_size <= 0:
        raise ValueError(f"missing or empty shower artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"shower artifact changed while reading: {path}")
    return {"path": str(path), "size": after.st_size, "mtime_ns": after.st_mtime_ns,
            "sha256": digest.hexdigest()}


def _validate_artifact(record, *, verify_file=False):
    if not isinstance(record, dict) or not isinstance(record.get("path"), str):
        raise ValueError("invalid shower artifact provenance")
    if not Path(record["path"]).is_absolute():
        raise ValueError("shower artifact path must be absolute")
    _integer(record.get("size"), "artifact size", minimum=1)
    _integer(record.get("mtime_ns"), "artifact mtime_ns")
    if not re.fullmatch(r"[0-9a-f]{64}", str(record.get("sha256", ""))):
        raise ValueError("invalid shower artifact SHA-256")
    if verify_file and _artifact(record["path"]) != record:
        raise ValueError(f"shower artifact changed: {record['path']}")


def _moments(record, *, positive=False):
    values = {key: _finite(record.get(key), key) for key in MOMENTS}
    if values["sum_abs_weight"] < 0 or values["sum_weight_squared"] <= 0:
        raise ValueError("invalid signed-weight moments")
    tolerance = max(1.e-12, REL_TOL * values["sum_abs_weight"])
    if values["sum_abs_weight"] + tolerance < abs(values["sum_weight"]):
        raise ValueError("sum_abs_weight is smaller than abs(sum_weight)")
    if positive and values["sum_weight"] <= 0:
        raise ValueError("non-positive inclusive signed weight population")
    return values


def _root_summary(metadata):
    if not isinstance(metadata, list) or not metadata:
        raise ValueError("missing finalized ROOT metadata")
    paths = set()
    for entry in metadata:
        if not isinstance(entry, dict) or not isinstance(entry.get("path"), str):
            raise ValueError("invalid ROOT metadata")
        path = str(Path(entry["path"]).resolve())
        if not Path(entry["path"]).is_absolute() or path in paths:
            raise ValueError("duplicate or relative ROOT path")
        paths.add(path)
        _integer(entry.get("entries"), "ROOT entries")
        _integer(entry.get("size"), "ROOT size", minimum=1)
        _integer(entry.get("mtime_ns"), "ROOT mtime_ns")
        if entry["entries"]:
            _moments(entry)
        elif any(_finite(entry.get(key), key) != 0 for key in MOMENTS):
            raise ValueError("empty ROOT tree has nonzero weight moments")
    result = {"entries": sum(entry["entries"] for entry in metadata)}
    result.update({key: math.fsum(_finite(entry[key], key) for entry in metadata)
                   for key in MOMENTS})
    if result["entries"] <= 0:
        raise ValueError("no saved ROOT events")
    _moments(result, positive=True)
    return result


def _exception_summary(text):
    marker = "The following exception classes were reported in this run:"
    no_exceptions = "No exceptions reported in this run."
    counts = dict.fromkeys(SEVERITIES, 0)
    if text.count(marker) + text.count(no_exceptions) != 1:
        raise ValueError("missing or repeated final Herwig exception summary")
    classes = []
    if marker in text:
        tail = text.split(marker, 1)[1].lstrip("\r\n")
        body = re.split(r"\r?\n\s*\r?\n", tail, maxsplit=1)[0]
        for line in body.splitlines():
            match = re.fullmatch(r"\s*(\S+)\s+(\w+)\s+\((\d+) times\)\s*", line)
            if not match or match[2] not in counts or int(match[3]) <= 0:
                raise ValueError("malformed final Herwig exception summary")
            classes.append({"class": match[1], "severity": match[2], "count": int(match[3])})
            counts[match[2]] += int(match[3])
        if not classes:
            raise ValueError("empty final Herwig exception summary")
    if any(counts[key] for key in ("setuperror", "maybeabort", "abortnow", "unknown")):
        raise ValueError("fatal Herwig exception prevents shower completion")
    return counts, classes


def _parse_logs(herwig_dir, sample_name, termination):
    paths = {"run_log": herwig_dir / "run.log", "herwig_out": herwig_dir / f"{sample_name}.out",
             "herwig_log": herwig_dir / f"{sample_name}.log"}
    artifacts = {key: _artifact(path) for key, path in paths.items()}
    texts = {key: path.read_text(errors="strict").replace("\r", "\n")
             for key, path in paths.items()}
    out = texts["herwig_out"]
    if out.count("Statistics for Les Houches event handler") != 1:
        raise ValueError("missing or repeated Les Houches statistics")
    section = out.split("Statistics for Les Houches event handler", 1)[1]
    section = section.split("Per Les Houches Reader", 1)[0]
    totals = re.findall(r"(?m)^Total:\s+(\S+)\s+(\S+)\s+([^\n]+)$", section)
    if len(totals) != 1:
        raise ValueError("missing or ambiguous final Les Houches Total")
    # ThePEG switches large counts to scientific notation (e.g. 1e+06).
    # Still require finite integer counts; completion validation below checks
    # them against the exact ROOT, LHE and discarded-event populations.
    generated = _lhe_integer(totals[0][0], "generated event count")
    attempted = _lhe_integer(totals[0][1], "attempted event count")
    # Herwig writes uncertainties as e.g. 0.214(1)e+00.
    cross_section = re.sub(r"\([^)]*\)", "", totals[0][2].strip()).split()[0]
    _number(cross_section)
    run = texts["run_log"]
    footer = re.findall(r"(?m)^Number of events that pass basic cuts:\s*(\d+)\s*$", run)
    weights = re.findall(r"(?m)^Weight of events that pass basic cuts:\s*(\S+)\s*$", run)
    writers = re.findall(r"(?m)^A root tree has been written to the file:\s*(\S+)\s*$", run)
    if len(footer) != 1 or len(weights) != 1 or not writers or len(set(writers)) != len(writers):
        raise ValueError("missing or repeated finalized HwSim output/footer")
    _number(weights[0])
    if run.index("A root tree has been written to the file:") > run.index("Number of events that pass basic cuts:"):
        raise ValueError("HwSim footer precedes finalized ROOT output")
    counts, classes = _exception_summary(texts["herwig_log"])
    all_text = "\n".join(texts.values())
    if re.search(r"Reopening LesHouchesReader|Segmentation fault|segmentation fault|"
                 r"core dumped|Traceback \(most recent call last\)|terminate called|"
                 r"(?:^|\n)(?:Aborted|Killed)(?:\s|$)", all_text):
        raise ValueError("replayed input or unrelated process failure in shower logs")
    caught = re.findall(r"(?m)^Herwig:.*caught\.?\s*$", run)
    eof_occurrences = run.count(EOF_MESSAGE)
    if termination == "source_exhausted":
        if eof_occurrences != 1 or counts["runerror"] != 1 or len(caught) != 1:
            raise ValueError("source exhaustion requires one terminal EOF runerror")
        terminal = run.split(EOF_MESSAGE, 1)[1].strip()
        if not re.fullmatch(r"\S+\s+See logfile for details\.", terminal):
            raise ValueError("unexpected diagnostic after terminal Les Houches EOF")
    elif eof_occurrences or counts["runerror"] or caught:
        raise ValueError("failed Herwig run cannot be ordinary shower completion")
    violations = re.findall(r"BasicConsistency: maximum 4-momentum violation:\s*(\S+)\s+MeV", run)
    if len(violations) > 1:
        raise ValueError("repeated BasicConsistency finalization")
    maximum = _number(violations[0]) if violations else None
    if maximum is not None and maximum < 0:
        raise ValueError("negative maximum momentum violation")
    return {"generated_events": generated, "attempted_events": attempted,
            "saved_events": int(footer[0]), "exception_counts": counts,
            "exception_classes": classes, "max_momentum_violation_mev": maximum,
            "root_writers": [str((herwig_dir / path).resolve()) for path in writers],
            "artifacts": artifacts}


def _audit_source_card(path):
    """Guard the prefix and weight units used by source_weight_summary."""
    text = path.read_text()
    commands = [line.split("#", 1)[0].strip() for line in text.splitlines()]
    commands = [line for line in commands if line]
    required = {"WeightOption": "VarNegWeight", "AllowedToReOpen": "No",
                "OnTheFlyAnalysis": "No"}
    for key, value in required.items():
        settings = [line.split()[-1] for line in commands if re.search(rf":{key}\s", line)]
        if settings != [value]:
            raise ValueError(f"source weight proof requires {key} {value}")
    readers = [line for line in commands if line.startswith("create ThePEG::LesHouchesFileReader ")]
    reader_insertions = [line for line in commands if re.search(r":LesHouchesReaders\s", line)]
    if len(readers) != 1 or len(reader_insertions) != 1:
        raise ValueError("source weight proof requires exactly one LHE reader")
    reader = shlex.split(readers[0])[-1]
    filenames = [shlex.split(line)[2:] for line in commands
                 if line.startswith(f"set {reader}:FileName ")]
    if len(filenames) != 1 or len(filenames[0]) != 1:
        raise ValueError("source weight proof requires exactly one reader FileName")
    source_path = Path(filenames[0][0])
    if not source_path.is_absolute():
        source_path = path.parent / source_path
    for key, allowed in (("NormalizeWeights", {"No", "0"}), ("MaxScan", {"-1"}),
                         ("MaxFactor", {"1", "1.0"}), ("ReweightPDF", {"No", "0"}),
                         ("WeightNormalization", {"Normalized", "0"})):
        settings = [line.split()[-1] for line in commands if re.search(rf":{key}\s", line)]
        if any(value not in allowed for value in settings):
            raise ValueError(f"unsupported {key} for signed source weights")
    if any(re.search(r":(?:Reweights|Preweights|CKKWHandler|CacheFileName)\s", line) for line in commands):
        raise ValueError("reweighted or cached LHE input is unsupported")
    # The class declaration "create ThePEG::Cuts ..." is not a Cuts setting.
    cuts = [line.split()[-1] for line in commands if re.match(r"set\s+\S+:Cuts\s", line)]
    if len(cuts) != 2 or any(value != "/Herwig/Cuts/NoCuts" for value in cuts):
        raise ValueError("source weight proof requires NoCuts on handler and reader")
    return str(source_path.resolve())


def build_completion(herwig_dir, sample_name, requested_events, lhe_events,
                     root_metadata, termination, *, source_weights=None):
    herwig_dir = Path(herwig_dir).resolve()
    if termination not in KINDS:
        raise ValueError("unknown shower termination")
    _integer(requested_events, "requested_events", minimum=1)
    _integer(lhe_events, "lhe_events", minimum=1)
    parsed = _parse_logs(herwig_dir, sample_name, termination)
    metadata = copy.deepcopy(root_metadata)
    summary = _root_summary(metadata)
    root_artifacts = []
    for entry in metadata:
        path = Path(entry["path"]).resolve()
        if not path.is_relative_to(herwig_dir / "events"):
            raise ValueError("ROOT artifact lies outside this shower's events directory")
        artifact = _artifact(path)
        if (artifact["size"], artifact["mtime_ns"]) != (entry["size"], entry["mtime_ns"]):
            raise ValueError("ROOT artifact changed since metadata inspection")
        root_artifacts.append(artifact)
    if set(parsed.pop("root_writers")) != {entry["path"] for entry in root_artifacts}:
        raise ValueError("HwSim finalized ROOT paths differ from inspected ROOT inventory")
    parsed["artifacts"]["roots"] = root_artifacts
    card = herwig_dir / f"{sample_name}.in"
    reader_lhe_path = None
    if card.is_file():
        reader_lhe_path = _audit_source_card(card)
        parsed["artifacts"]["herwig_card"] = _artifact(card)
    elif source_weights is not None:
        raise ValueError("signed source weights require the persisted Herwig input card")
    record = {"schema_version": 1, "kind": "ho_shower_completion", "sample": sample_name,
              "herwig_dir": str(herwig_dir), "termination": termination,
              "requested_events": requested_events, "lhe_events": lhe_events,
              "discarded_events": parsed["attempted_events"] - parsed["generated_events"],
              "root_metadata": metadata, "root_summary": summary,
              "normalization_population": "consumed_signed_source_weights",
              "reader_lhe_path": reader_lhe_path,
              "failed_event_response": "zero_unsimulated_response", **parsed}
    if source_weights is not None:
        record["source_weights"] = copy.deepcopy(source_weights)
    record["fingerprint"] = completion_fingerprint(record)
    validate_completion(record, {"sample": sample_name, "nevents_requested": requested_events,
                                 "lhe": {"events": lhe_events}},
                        require_source=source_weights is not None)
    return record


def _data_lines(text):
    return [line.split("#", 1)[0].strip().split() for line in (text or "").splitlines()
            if line.split("#", 1)[0].strip()]


def source_weight_summary(lhe_path, attempts):
    """Stream a complete one-process LHE and normalize its consumed prefix.

    ThePEG's full MaxScan replaces XMAXUP with the largest absolute event weight
    for IDWTUP +/-3 and +/-4. VarNegWeight and NormalizeWeights No then deliver
    XWGTUP/max(abs(XWGTUP)); no additional BR belongs in these moments.
    """
    _integer(attempts, "attempts", minimum=1)
    path = Path(lhe_path).resolve()
    artifact = _artifact(path)
    opener = gzip.open if path.suffix == ".gz" else open
    count, maximum, prefix = 0, 0.0, []
    init_seen = False
    event_seen = False
    root = None
    idwtup = process_id = None
    try:
        with opener(path, "rb") as stream:
            for action, node in ET.iterparse(stream, events=("start", "end")):
                if root is None:
                    root = node
                    if action != "start" or root.tag != "LesHouchesEvents":
                        raise ValueError("missing LesHouchesEvents root")
                if action != "end":
                    continue
                if node.tag == "init":
                    if init_seen or event_seen:
                        raise ValueError("repeated or misplaced LHE init")
                    rows = _data_lines(node.text)
                    if len(rows) != 2 or len(rows[0]) != 10 or len(rows[1]) != 4:
                        raise ValueError("one complete LHE subprocess is required")
                    header = rows[0]
                    for index in (0, 1, 4, 5, 6, 7, 8, 9):
                        _lhe_integer(header[index], "init header")
                    if any(_number(header[index]) <= 0 for index in (2, 3)):
                        raise ValueError("invalid LHE beam energy")
                    idwtup = _lhe_integer(header[8], "IDWTUP")
                    if idwtup not in (-4, -3, 3, 4) or _lhe_integer(header[9], "NPRUP") != 1:
                        raise ValueError("source weights require IDWTUP +/-3 or +/-4 and one subprocess")
                    for value in rows[1][:3]:
                        _number(value)
                    process_id = _lhe_integer(rows[1][3], "LPRUP")
                    init_seen = True
                    node.clear()
                elif node.tag == "event":
                    if not init_seen:
                        raise ValueError("LHE event precedes init")
                    event_seen = True
                    rows = _data_lines(node.text)
                    if not rows or len(rows[0]) != 6:
                        raise ValueError("malformed LHE event header")
                    nup = _lhe_integer(rows[0][0], "NUP")
                    if nup <= 0 or len(rows) != nup + 1:
                        raise ValueError("LHE particle count does not match NUP")
                    if _lhe_integer(rows[0][1], "IDPRUP") != process_id:
                        raise ValueError("undeclared LHE subprocess")
                    weight = _number(rows[0][2])
                    for value in rows[0][3:]:
                        _number(value)
                    if weight == 0:
                        raise ValueError("zero LHE weights make source-attempt bookkeeping ambiguous")
                    for row in rows[1:]:
                        if len(row) != 13:
                            raise ValueError("malformed LHE particle row")
                        for value in row[:6]:
                            _lhe_integer(value, "particle header")
                        for index in (2, 3):
                            if not 0 <= _lhe_integer(row[index], "mother index") <= nup:
                                raise ValueError("invalid LHE mother index")
                        for value in row[6:]:
                            _number(value)
                    count += 1
                    maximum = max(maximum, abs(weight))
                    if count <= attempts:
                        prefix.append(weight)
                    node.clear()
                    root.clear()
        if not init_seen or count < attempts or root is None:
            raise ValueError("incomplete LHE source population")
    except (ET.ParseError, OSError, EOFError) as error:
        raise ValueError(f"malformed or unclosed LHE source: {path}") from error
    if _artifact(path) != artifact:
        raise ValueError("LHE changed while computing source weight proof")
    normalized = [weight / maximum for weight in prefix]
    result = {"schema_version": 1, "kind": "ho_source_weight_summary", "lhe_events": count,
              "attempted_events": attempts, "idwtup": idwtup, "process_id": process_id,
              "effective_max_weight": maximum, "normalization": "XWGTUP/full_file_max_abs_XWGTUP",
              "sampling": "first_attempted_events_without_replay", "zero_weight_events": 0,
              "negative_events": sum(weight < 0 for weight in normalized),
              "sum_weight": math.fsum(normalized),
              "sum_abs_weight": math.fsum(abs(weight) for weight in normalized),
              "sum_weight_squared": math.fsum(weight * weight for weight in normalized),
              "artifact": artifact}
    _moments(result, positive=True)
    result["fingerprint"] = completion_fingerprint(result)
    return result


def _validate_source(record, root_summary, *, require_source=True):
    source = record.get("source_weights")
    if source is None:
        if require_source and record["discarded_events"]:
            raise ValueError("discarded shower events require signed source-weight proof")
        return
    if (not isinstance(source, dict) or source.get("kind") != "ho_source_weight_summary" or source.get("schema_version") != 1
            or source.get("fingerprint") != completion_fingerprint(source)):
        raise ValueError("invalid signed source-weight provenance")
    _integer(source.get("lhe_events"), "source lhe_events", minimum=1)
    _integer(source.get("attempted_events"), "source attempted_events", minimum=1)
    _integer(source.get("negative_events"), "source negative_events")
    _integer(source.get("zero_weight_events"), "source zero_weight_events")
    if source["negative_events"] > source["attempted_events"]:
        raise ValueError("source negative-weight count exceeds attempted events")
    if (source.get("lhe_events") != record["lhe_events"]
            or source.get("attempted_events") != record["attempted_events"]
            or source.get("zero_weight_events") != 0 or source.get("idwtup") not in (-4, -3, 3, 4)
            or source.get("normalization") != "XWGTUP/full_file_max_abs_XWGTUP"
            or source.get("sampling") != "first_attempted_events_without_replay"):
        raise ValueError("source-weight population differs from shower bookkeeping")
    if _finite(source.get("effective_max_weight"), "effective_max_weight") <= 0:
        raise ValueError("invalid effective source maximum weight")
    _validate_artifact(source.get("artifact"))
    moments = _moments(source, positive=True)
    tolerance = max(1.e-12, REL_TOL * moments["sum_abs_weight"])
    lost_abs = moments["sum_abs_weight"] - root_summary["sum_abs_weight"]
    lost_signed = moments["sum_weight"] - root_summary["sum_weight"]
    lost_squared = moments["sum_weight_squared"] - root_summary["sum_weight_squared"]
    if lost_abs < -tolerance or abs(lost_signed) > lost_abs + tolerance or lost_squared < -tolerance:
        raise ValueError("saved/source signed-weight bookkeeping does not close")
    if not record["discarded_events"] and any(not _close(moments[key], root_summary[key]) for key in MOMENTS):
        raise ValueError("complete shower source/ROOT weights do not close")
    if record["discarded_events"]:
        discarded = record["discarded_events"]
        if lost_abs > discarded + tolerance or lost_squared > discarded + tolerance:
            raise ValueError("discarded source weights exceed normalized event bounds")
        if lost_abs * lost_abs > discarded * max(0.0, lost_squared) + tolerance:
            raise ValueError("discarded signed-weight moments are inconsistent")


def validate_completion(record, manifest, directory=None, *, require_source=True):
    """Check a sealed completion record; optionally rehash files in its sample directory."""
    if not isinstance(record, dict) or record.get("schema_version") != 1 or record.get("kind") != "ho_shower_completion":
        raise ValueError("invalid HO shower completion record")
    if record.get("fingerprint") != completion_fingerprint(record):
        raise ValueError("HO shower completion fingerprint changed")
    if record.get("termination") not in KINDS:
        raise ValueError("unknown shower termination")
    if (record.get("normalization_population") != "consumed_signed_source_weights"
            or record.get("failed_event_response") != "zero_unsimulated_response"):
        raise ValueError("unsupported shower normalization population")
    maximum = record.get("max_momentum_violation_mev")
    if maximum is not None and _finite(maximum, "max_momentum_violation_mev") < 0:
        raise ValueError("negative maximum momentum violation")
    for key in ("requested_events", "lhe_events", "attempted_events", "generated_events", "saved_events"):
        _integer(record.get(key), key, minimum=1)
    _integer(record.get("discarded_events"), "discarded_events")
    if (record["saved_events"] != record["generated_events"]
            or record["attempted_events"] > record["lhe_events"]
            or record["generated_events"] > record["requested_events"]
            or record["discarded_events"] != record["attempted_events"] - record["generated_events"]):
        raise ValueError("inconsistent source/generated/saved shower counts")
    counts = record.get("exception_counts", {})
    if not isinstance(counts, dict) or set(counts) != set(SEVERITIES):
        raise ValueError("incomplete shower exception bookkeeping")
    for value in counts.values():
        _integer(value, "exception count")
    classes = record.get("exception_classes")
    class_counts = dict.fromkeys(SEVERITIES, 0)
    if not isinstance(classes, list):
        raise ValueError("missing final exception class provenance")
    for entry in classes:
        if (not isinstance(entry, dict) or not isinstance(entry.get("class"), str)
                or entry.get("severity") not in class_counts):
            raise ValueError("invalid final exception class provenance")
        class_counts[entry["severity"]] += _integer(entry.get("count"), "exception class count", minimum=1)
    if class_counts != counts:
        raise ValueError("final exception classes differ from recorded exception counts")
    if (counts["eventerror"] != record["discarded_events"]
            or any(counts[key] for key in ("setuperror", "maybeabort", "abortnow", "unknown"))):
        raise ValueError("unaccounted discarded events or fatal shower exceptions")
    if record["termination"] == "source_exhausted":
        if (record["generated_events"] >= record["requested_events"]
                or record["attempted_events"] != record["lhe_events"] or counts["runerror"] != 1):
            raise ValueError("source exhaustion counts/runerror do not close")
    elif record["generated_events"] != record["requested_events"] or counts["runerror"]:
        raise ValueError("requested-event completion counts/runerror do not close")
    if record.get("sample") != manifest.get("sample") or record["requested_events"] != manifest.get("nevents_requested"):
        raise ValueError("shower completion differs from campaign identity")
    if record["lhe_events"] != manifest.get("lhe", {}).get("events"):
        raise ValueError("shower completion differs from recorded LHE event count")
    roots = _root_summary(record.get("root_metadata"))
    if roots != record.get("root_summary") or roots["entries"] != record["saved_events"]:
        raise ValueError("ROOT metadata differs from completed shower count/weights")
    artifacts = record.get("artifacts", {})
    if not isinstance(artifacts, dict) or any(key not in artifacts for key in ("run_log", "herwig_out", "herwig_log", "roots")):
        raise ValueError("missing shower artifact provenance")
    root_artifacts = artifacts["roots"]
    if (not isinstance(root_artifacts, list) or not all(isinstance(entry, dict) for entry in root_artifacts)
            or len(root_artifacts) != len(record["root_metadata"])
            or [(entry.get("path"), entry.get("size"), entry.get("mtime_ns")) for entry in root_artifacts]
            != [(entry["path"], entry["size"], entry["mtime_ns"]) for entry in record["root_metadata"]]):
        raise ValueError("ROOT artifact provenance differs from inspected metadata")
    inventory = manifest.get("root_files")
    if inventory is not None and inventory != [
            {key: entry[key] for key in ("path", "size", "mtime_ns")}
            for entry in record["root_metadata"]]:
        raise ValueError("campaign ROOT inventory differs from the completed shower")
    for key, artifact in artifacts.items():
        for item in artifact if key == "roots" else [artifact]:
            _validate_artifact(item, verify_file=directory is not None)
    _validate_source(record, roots, require_source=require_source)
    source = record.get("source_weights")
    if source is not None:
        if "herwig_card" not in artifacts:
            raise ValueError("signed source weights require audited Herwig input provenance")
        if record.get("reader_lhe_path") != source["artifact"]["path"]:
            raise ValueError("Herwig reader FileName differs from the signed source-weight proof")
        card_text = manifest.get("configuration", {}).get("herwig_card")
        if card_text is not None and hashlib.sha256(card_text.encode()).hexdigest() != artifacts["herwig_card"]["sha256"]:
            raise ValueError("persisted Herwig input differs from the recorded campaign card")
        native = manifest.get("lhe", {})
        source_path = manifest.get("lhe_file", native.get("source"))
        if source_path is not None and Path(source_path).resolve() != Path(source["artifact"]["path"]).resolve():
            raise ValueError("signed source-weight proof refers to a different LHE")
        if (native.get("source_size", source["artifact"]["size"]) != source["artifact"]["size"]
                or native.get("source_mtime_ns", source["artifact"]["mtime_ns"]) != source["artifact"]["mtime_ns"]):
            raise ValueError("signed source-weight proof differs from native LHE inventory")
        _validate_artifact(source["artifact"], verify_file=directory is not None)
    if directory is not None:
        directory = Path(directory).resolve()
        herwig_dir = Path(record.get("herwig_dir", "")).resolve()
        if herwig_dir != directory / "herwig" and herwig_dir != directory:
            raise ValueError("shower completion belongs to a different sample directory")
        parsed = _parse_logs(herwig_dir, record["sample"], record["termination"])
        for key in ("generated_events", "attempted_events", "saved_events", "exception_counts", "exception_classes",
                    "max_momentum_violation_mev"):
            if parsed[key] != record.get(key):
                raise ValueError("shower completion differs from finalized logs")
        if set(parsed["root_writers"]) != {entry["path"] for entry in root_artifacts}:
            raise ValueError("finalized ROOT writer inventory changed")
        if source is not None:
            if _audit_source_card(herwig_dir / f"{record['sample']}.in") != record["reader_lhe_path"]:
                raise ValueError("Herwig reader FileName changed")
    return record


def expected_analysis_events(manifest):
    record = manifest.get("shower_completion")
    if record is None:
        return _integer(manifest.get("nevents_requested"), "nevents_requested", minimum=1)
    return validate_completion(record, manifest)["saved_events"]


def normalization_denominator(manifest, analysis_sum_weight, weight_scale):
    analysis_sum_weight = _finite(analysis_sum_weight, "analysis_sum_weight")
    weight_scale = _finite(weight_scale, "weight_scale")
    if analysis_sum_weight <= 0 or weight_scale <= 0:
        raise ValueError("normalization requires positive signed weights and weight_scale")
    record = manifest.get("shower_completion")
    if record is None:
        return analysis_sum_weight
    validate_completion(record, manifest)
    if not _close(analysis_sum_weight / weight_scale, record["root_summary"]["sum_weight"]):
        raise ValueError("analysis weights differ from saved ROOT shower weights")
    source = record.get("source_weights")
    denominator = analysis_sum_weight if source is None else source["sum_weight"] * weight_scale
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("non-positive signed source normalization denominator")
    return denominator
