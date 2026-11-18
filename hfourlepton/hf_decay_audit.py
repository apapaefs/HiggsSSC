"""Independent audit of native Herwig full-support B-decay importance weights.

The unmodified table defines p, the explicitly multiplied direct e/mu modes
define q, and the post-decay handler supplies actual selected mode names. Every
logged p/q is recomputed; finalized ROOT moments must close to the audit. The
pre-decay LHE population always remains the cross-section denominator.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import re

from hgammagamma import ho_shower_completion as completion
from hfourlepton import run_four_lepton_campaign as four

PREFIX = "HIGGSSSC_HF_AUDIT "
PARENTS = {511: "B0", -511: "Bbar0", 521: "B+", -521: "B-", 531: "B_s0", -531: "B_sbar0"}


def mode_key(name):
    name = name.rsplit("/", 1)[-1].rstrip(";")
    parent, children = name.split("->", 1)
    # ThePEG can rewrite colour connections, e.g. d=dbar,dbar=d becomes
    # d=dbar=dbar,d. All these partonic modes have unit proposal bias;
    # the p/q classification depends on their products, not this spelling.
    return parent, tuple(sorted(children.replace("=", ",").split(",")))


def _table(path):
    selected, _ = four.parse_semileptonic_decay_modes(Path(path).read_text())
    commands, _, proposals = four.semileptonic_decay_commands(Path(path))
    keys = {parent: {mode_key(mode) for mode, _ in modes} for parent, modes in selected.items()}
    all_keys = set()
    for line in Path(path).read_text().splitlines():
        match = re.match(r"\s*decaymode\s+(\S+)\s+([\d.eEdD+-]+)\s+1\s+", line)
        if match and float(match[2].replace("D", "e").replace("d", "e")) > 0:
            all_keys.add(mode_key(match[1]))
    return keys, all_keys, proposals, commands


def inspect_log(log, table):
    selected, all_modes, proposals, _ = _table(table)
    values, before_values, events = [], [], set()
    parent_counts = {parent: {"all": 0, "selected": 0} for parent in proposals}
    last_event = -1
    for line in Path(log).read_text().splitlines():
        if PREFIX not in line:
            continue
        row = json.loads(line.split(PREFIX, 1)[1])
        event = completion._integer(row.get("event"), "HF audit event")
        if event in events or event <= last_event:
            raise ValueError("repeated or unordered HF audit event")
        events.add(event)
        last_event = event
        factor = 1.
        for parent in row["parents"]:
            name = PARENTS.get(parent["pdg"])
            if name is None or mode_key(parent["mode"])[0] != name:
                raise ValueError("HF parent/mode identity differs")
            key = mode_key(parent["mode"])
            if key not in all_modes or "MatchAny" in parent["mode"]:
                raise ValueError("HF audit does not use an active native decay mode")
            proposal = proposals[name]
            bias = proposal["bias_factor"] if key in selected[name] else 1.
            expected = proposal["proposal_normalization"] / (proposal["original_active_branching_sum"] * bias)
            for field, desired in (("original_sum", proposal["original_active_branching_sum"]),
                                   ("proposal_sum", proposal["proposal_normalization"]),
                                   ("mode_bias", bias), ("factor", expected)):
                if not completion._close(completion._finite(parent[field], field), desired):
                    raise ValueError(f"HF runtime {field} differs from the independently parsed native table")
            parent_counts[name]["all"] += 1
            parent_counts[name]["selected"] += int(bias != 1.)
            factor *= expected
        before, after = (completion._finite(row[key], key) for key in ("before", "after"))
        if before <= 0 or after <= 0 or not completion._close(row["factor"], factor) or not completion._close(after, before * factor):
            raise ValueError("HF event weight differs from the independently recomputed product of p/q")
        before_values.append(before)
        values.append(after)
    if not values:
        raise ValueError("no native heavy-flavour weight audit events")
    def moments(items):
        return dict(zip(completion.MOMENTS, (math.fsum(items), math.fsum(map(abs, items)), math.fsum(x*x for x in items))))
    return {"events": len(events), "first_event": min(events), "last_event": max(events),
            "root_moments": moments(values), "before_moments": moments(before_values),
            "parent_counts": parent_counts, "parent_proposals": proposals}


def build_audit(log, table, card, plugin_source, plugin):
    record = {"schema_version": 1, "kind": "native_b_full_support_weight_audit",
              "scheme": "full_support_direct_b_semileptonic_importance_sampling",
              "full_tail_support": True,
              "artifacts": {key: completion._artifact(value) for key, value in
                            (("run_log", log), ("table", table), ("card", card),
                             ("plugin_source", plugin_source), ("plugin", plugin))},
              "summary": inspect_log(log, table)}
    _check_card(record)
    record["fingerprint"] = completion.completion_fingerprint(record)
    return record


def _check_card(audit):
    artifacts = audit["artifacts"]
    text = Path(artifacts["card"]["path"]).read_text()
    commands = [line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    table = artifacts["table"]["path"]
    _, _, _, bias_commands = _table(table)
    if commands.count(f"read {table}") != 1:
        raise ValueError("HF card must load exactly the audited native decay table")
    for command in bias_commands:
        if commands.count(command) != 1:
            raise ValueError("HF card differs from the declared full-support proposal")
    if ("HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter" not in text
            or "PostDecayHandlers" not in text
            or f"globallibrary {artifacts['plugin']['path']}" not in commands):
        raise ValueError("HF card lacks the audited late decay reweighter")


def validate_audit(record, source, root_summary, *, verify_file=False):
    audit = record.get("heavy_flavour_weight_audit")
    if (not isinstance(audit, dict) or audit.get("schema_version") != 1
            or audit.get("kind") != "native_b_full_support_weight_audit"
            or audit.get("scheme") != "full_support_direct_b_semileptonic_importance_sampling"
            or audit.get("full_tail_support") is not True
            or record.get("sample") not in {name + "_hfbiased" for name in ("bkg_zbb_ee", "bkg_zbb_mumu", "bkg_ttbar")}
            or audit.get("fingerprint") != completion.completion_fingerprint(audit)):
        raise ValueError("unsupported or changed heavy-flavour weight audit")
    if source.get("negative_events") != 0 or any(not completion._close(source[key], record["attempted_events"]) for key in completion.MOMENTS):
        raise ValueError("HF audit requires positive unit pre-decay LHE weights")
    artifacts = audit["artifacts"]
    if artifacts["card"] != record["artifacts"].get("herwig_card") or artifacts["run_log"] != record["artifacts"].get("run_log"):
        raise ValueError("HF audit belongs to a different shower")
    for artifact in artifacts.values():
        completion._validate_artifact(artifact, verify_file=verify_file)
    summary = audit["summary"]
    if verify_file:
        _check_card(audit)
        if inspect_log(artifacts["run_log"]["path"], artifacts["table"]["path"]) != summary:
            raise ValueError("HF audit summary differs from its actual event records")
    if summary["events"] != record["saved_events"]:
        raise ValueError("HF audit must account for every saved event and no event lost after weighting")
    if any(not completion._close(summary["root_moments"][key], root_summary[key]) for key in completion.MOMENTS):
        raise ValueError("HF recomputed weights do not close to finalized ROOT moments")
    if record["sample"] == "bkg_ttbar_hfbiased":
        if record.get("source_weight_transform") is None:
            raise ValueError("weighted ttbar requires independent prompt W branching proof")
        completion._prompt_weight_factor({**record, "sample": "bkg_ttbar"}, source, summary["before_moments"], verify_file=verify_file)
    else:
        if record.get("source_weight_transform") is not None or any(not completion._close(summary["before_moments"][key], record["saved_events"]) for key in completion.MOMENTS):
            raise ValueError("weighted Zbb must enter late decay weighting with unit source weights")
