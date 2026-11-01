#!/usr/bin/env python3
"""Managed minimal runner with a hard-event count gate before any workers.

The compatibility wrapper adds validation without changing the generation
identities of existing campaigns. Validation receipts are separate from physics
stage markers; no producer, cuts, normalization or event budgets are replaced.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hfourlepton import run_four_lepton_campaign as four


def count_events(path):
    """Count closed event blocks and reject a truncated LHE, including gzip."""
    count, inside, closed = 0, False, False
    with four.open_maybe_gzip(path) as stream:
        for line in stream:
            line = line.strip()
            if re.fullmatch(r"<event(?:\s+[^>]*)?>", line):
                if inside or closed:
                    raise ValueError(f"malformed LHE event structure: {path}")
                inside = True
            elif line == "</event>":
                if not inside or closed:
                    raise ValueError(f"malformed LHE event structure: {path}")
                inside = False
                count += 1
            elif line == "</LesHouchesEvents>":
                if inside or closed:
                    raise ValueError(f"malformed LHE event structure: {path}")
                closed = True
    if inside or not closed:
        raise ValueError(f"truncated LHE: {path}")
    return count


def require_events(path, needed):
    count = count_events(path)
    if count < needed:
        raise ValueError(f"{path} has {count} events, expected at least {needed}; "
                         "hard production is incomplete. Preserve the output and inspect the MG5 log.")
    return count


def validate_retry_receipt(marker):
    """An exceptional numerical retry must keep its additional provenance."""
    payload = json.loads(marker.read_text())
    receipt = marker.with_name("numerical_retry.json")
    expected = payload.get("numerical_retry")
    if expected is None and not receipt.exists():
        return
    if expected != four.artifact_stage_record(receipt):
        raise ValueError(f"numerical retry provenance is missing or changed: {receipt}")
    record = json.loads(receipt.read_text())
    if record.get("status") != "complete" or record.get("artifact") != payload["artifact"]:
        raise ValueError(f"numerical retry is not complete: {receipt}")


@contextmanager
def production_checks():
    """Scope the compatibility hook to the sequential production coordinator."""
    producer = four.obtain_production_lhe

    def checked(sample, selected, cfg):
        path = producer(sample, selected, cfg)
        if cfg.preset != "minimal" or cfg.dry_run:
            return path
        # A signal producer returns one decay-channel slice. Check the native
        # group before comparing with the sum of all three channel budgets.
        from hfourlepton import repair_loop_production as production
        record = dict(production.completed_production(sample, selected, cfg),
                      validator=four.artifact_stage_record(Path(__file__).resolve()),
                      inspector=four.artifact_stage_record(Path(production.__file__).resolve()))
        destination = cfg.run_dir / "production_validation" / (sample.production_group_id + ".json")
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(".tmp")
        temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        temporary.replace(destination)
        return path

    four.obtain_production_lhe = checked
    try:
        yield
    finally:
        four.obtain_production_lhe = producer


def main(argv=None):
    with production_checks():
        return four.main(sys.argv[1:] if argv is None else argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(f"ERROR: {error}")
