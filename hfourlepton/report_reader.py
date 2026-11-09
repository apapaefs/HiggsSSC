"""Read large detector-hypothesis trees in bounded chunks.

Rejected hypotheses with the same source, channel and cut mask contribute only
additive probabilities and weights to the report. Combine those rows before
creating Python records; keep every selected hypothesis separately for shapes
and ML. Statistical sums of squares are still formed *after* grouping by source.
"""
from __future__ import annotations

from pathlib import Path
import time

import numpy as np


class RejectedHypotheses:
    def __init__(self):
        self.rows = {}

    def add(self, columns, rejected):
        if not np.any(rejected):
            return
        keys = np.rec.fromarrays([columns[name][rejected] for name in
                                  ("source_index", "channel", "cut_mask")])
        unique, inverse = np.unique(keys, return_inverse=True)
        sums = {}
        for name in ("event_weight", "pretrigger_event_weight", "hypothesis_probability"):
            values = columns.get(name, columns["event_weight"])[rejected]
            if name == "hypothesis_probability":
                # Match response_closure_audit: absent/nonfinite probabilities
                # do not contribute, including when every value is absent.
                finite = np.isfinite(values)
                sums["probability_count"] = np.bincount(inverse, weights=finite, minlength=len(unique))
                values = np.where(finite, values, 0.)
            elif not np.all(np.isfinite(values)):
                raise ValueError(f"nonfinite {name} in rejected detector hypotheses")
            sums[name] = np.bincount(inverse, weights=values, minlength=len(unique))
        for i, key in enumerate(unique):
            identity = tuple(int(value) for value in key)
            row = self.rows.setdefault(identity, dict(source_index=identity[0], channel=identity[1],
                cut_mask=identity[2], event_weight=0., pretrigger_event_weight=0.))
            for name in ("event_weight", "pretrigger_event_weight"):
                row[name] += float(sums[name][i])
            if sums["probability_count"][i]:
                row["hypothesis_probability"] = row.get("hypothesis_probability", 0.) + float(
                    sums["hypothesis_probability"][i])


def load_report_records(path: Path, *, selected_mask=1023, step_size="32 MB", progress=False):
    """Keep memory proportional to source events and selected hypotheses.

    ``four_lepton_entries`` records the original row count, including combined
    rejected rows, so detector-provenance checks continue to cover the full tree.
    No ROOT file or manifest is modified.
    """
    try:
        import awkward as ak
        import uproot
    except ImportError as error:
        raise RuntimeError("Chunked four-lepton reporting requires uproot and awkward; "
                           "run run_pipeline.py setup --workflow fourlepton-ho") from error
    from hfourlepton import make_fourlepton_report as report

    sources, selected = [], []
    rejected = RejectedHypotheses()
    audit = report.Audit()
    metadata = {}
    name_lists = {}
    entries = 0
    started = time.monotonic()
    with uproot.open(path, array_cache=None) as root:
        for tree_name, aliases in (("SourceEvents", report.SOURCE_BRANCH_ALIASES),
                                   ("FourLepton", report.FOUR_LEPTON_BRANCH_ALIASES)):
            if tree_name not in root:
                audit.errors.append(f"missing ROOT tree {tree_name}")
                continue
            tree = root[tree_name]
            chosen = report._chosen_branches(tree.keys(), aliases)
            missing = report.REQUIRED_REPORT_BRANCHES[tree_name] - chosen.keys()
            if missing:
                audit.errors.append(f"ROOT tree {tree_name} lacks required report branches: " +
                                    ", ".join(sorted(missing)))
                continue
            if tree_name == "SourceEvents":
                for batch in tree.iterate(list(chosen.values()), step_size=step_size, library="ak"):
                    for values in ak.to_list(batch):
                        record = {canonical: values[actual] for canonical, actual in chosen.items()}
                        names = record.get("optional_weight_names")
                        if names is not None:
                            # Verify all rows later in optional_weight_audit;
                            # sharing identical labels does not hide mismatches.
                            record["optional_weight_names"] = name_lists.setdefault(tuple(names), names)
                        sources.append(record)
                continue
            entries = int(tree.num_entries)
            read = 0
            next_progress = 1_000_000
            for batch in tree.iterate(list(chosen.values()), step_size=step_size, library="np"):
                columns = {canonical: batch[actual] for canonical, actual in chosen.items()}
                mask = (columns["cut_mask"] & selected_mask) == selected_mask
                # The analysis previously checked these kinematics even for
                # rejected rows; combining rows must not conceal invalid data.
                for name in ("m4l", "mZ1", "mZ2", "pt4l", "y4l", "min_delta_r", "lepton_pt"):
                    if name in columns and not np.all(np.isfinite(columns[name])):
                        audit.errors.append(f"nonfinite {name} in FourLepton entries starting at {read}")
                for i in np.flatnonzero(mask):
                    selected.append({name: report._python_value(values[i]) for name, values in columns.items()})
                rejected.add(columns, ~mask)
                read += len(mask)
                if progress and (read >= next_progress or read == entries):
                    print(f"  FourLepton {read:,}/{entries:,}; {len(selected):,} selected hypotheses; "
                          f"{time.monotonic() - started:.1f}s", flush=True)
                    next_progress = read + 1_000_000
            if read != entries:
                raise ValueError(f"incomplete FourLepton read in {path}: {read}/{entries}")
        present = "AnalysisMetadata" in root
        if present:
            tree = root["AnalysisMetadata"]
            chosen = report._chosen_branches(tree.keys(), report.METADATA_BRANCH_ALIASES)
            if chosen and tree.num_entries:
                arrays = tree.arrays(list(chosen.values()), entry_stop=1, library="np")
                metadata.update({name: report._python_value(arrays[actual][0]) for name, actual in chosen.items()})
        else:
            audit.warnings.append("missing optional ROOT tree AnalysisMetadata")
    metadata["analysis_metadata_present"] = present
    metadata["analysis_metadata"] = {key: value for key, value in metadata.items()
                                      if key not in ("analysis_metadata_present", "analysis_metadata")}
    rows = list(rejected.rows.values()) + selected
    records = report.RootRecords(sources, rows, audit, metadata, four_lepton_entries=entries)
    return report._add_summary_sidecar(path, records)
