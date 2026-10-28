#!/usr/bin/env python3
"""Shared HO mass plots and replotting of saved classifier reports.

Older scores.csv files are joined to ROOT by (sample, tree_entry), checking
sourceevent and signed weights. New reports export all reconstructed variables.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from read_root_varfiles import FEATURE_NAMES


# The signal core is about 0.85--0.9 GeV wide in the SSC response model.
MASS_BIN_WIDTH_GEV = 1.0
MASS_VIEWS = (("overview", np.arange(100., 151., MASS_BIN_WIDTH_GEV)),
              ("peak", np.arange(120., 131., MASS_BIN_WIDTH_GEV)))


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _keys(*arrays):
    return np.rec.fromarrays(arrays)


def _root_variables(path, provenance):
    """Read the complete eligible population and check full-tree closure."""
    from read_root_varfiles import _open_tree
    import ROOT

    before = path.stat()
    if before.st_size != provenance["input_size"]:
        raise ValueError(f"ROOT size differs from the saved classifier input: {path}")
    _, handle, tree = _open_tree(path, strict=True)
    try:
        leaf = tree.GetLeaf("sourceevent")
        if (not leaf or leaf.GetLenStatic() != 1 or leaf.GetLeafCount()
                or leaf.GetTypeName() != "Long64_t"):
            raise ValueError(f"Invalid sourceevent branch: {path}")
    finally:
        handle.Close()
    frame = ROOT.RDataFrame("Data2", str(path))
    # ROOT versions differ in whether a fixed one-element leaf is exposed as
    # a scalar or RVec by RDataFrame (the persisted branch is unchanged).
    weight = "eventweight[0]" if "RVec" in frame.GetColumnType("eventweight") else "eventweight"
    source = "sourceevent[0]" if "RVec" in frame.GetColumnType("sourceevent") else "sourceevent"
    frame = frame.Define("raw_signed_weight", weight).Define("absolute_weight", "std::abs(raw_signed_weight)")
    count, total, absolute = frame.Count(), frame.Sum("raw_signed_weight"), frame.Sum("absolute_weight")
    selected = (frame.Filter("variables[9] >= 2")
                .Define("tree_entry", "static_cast<Long64_t>(rdfentry_)")
                .Define("source_id", source))
    for index, name in enumerate(FEATURE_NAMES):
        selected = selected.Define(name, f"variables[{index}]")
    arrays = selected.AsNumpy([*FEATURE_NAMES, "tree_entry", "source_id", "raw_signed_weight"])
    if count.GetValue() != provenance["tree_entries"]:
        raise ValueError(f"ROOT entry count differs from the saved report: {path}")
    if not math.isclose(total.GetValue(), provenance["sum_tree_weight"], rel_tol=0,
                        abs_tol=1e-8 * max(absolute.GetValue(), 1e-300)):
        raise ValueError(f"ROOT signed weight sum differs from the saved report: {path}")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"ROOT input changed while reading: {path}")
    return arrays


def load_events(report_dir, input_dir=None):
    metadata_path, scores_path = report_dir / "model_metadata.json", report_dir / "scores.csv"
    metadata = json.loads(metadata_path.read_text())
    with scores_path.open(newline="") as stream:
        records = list(csv.DictReader(stream))
    if not records:
        raise ValueError("No eligible events in scores.csv")
    data = {name: np.asarray([row[name] for row in records]) for name in ("sample", "partition")}
    for name in ("sourceevent", "tree_entry", "label"):
        data[name] = np.asarray([int(row[name]) for row in records], dtype=np.int64)
    for name in ("score", "raw_signed_weight", "physical_weight", "partition_fraction", "projected_physical_weight"):
        data[name] = np.asarray([float(row[name]) for row in records])
        if not np.all(np.isfinite(data[name])):
            raise ValueError(f"Nonfinite {name} in scores.csv")
    if (not np.all(np.isin(data["label"], [0, 1])) or np.any(data["sourceevent"] < 0)
            or np.any(data["tree_entry"] < 0) or np.any((data["score"] < 0) | (data["score"] > 1))):
        raise ValueError("Invalid labels, identities, or scores")
    if len(np.unique(_keys(data["sample"], data["tree_entry"]))) != len(records):
        raise ValueError("Duplicate (sample, tree_entry) in scores.csv")
    groups = _keys(data["sample"], data["sourceevent"])
    if len(np.unique(groups)) != len(np.unique(_keys(data["sample"], data["sourceevent"], data["partition"]))):
        raise ValueError("A source event appears in multiple partitions")
    for partition in ("train", "validation", "test"):
        if np.sum(data["partition"] == partition) != metadata["eligible_rows_per_partition"][partition]:
            raise ValueError(f"Incomplete {partition} scores")
    if not np.all(np.isin(data["partition"], ["train", "validation", "test"])):
        raise ValueError("Unknown score partition")

    embedded = all(name in records[0] for name in FEATURE_NAMES)
    for name in FEATURE_NAMES:
        data[name] = (np.asarray([float(row[name]) for row in records]) if embedded
                      else np.full(len(records), np.nan))
    provenance = {item["sample"]: item for item in metadata["normalization_provenance"]}
    if set(data["sample"]) - set(provenance):
        raise ValueError("Scores contain samples absent from normalization provenance")
    sources = {"metadata_sha256": _sha256(metadata_path), "scores_sha256": _sha256(scores_path),
               "variables_source": "scores.csv" if embedded else "ROOT join", "root_files": {}}
    for sample, original in provenance.items():
        mask = data["sample"] == sample
        factor = (1000 * metadata["luminosity_fb"] * original["cross_section_pb"]
                  * original["weight_scale"] / original["normalization_sum_weight"])
        if not math.isfinite(factor) or factor <= 0:
            raise ValueError(f"Invalid physical normalization: {sample}")
        fractions = np.asarray([metadata["sample_split_fractions"][sample][p] for p in data["partition"][mask]])
        if (np.any(fractions <= 0) or not np.allclose(data["partition_fraction"][mask], fractions, rtol=1e-12, atol=0)
                or not np.allclose(data["physical_weight"][mask], data["raw_signed_weight"][mask] * factor, rtol=1e-10, atol=1e-14)
                or not np.allclose(data["projected_physical_weight"][mask], data["physical_weight"][mask] / fractions, rtol=1e-10, atol=1e-14)):
            raise ValueError(f"Saved physical weights or partition fractions do not close: {sample}")
        if embedded:
            continue
        path = input_dir / f"{sample}.root" if input_dir else Path(original["input_file"])
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}; use --input-dir containing SAMPLE.root copies of the original inputs")
        print(f"Joining saved scores to {path}", flush=True)
        arrays = _root_variables(path, original)
        indices = np.flatnonzero(mask)
        indices = indices[np.argsort(data["tree_entry"][indices])]
        order = np.argsort(arrays["tree_entry"])
        if (not np.array_equal(data["tree_entry"][indices], arrays["tree_entry"][order])
                or not np.array_equal(data["sourceevent"][indices], arrays["source_id"][order])):
            raise ValueError(f"ROOT/score event identities or eligible populations differ: {sample}")
        if not np.allclose(data["raw_signed_weight"][indices], arrays["raw_signed_weight"][order], rtol=1e-12, atol=1e-18):
            raise ValueError(f"ROOT/score event weights differ: {sample}")
        for name in FEATURE_NAMES:
            data[name][indices] = arrays[name][order]
        sources["root_files"][sample] = {"path": str(path.resolve()), "sha256": _sha256(path)}
    for name in FEATURE_NAMES:
        if not np.all(np.isfinite(data[name])) or np.any(data[name] <= -900):
            raise ValueError(f"Invalid reconstructed variable: {name}")
    if np.any(data["n_selected_photons"] < 2):
        raise ValueError("scores.csv must contain eligible two-photon rows only")
    return metadata, data, sources


def selection_masks(data, cuts, threshold):
    eligible = data["n_selected_photons"] >= 2
    baseline, other = eligible.copy(), eligible.copy()
    for specification in cuts:
        variable = specification["variable"]
        if variable not in FEATURE_NAMES:
            raise ValueError(f"Unknown cut variable: {variable}")
        minimum, maximum = specification.get("min"), specification.get("max")
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValueError(f"Cut minimum exceeds maximum: {variable}")
        values = data[variable]
        accepted = np.isfinite(values) & (values > -900)
        if minimum is not None:
            accepted &= values >= minimum
        if maximum is not None:
            accepted &= values <= maximum
        baseline &= accepted
        if variable != "m_gg":
            other &= accepted
    masks = {"two_photons": eligible, "other_cuts": other, "all_cuts": baseline}
    if threshold is not None:
        if not math.isfinite(threshold) or not 0 <= threshold <= 1:
            raise ValueError("Invalid saved threshold")
        masks["xgboost"] = eligible & (data["score"] >= threshold)
    return masks


def histogram(mass, weights, samples, sources, edges):
    """Signed counts; combine each source's contributions in a bin before squaring."""
    bins = np.searchsorted(edges, mass, side="right") - 1
    bins[mass == edges[-1]] = len(edges) - 2
    valid = (bins >= 0) & (bins < len(edges) - 1)
    keys = _keys(samples[valid], sources[valid], bins[valid])
    unique, inverse = np.unique(keys, return_inverse=True)
    grouped = np.bincount(inverse, weights=weights[valid])
    counts = np.bincount(bins[valid], weights=weights[valid], minlength=len(edges) - 1)
    variance = np.bincount(unique[unique.dtype.names[2]].astype(int), weights=grouped**2, minlength=len(edges) - 1)
    return counts, np.sqrt(variance)


def mass_summary(mass, weights):
    """Signed moments and explicitly absolute-weight central quantiles, without cuts."""
    total, absolute = float(np.sum(weights)), float(np.sum(np.abs(weights)))
    output = {"entries": len(mass), "signed_events": total, "mean_gev": None, "rms_gev": None,
              "q16_gev": None, "median_gev": None, "q84_gev": None, "sigma68_gev": None,
              "quantile_weights": "absolute physical weights", "negative_absolute_weight_fraction": None,
              "minimum_gev": None, "maximum_gev": None}
    if absolute <= 0:
        return output
    output.update(minimum_gev=float(np.min(mass)), maximum_gev=float(np.max(mass)),
                  negative_absolute_weight_fraction=float(np.sum(np.abs(weights[weights < 0])) / absolute))
    order = np.argsort(mass)
    ordered_weights = np.abs(weights[order])
    nonzero = ordered_weights > 0
    cumulative = (np.cumsum(ordered_weights[nonzero]) - .5 * ordered_weights[nonzero]) / absolute
    q16, median, q84 = np.interp([.16, .5, .84], cumulative, mass[order][nonzero])
    output.update(q16_gev=float(q16), median_gev=float(median), q84_gev=float(q84), sigma68_gev=float((q84-q16)/2))
    if total > 0:
        mean = float(np.sum(mass * weights) / total)
        variance = float(np.sum((mass - mean)**2 * weights) / total)
        output.update(mean_gev=mean, rms_gev=math.sqrt(variance) if variance >= 0 else None)
    return output


def _spectrum_plot(output, stem, data, selections, weights, scope, title, luminosity, cuts, plt):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    table = []
    for axes_row, (view, edges) in zip(axes, MASS_VIEWS):
        for ax, label, category in zip(axes_row, [1, 0], ["Signal", "Total background"]):
            for key, name, color, selection in selections:
                mask = scope & selection & (data["label"] == label)
                values, errors = histogram(data["m_gg"][mask], weights[mask], data["sample"][mask], data["sourceevent"][mask], edges)
                centers = (edges[1:] + edges[:-1]) / 2
                ax.stairs(values, edges, label=name, color=color, linewidth=1.6)
                ax.fill_between(centers, values-errors, values+errors, step="mid", color=color, alpha=.16)
                for low, high, value, error in zip(edges[:-1], edges[1:], values, errors):
                    table.append(dict(view=view, category=category, selection=key, low_gev=low, high_gev=high,
                                      expected_events=value, mc_error_events=error))
            for cut in cuts:
                if cut["variable"] == "m_gg":
                    for boundary in (cut.get("min"), cut.get("max")):
                        if boundary is not None:
                            ax.axvline(boundary, color="0.5", linestyle=":", linewidth=1)
            ax.axhline(0, color="0.3", linewidth=.6)
            ax.set(xlim=(edges[0], edges[-1]), xlabel=r"$m_{\gamma\gamma}$ [GeV]",
                   ylabel=f"Expected events / {edges[1]-edges[0]:g} GeV", title=f"{category} · {view}")
            ax.legend(fontsize=9)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
    fig.suptitle(f"{title}\nSSC 40 TeV · {luminosity:g} fb$^{{-1}}$ · bands: source-event MC uncertainty", fontsize=13)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"{stem}.{extension}", dpi=180)
    plt.close(fig)
    with (output / f"{stem}.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)


def _pyplot():
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def sample_plot_data(sample, rows, weights, source_events, luminosity):
    """Build compact two-candidate plotting arrays from a validated HO sample."""
    if not (len(rows) == len(weights) == len(source_events)):
        raise ValueError(f"Inconsistent mass-plot population: {sample.name}")
    indices = [i for i, row in enumerate(rows) if row["n_selected_photons"] >= 2]
    data = {name: np.asarray([rows[i][name] for i in indices], dtype=float) for name in FEATURE_NAMES}
    for name, values in data.items():
        if not np.all(np.isfinite(values)) or np.any(values <= -900):
            raise ValueError(f"Invalid selected {name}: {sample.name}")
    factor = 1000 * luminosity * sample.cross_section_pb * sample.weight_scale / sample.normalization_sum_weight
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError(f"Invalid physical mass-plot normalization: {sample.name}")
    data.update(sample=np.full(len(indices), sample.name),
                label=np.full(len(indices), int(sample.category == "Signal")),
                sourceevent=np.asarray([source_events[i] for i in indices], dtype=np.int64),
                physical_weight=np.asarray([weights[i] * factor for i in indices]))
    return data


def write_cut_mass_plots(output, data, cuts, luminosity, plt=None):
    """Write the full-sample mass spectrum after all cuts except mass, plus the baseline."""
    plt = plt if plt is not None else _pyplot()
    output.mkdir(parents=True, exist_ok=True)
    masks = selection_masks(data, cuts, None)
    _spectrum_plot(output, "m_gg_cut_based", data,
                   [("other_cuts", "All cuts except mass", "#236192", masks["other_cuts"]),
                    ("all_cuts", "All cuts including mass", "#d47500", masks["all_cuts"])],
                   data["physical_weight"], np.ones(len(data["label"]), dtype=bool),
                   "Cut-based selection · full sample", luminosity, cuts, plt)
    return [output / f"m_gg_cut_based.{extension}" for extension in ("png", "pdf", "csv")]


def write_xgboost_mass_plots(output, data, threshold, luminosity, plt=None):
    """Use only independent test events and the frozen validation threshold."""
    plt = plt if plt is not None else _pyplot()
    output.mkdir(parents=True, exist_ok=True)
    masks = selection_masks(data, [], threshold)
    test = data["partition"] == "test"
    projected = data["projected_physical_weight"]
    selections = [("two_photons", "Before score selection", "#236192", masks["two_photons"])]
    title = "XGBoost selection · independent test sample, projected to full luminosity"
    if threshold is not None:
        selections.append(("xgboost", f"Score ≥ {threshold:g}", "#20936e", masks["xgboost"]))
    else:
        title = "Independent test sample · no supported validation threshold; no score selection"
    _spectrum_plot(output, "m_gg_xgboost", data, selections, projected, test, title, luminosity, [], plt)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    edges = MASS_VIEWS[1][1]
    table = []
    for ax, label, category in zip(axes, [1, 0], ["Signal", "Background"]):
        before = test & masks["two_photons"] & (data["label"] == label)
        denominator = np.histogram(data["m_gg"][before], edges, weights=projected[before])[0]
        numerator = np.full(len(edges)-1, np.nan)
        ratio = np.full(len(edges)-1, np.nan)
        if threshold is not None:
            after = before & masks["xgboost"]
            numerator = np.histogram(data["m_gg"][after], edges, weights=projected[after])[0]
            ratio = np.divide(numerator, denominator, out=ratio, where=denominator > 0)
            ax.stairs(ratio, edges, color="#20936e", linewidth=1.8)
        else:
            ax.text(.5, .5, "No supported validation threshold", ha="center", transform=ax.transAxes)
        ax.axhline(0, color="0.5", linewidth=.5)
        ax.set(xlabel=r"$m_{\gamma\gamma}$ [GeV]", ylabel="Signed yield after / before score cut",
               title=category, xlim=(edges[0], edges[-1]))
        for low, high, total, selected, efficiency in zip(edges[:-1], edges[1:], denominator, numerator, ratio):
            table.append(dict(category=category, low_gev=low, high_gev=high, before_events=total,
                              after_events=float(selected) if np.isfinite(selected) else None,
                              efficiency=float(efficiency) if np.isfinite(efficiency) else None))
    title = f"Mass dependence of score ≥ {threshold:g}" if threshold is not None else "Mass acceptance unavailable"
    fig.suptitle(f"{title} · independent test events\n{MASS_BIN_WIDTH_GEV:g} GeV bins; bins with nonpositive signed denominator are omitted", fontsize=12)
    for extension in ("png", "pdf"):
        fig.savefig(output / f"mass_acceptance.{extension}", dpi=180)
    plt.close(fig)
    with (output / "mass_acceptance.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    return [output / f"{stem}.{extension}" for stem in ("m_gg_xgboost", "mass_acceptance")
            for extension in ("png", "pdf", "csv")]


def make_report(report_dir, output, input_dir=None):
    metadata, data, sources = load_events(report_dir, input_dir)
    threshold = metadata["best_threshold"]
    cuts = metadata["baseline_cuts"]
    masks = selection_masks(data, cuts, threshold)
    test = data["partition"] == "test"
    weights, projected = data["physical_weight"], data["projected_physical_weight"]
    summary = {"source_report": str(report_dir.resolve()), "run_tag": metadata["run_tag"],
               "luminosity_fb": metadata["luminosity_fb"], "threshold": threshold,
               "mass_is_classifier_input": "m_gg" in metadata["feature_names"], "baseline_cuts": cuts,
               "provenance": sources, "selections": {},
               "resolution_definition": "Preselection signal sigma68 = (q84-q16)/2 using absolute physical weights; no mass window. Post-selection widths are conditional shapes, not detector resolution.",
               "mc_error_definition": "Square the sum of each source event's signed contributions within each bin; fixed normalization, no systematic uncertainties."}
    for scope_name, scope, scope_weights in (("full_sample", np.ones(len(test), dtype=bool), weights), ("test", test, projected)):
        results = {}
        for key, mask in masks.items():
            # Model evaluation is restricted to events that were not used to fit or tune it.
            if key == "xgboost" and scope_name != "test":
                continue
            results[key] = {category: mass_summary(data["m_gg"][scope & mask & (data["label"] == label)],
                                                  scope_weights[scope & mask & (data["label"] == label)])
                            for label, category in ((1, "signal"), (0, "background"))}
        summary["selections"][scope_name] = results
    signal_before = masks["other_cuts"] & (data["label"] == 1)
    core_mask = signal_before & (data["m_gg"] >= 120) & (data["m_gg"] <= 130)
    core = mass_summary(data["m_gg"][core_mask], weights[core_mask])
    signal_total = summary["selections"]["full_sample"]["other_cuts"]["signal"]["signed_events"]
    summary["signal_core_diagnostic"] = {"window_gev": [120, 130], **core,
                                         "signed_yield_fraction_outside": 1 - core["signed_events"] / signal_total if signal_total > 0 else None}
    if threshold is not None:
        summary["test_score_efficiencies"] = {}
        for category in ("signal", "background"):
            before = summary["selections"]["test"]["two_photons"][category]["signed_events"]
            after = summary["selections"]["test"]["xgboost"][category]["signed_events"]
            summary["test_score_efficiencies"][category] = after / before if before > 0 else None
    # Check the diagnostic reproduces the saved independent test yields.
    comparison_path = report_dir / "comparison.json"
    if comparison_path.exists():
        comparison = json.loads(comparison_path.read_text())
        for selection, key in (("baseline", "all_cuts"), ("xgboost", "xgboost")):
            if key not in summary["selections"]["test"]:
                continue
            for category in ("signal", "background"):
                actual = summary["selections"]["test"][key][category]["signed_events"]
                expected = comparison[selection]["totals"][f"{category}_expected_events"]
                if not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-7):
                    raise ValueError(f"Saved {selection} {category} yield does not match the score reconstruction")
        summary["saved_test_yields_verified"] = True
        summary["saved_test_comparison"] = {selection: comparison[selection]["totals"]
                                             for selection in ("baseline", "xgboost") if comparison[selection] is not None}

    plt = _pyplot()
    write_cut_mass_plots(output, data, cuts, metadata["luminosity_fb"], plt)
    write_xgboost_mass_plots(output, data, threshold, metadata["luminosity_fb"], plt)
    plots = ["m_gg_cut_based", "m_gg_xgboost", "mass_acceptance"]
    summary["mass_bin_width_gev"] = MASS_BIN_WIDTH_GEV

    (output / "mass_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    resolution = summary["selections"]["full_sample"]["other_cuts"]["signal"]
    def fmt(value):
        return "undefined" if value is None else f"{value:.4g}"
    pictures = "".join(f'<figure><img src="{name}.png" alt="{name}"><figcaption><a href="{name}.pdf">PDF</a>'
                       + f' · <a href="{name}.csv">CSV</a>'
                       + '</figcaption></figure>' for name in plots)
    cut_text = "; ".join(f"{c['variable']}: {c.get('min')} to {c.get('max')}" for c in cuts)
    (output / "index.html").write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8">
<title>HO diphoton mass diagnostics</title><style>body{{font:17px system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;line-height:1.5;color:#243240}}img{{max-width:100%}}figure{{margin:2rem 0}}a{{color:#236192}}</style>
<h1>HO diphoton mass diagnostics</h1>
<p>Source: {html.escape(str(report_dir))}. Campaign tag: {html.escape(metadata['run_tag'])}.
Threshold: {fmt(threshold)}. <a href="mass_summary.json">Measurements and provenance</a>.</p>
<p>Baseline: {html.escape(cut_text)}. The cut-based plot removes every explicit m_gg cut in its blue curve;
all detector selections and other analysis cuts still apply. The orange curve includes the mass window.</p>
<p>The classifier uses m_gg among its inputs and has no mass-window preselection. Its plot uses only the
independent test partition, with the saved threshold and each sample's recorded projection factor.</p>
<p>Signal before the mass cut: mean {fmt(resolution['mean_gev'])} GeV,
central-68% half-width {fmt(resolution['sigma68_gev'])} GeV.
Central quantiles use absolute physical weights; the negative fraction of absolute signal weight is
{fmt(resolution['negative_absolute_weight_fraction'])}. The signed RMS inside 120–130 GeV is
{fmt(core['rms_gev'])} GeV. The unrestricted signed RMS is {fmt(resolution['rms_gev'])} GeV;
rare far-tail events make it a poor description of the peak core. The signed yield fraction outside
120–130 GeV is {fmt(summary['signal_core_diagnostic']['signed_yield_fraction_outside'])}.
Widths after a mass-dependent score cut describe
the selected shape and cannot be interpreted as improved detector resolution.</p>
<p>Mass plots use {MASS_BIN_WIDTH_GEV:g} GeV bins. Physical histograms use signed weights. Bands include finite Monte Carlo statistics with response
hypotheses grouped by source event. The fixed normalization and detector/systematic uncertainties are
not varied. Viewing ranges do not impose selection cuts. No threshold has been retuned.</p>
{pictures}</html>''')
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xgboost-dir", type=Path, required=True, help="Existing HO classifier report")
    parser.add_argument("--input-dir", type=Path, help="Optional original ROOT copies named SAMPLE.root")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = make_report(args.xgboost_dir, args.output_dir, args.input_dir)
    print(json.dumps(summary["selections"], indent=2, allow_nan=False))
    print(f"Report: {(args.output_dir / 'index.html').resolve()}")


if __name__ == "__main__":
    main()
