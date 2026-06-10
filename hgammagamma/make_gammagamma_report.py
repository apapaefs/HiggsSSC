#!/usr/bin/env python3
"""Build an HTML report from LO h -> gamma gamma analysis .top files.

The report stacks backgrounds first and the Higgs signal on top.  By default
each histogram shape is normalized so its area is sigma * weight_scale *
analysis efficiency for that sample.  The analysis efficiency and weight scale
are read from the matching .dat file.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ANALYSIS_ROOT = SCRIPT_DIR / "LOAnalysis"


@dataclass
class Histogram:
    title: str
    x_min: float
    x_max: float
    x: list[float]
    y: list[float]


@dataclass
class SampleResult:
    name: str
    label: str
    category: str
    sample_dir: Path
    top_file: Path
    dat_file: Path
    cross_section_pb: float
    cross_section_error_pb: float | None
    events_read: float
    selected_events: float
    sum_weight: float
    sum_diphoton_weight: float
    weight_scale: float
    efficiency: float
    selected_cross_section_pb: float
    histograms: dict[str, Histogram]


@dataclass(frozen=True)
class PlotMeta:
    slug: str
    display_name: str
    x_label: str
    unit: str
    x_range: tuple[float, float] | None = None


PLOT_META: dict[str, PlotMeta] = {
    "number of selected photons": PlotMeta("n_selected_photons", "Selected Photon Multiplicity", r"$N_\gamma$", ""),
    "pT of selected photons": PlotMeta("pt_photons", "Selected Photon Transverse Momentum", r"$p_T^\gamma$ [GeV]", "GeV", (0.0, 250.0)),
    "eta of selected photons": PlotMeta("eta_photons", "Selected Photon Pseudorapidity", r"$\eta_\gamma$", ""),
    "pT of leading photon": PlotMeta("pt_gamma1", "Leading Photon Transverse Momentum", r"$p_T^{\gamma_1}$ [GeV]", "GeV", (0.0, 250.0)),
    "eta of leading photon": PlotMeta("eta_gamma1", "Leading Photon Pseudorapidity", r"$\eta_{\gamma_1}$", ""),
    "pT of subleading photon": PlotMeta("pt_gamma2", "Subleading Photon Transverse Momentum", r"$p_T^{\gamma_2}$ [GeV]", "GeV", (0.0, 250.0)),
    "eta of subleading photon": PlotMeta("eta_gamma2", "Subleading Photon Pseudorapidity", r"$\eta_{\gamma_2}$", ""),
    "diphoton invariant mass": PlotMeta("m_gg", "Diphoton Invariant Mass", r"$m_{\gamma\gamma}$ [GeV]", "GeV", (100.0, 160.0)),
    "DeltaR of two leading photons": PlotMeta("deltaR_gg", "Diphoton Separation", r"$\Delta R_{\gamma\gamma}$", ""),
    "DeltaPhi of two leading photons": PlotMeta("deltaPhi_gg", "Diphoton Azimuthal Separation", r"$\Delta\phi_{\gamma\gamma}$ [rad]", "rad"),
    "pT of diphoton system": PlotMeta("pt_gg", "Diphoton Transverse Momentum", r"$p_T^{\gamma\gamma}$ [GeV]", "GeV", (0.0, 250.0)),
    "rapidity of diphoton system": PlotMeta("y_gg", "Diphoton Rapidity", r"$y_{\gamma\gamma}$", ""),
}


SAMPLE_LABELS = {
    "signal_gg_h_aa": "Higgs signal",
    "bkg_prompt_aa": "Prompt gamma gamma",
    "bkg_loop_gg_aa": "Loop-induced gg -> gamma gamma",
    "bkg_gamma_j": "gamma + j fake",
    "bkg_jj": "jj fake",
    "bkg_ee": "e+e- fake",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--run-tag", default=os.environ.get("RUN_TAG", "run_01"))
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--samples", default="all", help="Comma-separated sample names, or 'all'.")
    parser.add_argument(
        "--normalization",
        choices=("selected_xsec", "event_xsec", "unit_area"),
        default="selected_xsec",
        help="selected_xsec gives each plotted histogram area sigma*weight_scale*efficiency.",
    )
    parser.add_argument(
        "--no-density",
        action="store_true",
        help="Do not divide by bin width.  The default gives differential pb/unit where possible.",
    )
    parser.add_argument("--title", default="LO h -> gamma gamma analysis")
    parser.add_argument("--allow-missing-xsec", action="store_true")
    return parser


def log(message: str) -> None:
    print(f"==> {message}")


def die(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def slugify(text: str) -> str:
    text = text.lower().replace("gamma", "g")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_") or "plot"


def parse_top_file(path: Path) -> dict[str, Histogram]:
    histograms: dict[str, Histogram] = {}
    current_lines: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("NEW FRAME"):
            if current_lines:
                add_top_block(current_lines, histograms)
            current_lines = [line]
        else:
            current_lines.append(line)
    if current_lines:
        add_top_block(current_lines, histograms)
    return histograms


def add_top_block(lines: Sequence[str], histograms: dict[str, Histogram]) -> None:
    title = None
    x_min = None
    x_max = None
    x_values: list[float] = []
    y_values: list[float] = []

    data_started = False
    for line in lines:
        title_match = re.match(r'TITLE TOP\s+"(.*)"', line)
        if title_match:
            title = title_match.group(1)
            continue
        limits_match = re.match(r"SET LIMITS X\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", line)
        if limits_match:
            x_min = float(limits_match.group(1))
            x_max = float(limits_match.group(2))
            data_started = True
            continue
        if line.strip().startswith("HIST"):
            data_started = False
            continue
        if not data_started:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0])
            y = float(parts[1])
        except ValueError:
            continue
        if math.isnan(y):
            y = 0.0
        x_values.append(x)
        y_values.append(y)

    if not title or title == "dummy histo" or x_min is None or x_max is None or not x_values:
        return
    histograms[title] = Histogram(title=title, x_min=x_min, x_max=x_max, x=x_values, y=y_values)


def parse_key_value_dat(path: Path) -> dict[str, float | str]:
    values: dict[str, float | str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            continue
        key, value = parts
        try:
            values[key] = float(value)
        except ValueError:
            values[key] = value
    return values


def parse_cross_section(sample_dir: Path, run_tag: str) -> tuple[float, float | None]:
    banner = sample_dir / "mg5_process" / "Events" / run_tag / f"{run_tag}_tag_1_banner.txt"
    if banner.exists():
        text = banner.read_text(errors="ignore")
        match = re.search(r"Integrated weight \(pb\)\s*:\s*([-+0-9.eE]+)", text)
        if match:
            return float(match.group(1)), None

    html_path = sample_dir / "mg5_process" / "HTML" / run_tag / "results.html"
    if html_path.exists():
        text = html_path.read_text(errors="ignore")
        match = re.search(r"s=\s*([-+0-9.eE]+)\s*&#177\s*([-+0-9.eE]+)\s*\(pb\)", text)
        if match:
            return float(match.group(1)), float(match.group(2))

    crossx = sample_dir / "mg5_process" / "crossx.html"
    if crossx.exists():
        text = crossx.read_text(errors="ignore")
        match = re.search(r"results\.html\">\s*([-+0-9.eE]+)\s*<font[^>]*>.*?</font>\s*([-+0-9.eE]+)", text)
        if match:
            return float(match.group(1)), float(match.group(2))

    raise FileNotFoundError(f"could not find MG5 cross section for {sample_dir}")


def discover_sample_dirs(analysis_root: Path) -> list[tuple[str, Path]]:
    discovered: list[tuple[str, Path]] = []
    for category in ("Backgrounds", "Signal"):
        events_dir = analysis_root / category / "events"
        if not events_dir.exists():
            continue
        for sample_dir in sorted(path for path in events_dir.iterdir() if path.is_dir()):
            discovered.append((category, sample_dir))
    return discovered


def requested_sample(name: str, requested: str) -> bool:
    if requested == "all":
        return True
    names = {part.strip() for part in requested.split(",") if part.strip()}
    return name in names


def find_analysis_files(sample_dir: Path, run_tag: str) -> tuple[Path, Path] | None:
    top_files = sorted(sample_dir.glob(f"*-{run_tag}.top"))
    if not top_files:
        return None
    top_file = top_files[-1]
    dat_file = top_file.with_suffix(".dat")
    if not dat_file.exists():
        return None
    return top_file, dat_file


def load_samples(args: argparse.Namespace) -> list[SampleResult]:
    samples: list[SampleResult] = []
    for category, sample_dir in discover_sample_dirs(args.analysis_root):
        name = sample_dir.name
        if not requested_sample(name, args.samples):
            continue
        files = find_analysis_files(sample_dir, args.run_tag)
        if files is None:
            continue
        top_file, dat_file = files
        dat = parse_key_value_dat(dat_file)
        try:
            cross_section_pb, cross_section_error_pb = parse_cross_section(sample_dir, args.run_tag)
        except FileNotFoundError as error:
            if not args.allow_missing_xsec:
                raise
            print(f"Warning: {error}; using 1 pb")
            cross_section_pb, cross_section_error_pb = 1.0, None

        events_read = float(dat.get("events_read", 0.0))
        selected_events = float(dat.get("events_with_two_selected_photons", 0.0))
        sum_weight = float(dat.get("sum_weight", events_read))
        sum_diphoton_weight = float(dat.get("sum_diphoton_weight", selected_events))
        weight_scale = float(dat.get("weight_scale", 1.0))
        efficiency = sum_diphoton_weight / sum_weight if sum_weight > 0 else 0.0
        selected_xsec = cross_section_pb * weight_scale * efficiency

        samples.append(
            SampleResult(
                name=name,
                label=SAMPLE_LABELS.get(name, name.replace("_", " ")),
                category=category,
                sample_dir=sample_dir,
                top_file=top_file,
                dat_file=dat_file,
                cross_section_pb=cross_section_pb,
                cross_section_error_pb=cross_section_error_pb,
                events_read=events_read,
                selected_events=selected_events,
                sum_weight=sum_weight,
                sum_diphoton_weight=sum_diphoton_weight,
                weight_scale=weight_scale,
                efficiency=efficiency,
                selected_cross_section_pb=selected_xsec,
                histograms=parse_top_file(top_file),
            )
        )
    return samples


def ordered_samples(samples: Sequence[SampleResult]) -> list[SampleResult]:
    backgrounds = sorted((sample for sample in samples if sample.category != "Signal"), key=lambda s: s.name)
    signals = sorted((sample for sample in samples if sample.category == "Signal"), key=lambda s: s.name)
    return backgrounds + signals


def bin_width(histogram: Histogram) -> float:
    if not histogram.x:
        return 1.0
    width = (histogram.x_max - histogram.x_min) / len(histogram.x)
    return width if width > 0 else 1.0


def bin_edges(histogram: Histogram) -> list[float]:
    width = bin_width(histogram)
    if not histogram.x:
        return [histogram.x_min, histogram.x_max]
    return [histogram.x[0] - 0.5 * width, *[x + 0.5 * width for x in histogram.x]]


def step_values(values: Sequence[float]) -> list[float]:
    if not values:
        return []
    return [*values, values[-1]]


def scaled_histogram(histogram: Histogram, sample: SampleResult, normalization: str, density: bool) -> list[float]:
    total = sum(histogram.y)
    if total <= 0:
        return [0.0 for _ in histogram.y]

    if normalization == "selected_xsec":
        target = sample.selected_cross_section_pb
        values = [value / total * target for value in histogram.y]
    elif normalization == "event_xsec":
        denominator = sample.sum_weight if sample.sum_weight > 0 else total
        values = [value / denominator * sample.cross_section_pb * sample.weight_scale for value in histogram.y]
    else:
        values = [value / total for value in histogram.y]

    if density:
        width = bin_width(histogram)
        values = [value / width for value in values]
    return values


def common_histogram_titles(samples: Sequence[SampleResult]) -> list[str]:
    titles: set[str] = set()
    for sample in samples:
        titles.update(sample.histograms.keys())
    return [title for title in PLOT_META if title in titles] + sorted(titles - set(PLOT_META))


def ensure_matplotlib(output_dir: Path):
    mpl_dir = output_dir / ".mplconfig"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "black",
            "axes.linewidth": 1.15,
            "axes.grid": False,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "legend.fontsize": 10,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "xtick.top": True,
            "ytick.right": True,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )
    import matplotlib.pyplot as plt

    return plt


def publication_color(sample: SampleResult, background_index: int) -> str:
    if sample.category == "Signal":
        return "#9ecae1"
    background_colors = ["#e41a1c", "#984ea3", "#ff7f00", "#4daf4a", "#a65628", "#f781bf"]
    return background_colors[background_index % len(background_colors)]


def add_publication_labels(ax, title: str, normalization: str) -> None:
    ax.text(0.06, 0.94, "SSCwf?", transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontstyle="italic", fontsize=14)
    ax.text(0.06, 0.875, r"LO simulation, $\sqrt{s}=40$ TeV", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    if title == "diphoton invariant mass":
        ax.text(0.06, 0.815, r"$H\rightarrow\gamma\gamma$", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    if normalization == "selected_xsec":
        ax.text(0.06, 0.065, r"Norm: $\sigma \times w \times \epsilon_{\gamma\gamma}$", transform=ax.transAxes, ha="left", va="bottom", fontsize=9)


def plot_histogram(
    title: str,
    samples: Sequence[SampleResult],
    output_dir: Path,
    normalization: str,
    density: bool,
) -> dict[str, str | float]:
    plt = ensure_matplotlib(output_dir)
    meta = PLOT_META.get(title, PlotMeta(slugify(title), title.title(), title, ""))

    reference = next(sample.histograms[title] for sample in samples if title in sample.histograms)
    x_values = reference.x
    width = bin_width(reference)
    edges = bin_edges(reference)
    bottoms = [0.0 for _ in x_values]

    fig, ax = plt.subplots(figsize=(6.8, 5.1))
    rows: list[dict[str, str | float]] = []
    background_index = 0

    for sample in ordered_samples(samples):
        if title not in sample.histograms:
            continue
        hist = sample.histograms[title]
        y_values = scaled_histogram(hist, sample, normalization, density)
        color = publication_color(sample, background_index)
        if sample.category != "Signal":
            background_index += 1

        tops = [bottom + value for bottom, value in zip(bottoms, y_values)]
        ax.fill_between(
            edges,
            step_values(bottoms),
            step_values(tops),
            step="post",
            label=sample.label,
            color=color,
            alpha=0.92 if sample.category == "Signal" else 0.96,
            linewidth=0.0,
        )
        ax.step(
            edges,
            step_values(tops),
            where="post",
            color="black" if sample.category != "Signal" else "#1f77b4",
            linewidth=0.7 if sample.category != "Signal" else 1.15,
            label="_nolegend_",
        )

        for bin_index, x in enumerate(x_values):
            rows.append(
                {
                    "plot": meta.slug,
                    "sample": sample.name,
                    "category": sample.category,
                    "x_center": x,
                    "raw_bin_weight": hist.y[bin_index],
                    "scaled_bin_value": y_values[bin_index],
                    "bin_width": width,
                    "cross_section_pb": sample.cross_section_pb,
                    "analysis_efficiency": sample.efficiency,
                    "selected_cross_section_pb": sample.selected_cross_section_pb,
                }
            )
        bottoms = tops

    if any(bottoms):
        ax.step(edges, step_values(bottoms), where="post", color="black", linewidth=1.05, label="_nolegend_")

    plot_x_min, plot_x_max = meta.x_range if meta.x_range else (reference.x_min, reference.x_max)
    ax.set_xlim(plot_x_min, plot_x_max)
    ax.set_xlabel(meta.x_label)
    if normalization == "unit_area":
        y_label = "Unit-normalized shape"
    elif density and meta.unit:
        y_label = rf"$\mathrm{{d}}\sigma_{{\rm sel}}/\mathrm{{d}}x$ [pb/{meta.unit}]"
    elif density:
        y_label = r"$\mathrm{d}\sigma_{\rm sel}/\mathrm{d}x$ [pb]"
    else:
        y_label = r"Bin cross section [pb]"
    ax.set_ylabel(y_label)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    visible_tops = [
        value
        for bin_index, value in enumerate(bottoms)
        if edges[bin_index] < plot_x_max and edges[bin_index + 1] > plot_x_min
    ]
    ymax = max(visible_tops) if visible_tops else (max(bottoms) if bottoms else 0.0)
    if ymax > 0:
        ax.set_ylim(0.0, ymax * 1.35)
    add_publication_labels(ax, title, normalization)
    ax.legend(frameon=False, loc="upper right", handlelength=1.6, borderaxespad=0.6)
    fig.tight_layout()

    plot_dir = output_dir / "plots"
    data_dir = output_dir / "data"
    plot_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    png_path = plot_dir / f"{meta.slug}.png"
    svg_path = plot_dir / f"{meta.slug}.svg"
    csv_path = data_dir / f"{meta.slug}.csv"
    fig.savefig(png_path, dpi=220)
    fig.savefig(svg_path)
    plt.close(fig)

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else ["plot"])
        writer.writeheader()
        writer.writerows(rows)

    return {
        "slug": meta.slug,
        "title": meta.display_name,
        "png": png_path.relative_to(output_dir).as_posix(),
        "svg": svg_path.relative_to(output_dir).as_posix(),
        "csv": csv_path.relative_to(output_dir).as_posix(),
        "x_min": plot_x_min,
        "x_max": plot_x_max,
        "unit": meta.unit,
    }


def write_summary_csv(samples: Sequence[SampleResult], output_dir: Path) -> str:
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / "sample_summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample",
                "category",
                "cross_section_pb",
                "cross_section_error_pb",
                "weight_scale",
                "events_read",
                "events_with_two_selected_photons",
                "sum_weight",
                "sum_diphoton_weight",
                "analysis_efficiency",
                "selected_cross_section_pb",
                "top_file",
                "dat_file",
            ],
        )
        writer.writeheader()
        for sample in samples:
            writer.writerow(
                {
                    "sample": sample.name,
                    "category": sample.category,
                    "cross_section_pb": sample.cross_section_pb,
                    "cross_section_error_pb": sample.cross_section_error_pb if sample.cross_section_error_pb is not None else "",
                    "weight_scale": sample.weight_scale,
                    "events_read": sample.events_read,
                    "events_with_two_selected_photons": sample.selected_events,
                    "sum_weight": sample.sum_weight,
                    "sum_diphoton_weight": sample.sum_diphoton_weight,
                    "analysis_efficiency": sample.efficiency,
                    "selected_cross_section_pb": sample.selected_cross_section_pb,
                    "top_file": sample.top_file,
                    "dat_file": sample.dat_file,
                }
            )
    return path.relative_to(output_dir).as_posix()


def write_html(
    output_dir: Path,
    title: str,
    run_tag: str,
    normalization: str,
    density: bool,
    samples: Sequence[SampleResult],
    plots: Sequence[dict[str, str | float]],
    summary_csv: str,
) -> Path:
    zip_name = "gammagamma_report_assets.zip"
    shutil.make_archive(str(output_dir / "gammagamma_report_assets"), "zip", root_dir=output_dir, base_dir="plots")

    sample_rows = "\n".join(
        f"""
        <tr>
          <td>{html.escape(sample.label)}</td>
          <td>{html.escape(sample.category)}</td>
          <td>{sample.cross_section_pb:.6g}</td>
          <td>{sample.weight_scale:.6g}</td>
          <td>{sample.efficiency:.6g}</td>
          <td>{sample.selected_cross_section_pb:.6g}</td>
          <td>{int(sample.events_read)}</td>
          <td>{int(sample.selected_events)}</td>
        </tr>
        """
        for sample in ordered_samples(samples)
    )

    plot_cards = "\n".join(
        f"""
        <section class="plot-card" id="{html.escape(str(plot['slug']))}">
          <div class="plot-head">
            <div>
              <h2>{html.escape(str(plot['title']))}</h2>
              <p>x range: {plot['x_min']:.6g} to {plot['x_max']:.6g}</p>
            </div>
            <div class="downloads">
              <a href="{html.escape(str(plot['png']))}" download>PNG</a>
              <a href="{html.escape(str(plot['svg']))}" download>SVG</a>
              <a href="{html.escape(str(plot['csv']))}" download>CSV</a>
            </div>
          </div>
          <a href="{html.escape(str(plot['png']))}" target="_blank">
            <img src="{html.escape(str(plot['png']))}" alt="{html.escape(str(plot['title']))}">
          </a>
        </section>
        """
        for plot in plots
    )

    normalization_text = {
        "selected_xsec": "Each plotted histogram area is normalized to cross section times weight scale times diphoton analysis efficiency for its sample.",
        "event_xsec": "Histogram bins are scaled by the event cross section and weight scale; selected-only histograms integrate to cross section times weight scale times efficiency.",
        "unit_area": "Each sample histogram is unit-area normalized before stacking.",
    }[normalization]
    if density:
        normalization_text += " Bin contents are divided by bin width."

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #1d2329;
      --muted: #65707c;
      --line: #d7dde4;
      --panel: #ffffff;
      --band: #f4f7fa;
      --accent: #b4232c;
    }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      color: var(--ink);
      background: var(--band);
    }}
    header {{
      padding: 28px 36px 22px;
      background: var(--panel);
      border-bottom: 1px solid var(--line);
    }}
    h1 {{
      margin: 0 0 8px;
      font-size: 26px;
      letter-spacing: 0;
    }}
    h2 {{
      margin: 0;
      font-size: 18px;
      letter-spacing: 0;
    }}
    p {{
      margin: 4px 0;
      color: var(--muted);
    }}
    main {{
      max-width: 1180px;
      margin: 0 auto;
      padding: 24px;
    }}
    .toolbar, .summary, .plot-card {{
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      margin-bottom: 18px;
    }}
    .toolbar {{
      padding: 14px 16px;
      display: flex;
      gap: 12px;
      align-items: center;
      flex-wrap: wrap;
    }}
    a {{
      color: var(--accent);
      text-decoration: none;
      font-weight: 600;
    }}
    .summary {{
      overflow-x: auto;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 14px;
    }}
    th, td {{
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      text-align: right;
      white-space: nowrap;
    }}
    th:first-child, td:first-child, th:nth-child(2), td:nth-child(2) {{
      text-align: left;
    }}
    th {{
      background: #eef2f6;
      font-weight: 700;
    }}
    .plot-card {{
      padding: 16px;
    }}
    .plot-head {{
      display: flex;
      align-items: flex-start;
      justify-content: space-between;
      gap: 16px;
      margin-bottom: 12px;
    }}
    .downloads {{
      display: flex;
      gap: 10px;
      flex-wrap: wrap;
    }}
    .downloads a, .toolbar a {{
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 6px 9px;
      background: #fff;
    }}
    img {{
      display: block;
      width: 100%;
      max-width: 920px;
      height: auto;
      margin: 0 auto;
    }}
  </style>
</head>
<body>
  <header>
    <h1>{html.escape(title)}</h1>
    <p>Run tag: {html.escape(run_tag)}. {html.escape(normalization_text)}</p>
  </header>
  <main>
    <nav class="toolbar">
      <a href="{html.escape(summary_csv)}" download>Download sample summary CSV</a>
      <a href="{html.escape(zip_name)}" download>Download plot images ZIP</a>
    </nav>
    <section class="summary">
      <table>
        <thead>
          <tr>
            <th>Sample</th>
            <th>Category</th>
            <th>&sigma; [pb]</th>
            <th>Weight scale</th>
            <th>Analysis efficiency</th>
            <th>&sigma;&times;&epsilon; [pb]</th>
            <th>Events read</th>
            <th>Diphoton events</th>
          </tr>
        </thead>
        <tbody>
          {sample_rows}
        </tbody>
      </table>
    </section>
    {plot_cards}
  </main>
</body>
</html>
"""
    index = output_dir / "index.html"
    index.write_text(html_text)
    return index


def main(argv: Sequence[str]) -> int:
    args = build_parser().parse_args(argv)
    analysis_root = args.analysis_root.resolve()
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = analysis_root / "plots" / f"gammagamma_{args.run_tag}"
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = load_samples(args)
    if not samples:
        die(f"no samples with .top/.dat files found under {analysis_root}")

    log(f"Loaded {len(samples)} samples")
    for sample in ordered_samples(samples):
        print(
            f"  {sample.name:18s} {sample.category:11s} "
            f"xsec={sample.cross_section_pb:.6g} pb weight={sample.weight_scale:.6g} "
            f"eff={sample.efficiency:.6g} xsec*weight*eff={sample.selected_cross_section_pb:.6g} pb"
        )

    titles = common_histogram_titles(samples)
    plots = [
        plot_histogram(title, samples, output_dir, args.normalization, density=not args.no_density)
        for title in titles
    ]
    summary_csv = write_summary_csv(samples, output_dir)
    index = write_html(
        output_dir,
        args.title,
        args.run_tag,
        args.normalization,
        density=not args.no_density,
        samples=samples,
        plots=plots,
        summary_csv=summary_csv,
    )
    log(f"Wrote report: {index}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
