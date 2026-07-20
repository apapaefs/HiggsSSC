#!/usr/bin/env python3
"""Analyze gamma-gamma LO campaign ``_var.root`` files.

Examples
--------
Run a rectangular-cut analysis from a YAML card:

    python3 analyze_lo_varfiles.py cuts --config cuts.yaml

Train and score an optional XGBoost analysis from a YAML card:

    python3 analyze_lo_varfiles.py xgboost --config xgboost.yaml
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from read_root_varfiles import FEATURE_NAMES, read_named_ROOT_varfile


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_ANALYSIS_ROOT = REPO_ROOT / "hgammagamma" / "LOAnalysis"
DEFAULT_LUMINOSITY_FB = 100.0
SENTINEL_VALUE = -900.0
YR4_BR_H_TO_GAMMAGAMMA = 2.27e-3
SIGNAL_GGH_K_FACTOR = 2.0
SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT = SIGNAL_GGH_K_FACTOR * YR4_BR_H_TO_GAMMAGAMMA
DEFAULT_SAMPLE_RATE_FACTORS = {
    "signal_gg_h_aa": SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT,
}


@dataclass(frozen=True)
class Cut:
    variable: str
    minimum: float | None = None
    maximum: float | None = None

    def __post_init__(self) -> None:
        if self.variable not in FEATURE_NAMES:
            raise ValueError(f"unknown cut variable '{self.variable}'. Known variables: {', '.join(FEATURE_NAMES)}")
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError(f"cut minimum is greater than maximum for {self.variable}")

    @classmethod
    def from_mapping(cls, mapping: dict[str, Any]) -> "Cut":
        return cls(
            variable=str(mapping["variable"]),
            minimum=_optional_float(mapping.get("min")),
            maximum=_optional_float(mapping.get("max")),
        )

    def accepts(self, row: dict[str, float]) -> bool:
        value = float(row[self.variable])
        if not math.isfinite(value) or value <= SENTINEL_VALUE:
            return False
        if self.minimum is not None and value < self.minimum:
            return False
        if self.maximum is not None and value > self.maximum:
            return False
        return True


@dataclass
class SampleInfo:
    name: str
    category: str
    sample_dir: Path
    var_file: Path
    dat_file: Path
    cross_section_pb: float
    cross_section_error_pb: float | None
    weight_scale: float
    events_read: float
    sum_weight: float
    analysis_name: str = "legacy_direct_photons"
    detector_response: str = "none"
    response_mode: str = "genuine"


@dataclass
class AnalysisRunResult:
    index_html: Path
    output_dir: Path
    metadata: dict[str, Any]
    rows: list[dict[str, Any]]
    assets: list[Path]


class ProgressBar:
    def __init__(
        self,
        total: int,
        label: str,
        stream: Any = None,
        enabled: bool = True,
        width: int = 28,
    ) -> None:
        self.total = max(int(total), 0)
        self.label = label
        self.stream = stream if stream is not None else sys.stderr
        self.enabled = enabled
        self.width = max(int(width), 1)
        self._finished = False

    def update(self, current: int, item: str | None = None) -> None:
        if not self.enabled:
            return
        denominator = max(self.total, 1)
        current = max(0, min(int(current), denominator))
        filled = int(round(self.width * current / denominator))
        bar = "#" * filled + "-" * (self.width - filled)
        suffix = f" {item}" if item else ""
        self.stream.write(f"\r{self.label} [{bar}] {current}/{self.total}{suffix}")
        self.stream.flush()

    def finish(self) -> None:
        if not self.enabled or self._finished:
            return
        self.stream.write("\n")
        self.stream.flush()
        self._finished = True


SUMMARY_FIELDS = [
    "sample",
    "category",
    "analysis",
    "detector_response",
    "response_mode",
    "input_file",
    "raw_cross_section_pb",
    "cross_section_pb",
    "weight_scale",
    "entries_read",
    "selected_entries",
    "sum_weight",
    "sum_selected_weight",
    "analysis_efficiency",
    "selected_cross_section_pb",
    "expected_events",
    "mc_events_after_analysis",
]


def _optional_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _parse_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if value[0:1] in {"'", '"'} and value[-1:] == value[0]:
        return value[1:-1]
    lowered = value.lower()
    if lowered in {"true", "yes"}:
        return True
    if lowered in {"false", "no"}:
        return False
    # YAML null spellings are ``null`` and ``~``.  Keep bare ``none`` as a
    # string because it is the detector-response profile used by this project.
    if lowered in {"null", "~"}:
        return None
    try:
        if re.search(r"[.eE]", value):
            return float(value)
        return int(value)
    except ValueError:
        return value


def _strip_yaml_comment(line: str) -> str:
    in_quote = None
    for index, char in enumerate(line):
        if char in {"'", '"'}:
            in_quote = None if in_quote == char else char
        if char == "#" and in_quote is None:
            return line[:index]
    return line


def _parse_simple_yaml(text: str) -> dict[str, Any]:
    """Parse the small YAML subset used by the analysis examples."""

    root: dict[str, Any] = {}
    current_section: dict[str, Any] | None = None
    current_nested: dict[str, Any] | None = None
    current_list: list[dict[str, Any]] | None = None
    current_item: dict[str, Any] | None = None

    for raw_line in text.splitlines():
        line = _strip_yaml_comment(raw_line).rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()

        if indent == 0:
            key, _, value = stripped.partition(":")
            if value.strip():
                root[key] = _parse_scalar(value)
                current_section = None
            else:
                current_section = {}
                root[key] = current_section
            current_nested = None
            current_list = None
            current_item = None
            continue

        if current_section is None:
            raise ValueError(f"unsupported YAML structure near: {raw_line}")

        if indent == 2 and stripped.startswith("- "):
            raise ValueError(f"top-level lists are not supported near: {raw_line}")

        if indent == 2:
            key, _, value = stripped.partition(":")
            if not _:
                raise ValueError(f"expected key/value YAML line near: {raw_line}")
            if value.strip():
                current_section[key] = _parse_scalar(value)
                current_nested = None
                current_list = None
            elif key == "cuts":
                current_list = []
                current_section[key] = current_list
                current_nested = None
            else:
                current_nested = {}
                current_section[key] = current_nested
                current_list = None
            current_item = None
            continue

        if indent == 4 and current_list is not None and stripped.startswith("- "):
            current_item = {}
            current_list.append(current_item)
            item_text = stripped[2:].strip()
            if item_text:
                key, _, value = item_text.partition(":")
                if not _:
                    raise ValueError(f"expected list item key/value near: {raw_line}")
                current_item[key] = _parse_scalar(value)
            continue

        if indent >= 4 and current_item is not None:
            key, _, value = stripped.partition(":")
            if not _:
                raise ValueError(f"expected list item key/value near: {raw_line}")
            current_item[key] = _parse_scalar(value)
            continue

        if indent >= 4 and current_nested is not None:
            key, _, value = stripped.partition(":")
            if not _:
                raise ValueError(f"expected nested key/value near: {raw_line}")
            current_nested[key] = _parse_scalar(value)
            continue

        raise ValueError(f"unsupported YAML structure near: {raw_line}")

    return root


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text()
    try:
        import yaml  # type: ignore
    except ImportError:
        config = _parse_simple_yaml(text)
    else:
        config = yaml.safe_load(text) or {}
    if "analysis" not in config or not isinstance(config["analysis"], dict):
        raise ValueError("config must contain an 'analysis' mapping")
    return config


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


def response_provenance(dat: dict[str, float | str]) -> tuple[str, str, str]:
    """Return analysis, detector profile, and response mode for new or historical files."""

    analysis_name = str(dat.get("analysis", "")).strip()
    explicit_response = str(dat.get("detector_response", "")).strip()
    if explicit_response:
        detector_response = explicit_response
    elif analysis_name == "SSC_GEM_weighted_response":
        detector_response = "ssc"
    elif not analysis_name or analysis_name == "legacy_direct_photons":
        detector_response = "none"
    else:
        detector_response = "unknown"
    if not analysis_name:
        analysis_name = "legacy_direct_photons"
    default_mode = "genuine" if detector_response in {"ssc", "none"} else "unknown"
    response_mode = str(dat.get("response_mode", default_mode))
    return analysis_name, detector_response, response_mode


def validate_detector_response(samples: Sequence[SampleInfo], expected: Any) -> str:
    actual = sorted({sample.detector_response for sample in samples})
    if expected is None:
        return actual[0] if len(actual) == 1 else ",".join(actual)
    expected_response = str(expected).strip().lower()
    if expected_response not in {"ssc", "none"}:
        raise ValueError("analysis.detector_response must be 'ssc' or 'none'")
    mismatched = [sample.name for sample in samples if sample.detector_response != expected_response]
    if mismatched:
        raise ValueError(
            f"analysis.detector_response={expected_response} does not match .dat metadata for "
            + ", ".join(mismatched)
        )
    return expected_response


def parse_rate_factors(value: Any) -> dict[str, float]:
    if value is None:
        return {}
    if not isinstance(value, dict):
        raise ValueError("analysis.rate_factors must be a mapping of sample or category names to numeric factors")
    return {str(key): float(factor) for key, factor in value.items() if factor is not None}


def resolve_weight_scale(
    sample_name: str,
    category: str,
    dat_weight_scale: float,
    configured_rate_factors: dict[str, float] | None,
) -> float:
    """Choose the total rate factor used above the LO campaign files.

    Precedence is sample-specific YAML, category YAML, nontrivial `.dat`
    metadata, then the known gamma-gamma signal default for legacy files.
    """

    configured_rate_factors = configured_rate_factors or {}
    if sample_name in configured_rate_factors:
        return float(configured_rate_factors[sample_name])
    if category in configured_rate_factors:
        return float(configured_rate_factors[category])
    dat_weight_scale = float(dat_weight_scale)
    if not math.isclose(dat_weight_scale, 1.0, rel_tol=0.0, abs_tol=1e-12):
        return dat_weight_scale
    return float(DEFAULT_SAMPLE_RATE_FACTORS.get(sample_name, dat_weight_scale))


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
        text = text.replace("&plusmn;", "+/-").replace("&#177;", "+/-")
        match = re.search(r"s=\s*([-+0-9.eE]+)\s*\+/-\s*([-+0-9.eE]+)\s*\(pb\)", text)
        if match:
            return float(match.group(1)), float(match.group(2))

    crossx = sample_dir / "mg5_process" / "crossx.html"
    if crossx.exists():
        text = crossx.read_text(errors="ignore")
        match = re.search(r"results\.html\">\s*([-+0-9.eE]+)\s*<font[^>]*>.*?</font>\s*([-+0-9.eE]+)", text)
        if match:
            return float(match.group(1)), float(match.group(2))

    raise FileNotFoundError(f"could not find MG5 cross section for {sample_dir}")


def discover_samples(
    analysis_root: Path,
    run_tag: str,
    requested: set[str] | None = None,
    rate_factors: dict[str, float] | None = None,
) -> list[SampleInfo]:
    samples: list[SampleInfo] = []
    for category in ("Backgrounds", "Signal"):
        events_dir = analysis_root / category / "events"
        if not events_dir.exists():
            continue
        for sample_dir in sorted(path for path in events_dir.iterdir() if path.is_dir()):
            if requested is not None and sample_dir.name not in requested:
                continue
            var_files = sorted(sample_dir.glob(f"*-{run_tag}_var.root"))
            if not var_files:
                continue
            var_file = var_files[-1]
            dat_file = Path(str(var_file).replace("_var.root", ".dat"))
            if not dat_file.exists():
                raise FileNotFoundError(f"missing .dat summary for {var_file}")
            dat = parse_key_value_dat(dat_file)
            analysis_name, detector_response, response_mode = response_provenance(dat)
            cross_section_pb, cross_section_error_pb = parse_cross_section(sample_dir, run_tag)
            dat_weight_scale = float(dat.get("weight_scale", 1.0))
            samples.append(
                SampleInfo(
                    name=sample_dir.name,
                    category=category,
                    sample_dir=sample_dir,
                    var_file=var_file,
                    dat_file=dat_file,
                    cross_section_pb=cross_section_pb,
                    cross_section_error_pb=cross_section_error_pb,
                    weight_scale=resolve_weight_scale(sample_dir.name, category, dat_weight_scale, rate_factors),
                    events_read=float(dat.get("events_read", 0.0)),
                    sum_weight=float(dat.get("sum_weight", 0.0)),
                    analysis_name=analysis_name,
                    detector_response=detector_response,
                    response_mode=response_mode,
                )
            )
    return samples


def apply_cuts(rows: Sequence[dict[str, float]], cuts: Sequence[Cut]) -> list[bool]:
    return [all(cut.accepts(row) for cut in cuts) for row in rows]


def summarize_sample(sample: SampleInfo, cuts: Sequence[Cut], luminosity_fb: float, max_events: int | None = None) -> dict[str, Any]:
    rows, weights = read_named_ROOT_varfile(sample.var_file, max_events=max_events)
    decisions = apply_cuts(rows, cuts)
    sum_weight = float(sum(weights))
    sum_selected_weight = float(sum(weight for weight, selected in zip(weights, decisions) if selected))
    selected_entries = int(sum(1 for selected in decisions if selected))
    efficiency = sum_selected_weight / sum_weight if sum_weight > 0.0 else 0.0
    cross_section_pb = sample.cross_section_pb * sample.weight_scale
    selected_cross_section_pb = cross_section_pb * efficiency
    return {
        "sample": sample.name,
        "category": sample.category,
        "analysis": sample.analysis_name,
        "detector_response": sample.detector_response,
        "response_mode": sample.response_mode,
        "input_file": str(sample.var_file),
        "raw_cross_section_pb": sample.cross_section_pb,
        "cross_section_pb": cross_section_pb,
        "weight_scale": sample.weight_scale,
        "entries_read": len(rows),
        "selected_entries": selected_entries,
        "sum_weight": sum_weight,
        "sum_selected_weight": sum_selected_weight,
        "analysis_efficiency": efficiency,
        "selected_cross_section_pb": selected_cross_section_pb,
        "expected_events": selected_cross_section_pb * luminosity_fb * 1000.0,
        "mc_events_after_analysis": selected_entries,
    }


def output_dir_for(analysis_root: Path, run_tag: str, name: str, configured: Any = None) -> Path:
    if configured:
        return Path(configured).expanduser()
    return analysis_root / "analyses" / run_tag / name


def write_summary_csv(rows: Sequence[dict[str, Any]], output_dir: Path) -> Path:
    path = output_dir / "summary.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})
    return path


def write_summary_json(metadata: dict[str, Any], rows: Sequence[dict[str, Any]], output_dir: Path) -> Path:
    path = output_dir / "summary.json"
    path.write_text(json.dumps({"metadata": metadata, "samples": list(rows)}, indent=2))
    return path


def compute_analysis_totals(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    signal_rows = [row for row in rows if row.get("category") == "Signal"]
    background_rows = [row for row in rows if row.get("category") != "Signal"]
    signal_expected = sum(float(row.get("expected_events", 0.0)) for row in signal_rows)
    background_expected = sum(float(row.get("expected_events", 0.0)) for row in background_rows)
    significance = signal_expected / math.sqrt(background_expected) if background_expected > 0.0 else None
    return {
        "signal_selected_cross_section_pb": sum(float(row.get("selected_cross_section_pb", 0.0)) for row in signal_rows),
        "background_selected_cross_section_pb": sum(float(row.get("selected_cross_section_pb", 0.0)) for row in background_rows),
        "signal_expected_events": signal_expected,
        "background_expected_events": background_expected,
        "signal_mc_events_after_analysis": sum(int(row.get("mc_events_after_analysis", 0)) for row in signal_rows),
        "background_mc_events_after_analysis": sum(int(row.get("mc_events_after_analysis", 0)) for row in background_rows),
        "approx_significance_s_over_sqrt_b": significance,
    }


def write_html_report(metadata: dict[str, Any], rows: Sequence[dict[str, Any]], output_dir: Path, assets: Sequence[Path] = ()) -> Path:
    def fmt(value: Any) -> str:
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)

    def fmt_optional(value: Any) -> str:
        if value is None:
            return "undefined (B <= 0)"
        return fmt(value)

    totals = metadata.get("totals")
    if not isinstance(totals, dict):
        totals = compute_analysis_totals(rows)
    row_html = "\n".join(
        "<tr>"
        + "".join(f"<td>{html.escape(fmt(row.get(field, '')))}</td>" for field in SUMMARY_FIELDS)
        + "</tr>"
        for row in rows
    )
    asset_links = "\n".join(
        f'<li><a href="{html.escape(path.relative_to(output_dir).as_posix())}">{html.escape(path.name)}</a></li>'
        for path in assets
        if path.exists()
    )
    cuts = metadata.get("cuts", [])
    cuts_html = "\n".join(
        f"<li>{html.escape(cut['variable'])}: "
        f"{'-inf' if cut.get('min') is None else html.escape(fmt(cut.get('min')))} to "
        f"{'inf' if cut.get('max') is None else html.escape(fmt(cut.get('max')))}</li>"
        for cut in cuts
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{html.escape(metadata['name'])}</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem; color: #1f2933; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
    th, td {{ border-bottom: 1px solid #d9e2ec; padding: 0.45rem 0.55rem; text-align: right; }}
    th:first-child, td:first-child, th:nth-child(2), td:nth-child(2), th:nth-child(3), td:nth-child(3) {{ text-align: left; }}
    th {{ background: #f0f4f8; }}
    .meta {{ display: flex; gap: 2rem; flex-wrap: wrap; margin: 1rem 0; }}
    a {{ color: #1f5eff; }}
  </style>
</head>
<body>
  <h1>{html.escape(metadata['name'])}</h1>
  <div class="meta">
    <div>Mode: {html.escape(metadata['mode'])}</div>
    <div>Run tag: {html.escape(metadata['run_tag'])}</div>
    <div>Detector response: {html.escape(str(metadata.get('detector_response', 'unspecified')))}</div>
    <div>Luminosity: {metadata['luminosity_fb']:.6g} fb<sup>-1</sup></div>
  </div>
  <p><a href="summary.csv">Download CSV</a> | <a href="summary.json">Download JSON</a></p>
  <h2>Selection</h2>
  <ul>{cuts_html}</ul>
  <h2>Samples</h2>
  <table>
    <thead><tr>{"".join(f"<th>{html.escape(field)}</th>" for field in SUMMARY_FIELDS)}</tr></thead>
    <tbody>{row_html}</tbody>
  </table>
  <h2>Totals</h2>
  <table>
    <tbody>
      <tr><th>signal_expected_events</th><td>{html.escape(fmt(totals['signal_expected_events']))}</td></tr>
      <tr><th>background_expected_events</th><td>{html.escape(fmt(totals['background_expected_events']))}</td></tr>
      <tr><th>approx_significance_s_over_sqrt_b</th><td>{html.escape(fmt_optional(totals['approx_significance_s_over_sqrt_b']))}</td></tr>
    </tbody>
  </table>
  <h2>Additional outputs</h2>
  <ul>{asset_links}</ul>
</body>
</html>
"""
    path = output_dir / "index.html"
    path.write_text(html_text)
    return path


def write_outputs(metadata: dict[str, Any], rows: Sequence[dict[str, Any]], output_dir: Path, assets: Sequence[Path] = ()) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_summary_csv(rows, output_dir)
    write_summary_json(metadata, rows, output_dir)
    return write_html_report(metadata, rows, output_dir, assets)


def requested_samples(value: Any) -> set[str] | None:
    if value is None or value == "all":
        return None
    if isinstance(value, str):
        return {part.strip() for part in value.split(",") if part.strip()}
    return {str(part) for part in value}


def resolve_run_tag(analysis: dict[str, Any], run_tag_override: str | None = None) -> str:
    if run_tag_override:
        return str(run_tag_override)
    return str(analysis.get("run_tag", os.environ.get("RUN_TAG", "run_01")))


def build_terminal_summary(
    metadata: dict[str, Any],
    rows: Sequence[dict[str, Any]],
    output_dir: Path,
    assets: Sequence[Path] = (),
) -> str:
    def fmt(value: float) -> str:
        return f"{float(value):.6g}"

    def category_line(label: str, selected_rows: Sequence[dict[str, Any]]) -> str:
        selected_xsec = sum(float(row.get("selected_cross_section_pb", 0.0)) for row in selected_rows)
        expected = sum(float(row.get("expected_events", 0.0)) for row in selected_rows)
        mc_events = sum(int(row.get("mc_events_after_analysis", 0)) for row in selected_rows)
        return (
            f"  {label}: samples={len(selected_rows)} "
            f"selected_xsec_pb={fmt(selected_xsec)} expected events={fmt(expected)} "
            f"mc_events={mc_events}"
        )

    totals = metadata.get("totals")
    if not isinstance(totals, dict):
        totals = compute_analysis_totals(rows)
    significance = totals["approx_significance_s_over_sqrt_b"]
    significance_text = fmt(significance) if significance is not None else "undefined (B <= 0)"
    signal_rows = [row for row in rows if row.get("category") == "Signal"]
    background_rows = [row for row in rows if row.get("category") != "Signal"]
    output_dir = Path(output_dir)
    lines = [
        "",
        "Analysis summary",
        f"  name: {metadata.get('name', 'analysis')}",
        f"  mode: {metadata.get('mode', 'analysis')}",
        f"  run tag: {metadata.get('run_tag', '')}",
        f"  detector response: {metadata.get('detector_response', 'unspecified')}",
        f"  luminosity: {fmt(float(metadata.get('luminosity_fb', 0.0)))} fb^-1",
        f"  samples: {len(rows)}",
        category_line("Signal", signal_rows),
        category_line("Backgrounds", background_rows),
        f"  approx_significance_s_over_sqrt_b={significance_text}",
        "  output files:",
        f"    {output_dir / 'summary.csv'}",
        f"    {output_dir / 'summary.json'}",
        f"    {output_dir / 'index.html'}",
    ]
    for asset in assets:
        lines.append(f"    {asset}")
    return "\n".join(lines)


def load_analysis_inputs(
    config_path: Path,
    run_tag_override: str | None = None,
) -> tuple[dict[str, Any], Path, str, str, float, list[SampleInfo], Path]:
    config = load_config(config_path)
    analysis = config["analysis"]
    name = str(analysis.get("name", "analysis"))
    run_tag = resolve_run_tag(analysis, run_tag_override)
    luminosity_fb = float(analysis.get("luminosity_fb", DEFAULT_LUMINOSITY_FB))
    analysis_root = Path(analysis.get("analysis_root", DEFAULT_ANALYSIS_ROOT)).expanduser()
    rate_factors = parse_rate_factors(analysis.get("rate_factors"))
    samples = discover_samples(analysis_root, run_tag, requested_samples(analysis.get("samples")), rate_factors)
    if not samples:
        raise RuntimeError(f"no gamma-gamma _var.root samples found under {analysis_root} for run tag {run_tag}")
    analysis["_resolved_detector_response"] = validate_detector_response(
        samples, analysis.get("detector_response")
    )
    output_dir = output_dir_for(analysis_root, run_tag, name, analysis.get("output_dir"))
    return analysis, analysis_root, name, run_tag, luminosity_fb, samples, output_dir


def run_cuts(config_path: Path, run_tag_override: str | None = None, progress_enabled: bool = True) -> AnalysisRunResult:
    analysis, _, name, run_tag, luminosity_fb, samples, output_dir = load_analysis_inputs(config_path, run_tag_override)
    cuts = [Cut.from_mapping(item) for item in analysis.get("cuts", [])]
    max_events = analysis.get("max_events")
    rows: list[dict[str, Any]] = []
    progress = ProgressBar(len(samples), "Analyzing samples", enabled=progress_enabled)
    try:
        for index, sample in enumerate(samples, start=1):
            rows.append(summarize_sample(sample, cuts, luminosity_fb, _optional_float(max_events)))
            progress.update(index, sample.name)
    finally:
        progress.finish()
    cut_metadata = [
        {"variable": cut.variable, "min": cut.minimum, "max": cut.maximum}
        for cut in cuts
    ]
    metadata = {
        "name": name,
        "mode": "cuts",
        "run_tag": run_tag,
        "luminosity_fb": luminosity_fb,
        "detector_response": analysis["_resolved_detector_response"],
        "cuts": cut_metadata,
    }
    metadata["totals"] = compute_analysis_totals(rows)
    index = write_outputs(metadata, rows, output_dir)
    return AnalysisRunResult(index_html=index, output_dir=output_dir, metadata=metadata, rows=rows, assets=[])


def run_xgboost(config_path: Path, run_tag_override: str | None = None, progress_enabled: bool = True) -> AnalysisRunResult:
    analysis, _, name, run_tag, luminosity_fb, samples, output_dir = load_analysis_inputs(config_path, run_tag_override)
    signal_samples = [sample for sample in samples if sample.category == "Signal"]
    background_samples = [sample for sample in samples if sample.category != "Signal"]
    if not signal_samples or not background_samples:
        raise RuntimeError("xgboost mode requires at least one Signal sample and one background sample")

    progress = ProgressBar(3, "XGBoost analysis", enabled=progress_enabled)
    try:
        progress.update(1, "prepared samples")
        from xgboost_root_varfiles_module import run_signal_background_analysis
    except ImportError as exc:
        raise RuntimeError(
            "xgboost mode requires optional packages: xgboost, scikit-learn, and tqdm. "
            "Install them before running this subcommand."
        ) from exc
    finally:
        if "run_signal_background_analysis" not in locals():
            progress.finish()

    xgb_cfg = analysis.get("xgboost", {}) or {}
    try:
        result = run_signal_background_analysis(
            signal_files=[sample.var_file for sample in signal_samples],
            background_files=[sample.var_file for sample in background_samples],
            output_dir=output_dir,
            signal_xsecs_fb=[sample.cross_section_pb * 1000.0 for sample in signal_samples],
            background_xsecs_fb=[sample.cross_section_pb * 1000.0 for sample in background_samples],
            signal_rate_factors=[sample.weight_scale for sample in signal_samples],
            background_rate_factors=[sample.weight_scale for sample in background_samples],
            signal_normalisation_weights=[sample.sum_weight for sample in signal_samples],
            background_normalisation_weights=[sample.sum_weight for sample in background_samples],
            signal_metadata=[
                {
                    "sample": sample.name,
                    "category": sample.category,
                    "analysis": sample.analysis_name,
                    "detector_response": sample.detector_response,
                    "response_mode": sample.response_mode,
                }
                for sample in signal_samples
            ],
            background_metadata=[
                {
                    "sample": sample.name,
                    "category": sample.category,
                    "analysis": sample.analysis_name,
                    "detector_response": sample.detector_response,
                    "response_mode": sample.response_mode,
                }
                for sample in background_samples
            ],
            luminosity=luminosity_fb,
            test_size=float(xgb_cfg.get("test_size", 0.35)),
            seed=int(xgb_cfg.get("seed", 12345)),
            systematics=float(xgb_cfg.get("systematics", 0.0)),
            max_events=xgb_cfg.get("max_events"),
            model_params=xgb_cfg.get("model_params"),
        )
        progress.update(2, "trained and scored")
        metadata = {
            "name": name,
            "mode": "xgboost",
            "run_tag": run_tag,
            "luminosity_fb": luminosity_fb,
            "detector_response": analysis["_resolved_detector_response"],
            "cuts": [{"variable": "xgboost_score", "min": result["metadata"]["best_threshold"], "max": None}],
            "xgboost": result["metadata"],
        }
        rows = result["summary_rows"]
        metadata["totals"] = compute_analysis_totals(rows)
        assets = [Path(path) for path in result["metadata"].get("outputs", {}).values() if str(path).endswith((".png", ".json", ".csv"))]
        index = write_outputs(metadata, rows, output_dir, assets)
        progress.update(3, "wrote outputs")
        return AnalysisRunResult(index_html=index, output_dir=output_dir, metadata=metadata, rows=list(rows), assets=assets)
    finally:
        progress.finish()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("cuts", "xgboost"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", type=Path, required=True)
        subparser.add_argument("--run-tag", help="Override analysis.run_tag from the YAML card.")
        subparser.add_argument("--no-progress", action="store_true", help="Disable terminal progress bars.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "cuts":
        result = run_cuts(args.config, run_tag_override=args.run_tag, progress_enabled=not args.no_progress)
    elif args.command == "xgboost":
        result = run_xgboost(args.config, run_tag_override=args.run_tag, progress_enabled=not args.no_progress)
    else:  # pragma: no cover - argparse prevents this.
        raise RuntimeError(f"unknown command {args.command}")
    print(build_terminal_summary(result.metadata, result.rows, result.output_dir, result.assets))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        raise SystemExit(f"ERROR: {error}") from error
