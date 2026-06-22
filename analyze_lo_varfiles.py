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


SUMMARY_FIELDS = [
    "sample",
    "category",
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
    if lowered in {"null", "none"}:
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


def write_html_report(metadata: dict[str, Any], rows: Sequence[dict[str, Any]], output_dir: Path, assets: Sequence[Path] = ()) -> Path:
    def fmt(value: Any) -> str:
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)

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


def load_analysis_inputs(config_path: Path) -> tuple[dict[str, Any], Path, str, str, float, list[SampleInfo], Path]:
    config = load_config(config_path)
    analysis = config["analysis"]
    name = str(analysis.get("name", "analysis"))
    run_tag = str(analysis.get("run_tag", os.environ.get("RUN_TAG", "run_01")))
    luminosity_fb = float(analysis.get("luminosity_fb", DEFAULT_LUMINOSITY_FB))
    analysis_root = Path(analysis.get("analysis_root", DEFAULT_ANALYSIS_ROOT)).expanduser()
    rate_factors = parse_rate_factors(analysis.get("rate_factors"))
    samples = discover_samples(analysis_root, run_tag, requested_samples(analysis.get("samples")), rate_factors)
    if not samples:
        raise RuntimeError(f"no gamma-gamma _var.root samples found under {analysis_root} for run tag {run_tag}")
    output_dir = output_dir_for(analysis_root, run_tag, name, analysis.get("output_dir"))
    return analysis, analysis_root, name, run_tag, luminosity_fb, samples, output_dir


def run_cuts(config_path: Path) -> Path:
    analysis, _, name, run_tag, luminosity_fb, samples, output_dir = load_analysis_inputs(config_path)
    cuts = [Cut.from_mapping(item) for item in analysis.get("cuts", [])]
    max_events = analysis.get("max_events")
    rows = [summarize_sample(sample, cuts, luminosity_fb, _optional_float(max_events)) for sample in samples]
    cut_metadata = [
        {"variable": cut.variable, "min": cut.minimum, "max": cut.maximum}
        for cut in cuts
    ]
    metadata = {
        "name": name,
        "mode": "cuts",
        "run_tag": run_tag,
        "luminosity_fb": luminosity_fb,
        "cuts": cut_metadata,
    }
    return write_outputs(metadata, rows, output_dir)


def run_xgboost(config_path: Path) -> Path:
    analysis, _, name, run_tag, luminosity_fb, samples, output_dir = load_analysis_inputs(config_path)
    signal_samples = [sample for sample in samples if sample.category == "Signal"]
    background_samples = [sample for sample in samples if sample.category != "Signal"]
    if not signal_samples or not background_samples:
        raise RuntimeError("xgboost mode requires at least one Signal sample and one background sample")

    try:
        from xgboost_root_varfiles_module import run_signal_background_analysis
    except ImportError as exc:
        raise RuntimeError(
            "xgboost mode requires optional packages: xgboost, scikit-learn, and tqdm. "
            "Install them before running this subcommand."
        ) from exc

    xgb_cfg = analysis.get("xgboost", {}) or {}
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
        signal_metadata=[{"sample": sample.name, "category": sample.category} for sample in signal_samples],
        background_metadata=[{"sample": sample.name, "category": sample.category} for sample in background_samples],
        luminosity=luminosity_fb,
        test_size=float(xgb_cfg.get("test_size", 0.35)),
        seed=int(xgb_cfg.get("seed", 12345)),
        systematics=float(xgb_cfg.get("systematics", 0.0)),
        max_events=xgb_cfg.get("max_events"),
        model_params=xgb_cfg.get("model_params"),
    )
    metadata = {
        "name": name,
        "mode": "xgboost",
        "run_tag": run_tag,
        "luminosity_fb": luminosity_fb,
        "cuts": [{"variable": "xgboost_score", "min": result["metadata"]["best_threshold"], "max": None}],
        "xgboost": result["metadata"],
    }
    rows = result["summary_rows"]
    assets = [Path(path) for path in result["metadata"].get("outputs", {}).values() if str(path).endswith((".png", ".json", ".csv"))]
    return write_outputs(metadata, rows, output_dir, assets)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("cuts", "xgboost"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "cuts":
        index = run_cuts(args.config)
    elif args.command == "xgboost":
        index = run_xgboost(args.config)
    else:  # pragma: no cover - argparse prevents this.
        raise RuntimeError(f"unknown command {args.command}")
    print(f"Wrote analysis report: {index}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        raise SystemExit(f"ERROR: {error}") from error
