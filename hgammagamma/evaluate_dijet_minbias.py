"""Evaluate SSC diphoton double-fake envelopes in inclusive Herwig events."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hgammagamma.dijet_response import load_model, response_rows


def generated_cross_section(out):
    """Read the unit-weight ThePEG event count, cross section and its error in nb."""
    table = Path(out).read_text()
    if "(nb)" not in table or "Events carry unit weights." not in table:
        raise ValueError("expected unit-weight events and a cross section in nb")
    match = re.search(
        r"Total \(from generated events\):\s*(\d+)\s+\d+\s+([\d.]+)\((\d+)\)(?:e([+-]?\d+))?",
        table,
    )
    if match is None:
        raise ValueError("missing generated-event cross section")
    count, mantissa, last_digit_error, exponent = match.groups()
    power = int(exponent or 0)
    decimals = len(mantissa.partition(".")[2])
    return int(count), float(mantissa) * 10**power, int(last_digit_error) * 10**(power-decimals)


def evaluate(directories, binary, model=None):
    rows = []
    seeds = set()
    for directory in map(Path, directories):
        card = (directory / "minbias.in").read_text()
        settings = [line.strip() for line in card.splitlines() if not line.lstrip().startswith("#")]
        if not {
            "set EventGenerator:EventHandler:LuminosityFunction:Energy 40000.0",
            "read snippets/MB.in",
            "read snippets/Diffraction.in",
        } <= set(settings):
            raise ValueError("expected 40 TeV inclusive minimum-bias card")
        seed = int(re.search(r"^set EventGenerator:RandomNumberGenerator:Seed (\d+)$", card, re.M)[1])
        if seed in seeds:
            raise ValueError("repeated generator seed")
        seeds.add(seed)
        n, sigma_nb, err_nb = generated_cross_section(directory / "minbias.out")
        count = int(re.search(r"^set EventGenerator:NumberOfEvents (\d+)$", card, re.M)[1])
        if n != count or n <= 0 or "A root tree has been written" not in (directory / "run.log").read_text():
            raise ValueError("incomplete native minimum-bias run")
        roots = list((directory / "events").glob("*.root"))
        if len(roots) != 1:
            raise ValueError("expected one finalized ROOT file")

        estimates = {}
        for entry in response_rows(binary, roots[0], seed + 100, model):
            if int(entry["events_read"]) != n:
                raise ValueError("ROOT entries differ from generated events")
            key = entry["scenario"] + "@" + entry["fake_scale"]
            estimates[key] = dict(
                rate_fb=sigma_nb * 1.e6 * float(entry["sum_probability"]) / n,
                contributing_sources=int(entry["contributing_events"]),
            )
        rows.append(dict(directory=str(directory), seed=seed, generated=n, generated_cross_section_nb=sigma_nb,
                         sigma_error_nb=err_nb, estimates=estimates))
    if len(rows) < 2 or len({r["generated"] for r in rows}) != 1:
        raise ValueError("supply at least two independent equal-size samples")

    totals = {}
    for key in rows[0]["estimates"]:
        series = [row["estimates"][key] for row in rows]
        rates = [value["rate_fb"] for value in series]
        totals[key] = dict(rate_fb=statistics.mean(rates),
                           mc_error_fb=statistics.stdev(rates) / math.sqrt(len(rows)),
                           contributing_sources=sum(value["contributing_sources"] for value in series))
    return dict(
        description="Herwig minimum-bias model at 40 TeV including diffraction. Each inclusive event is counted once. "
                    "Each seed is normalized with its own event-handler cross section. Errors are standard errors "
                    "between independent runs.",
        window_GeV=[120, 130],
        frozen_selection="Pair passes the saved classifier threshold across the full accepted mass range.",
        analyzer_sha256=hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        generated_events=sum(row["generated"] for row in rows),
        totals=totals, frozen_model=model[1] if model else None, shards=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", type=Path, required=True)
    parser.add_argument("--analyzer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frozen-dir", type=Path)
    args = parser.parse_args()
    result = evaluate(args.runs, args.analyzer.resolve(), load_model(args.frozen_dir))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
