"""Selected four-lepton observables with signed, source-grouped uncertainties."""
from __future__ import annotations
from collections import defaultdict
import csv
from dataclasses import dataclass
import html
import json
import math
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Observable:
    name: str
    label: str
    branch: str
    low: float
    high: float
    bins: int
    unit: str = ""
    index: int | None = None
    absolute: bool = False

    @property
    def edges(self):
        return np.linspace(self.low, self.high, self.bins + 1)

    def value(self, row):
        value = row[self.branch]
        if self.index is not None:
            value = value[self.index]
        value = abs(float(value)) if self.absolute else float(value)
        if not math.isfinite(value):
            raise ValueError(f"nonfinite {self.name}")
        return value


OBSERVABLES = (
    Observable("m4l", r"$m_{4\ell}$ [GeV]", "m4l", 70., 200., 65, "GeV"),
    Observable("mZ1", r"$m_{Z_1}$ [GeV]", "mZ1", 40., 120., 40, "GeV"),
    Observable("mZ2", r"$m_{Z_2}$ [GeV]", "mZ2", 0., 120., 40, "GeV"),
    Observable("pt4l", r"$p_T^{4\ell}$ [GeV]", "pt4l", 0., 300., 30, "GeV"),
    Observable("y4l", r"$y_{4\ell}$", "y4l", -5., 5., 40),
    Observable("abs_y4l", r"$|y_{4\ell}|$", "y4l", 0., 5., 25, absolute=True),
    *(Observable(f"pt{i+1}", rf"$p_T^{{\ell_{i+1}}}$ [GeV]", "lepton_pt", 0., high, 40, "GeV", i)
      for i, high in enumerate((200., 160., 100., 80.))),
    *(Observable(f"eta{i+1}", rf"$\eta_{{\ell_{i+1}}}$", "lepton_eta", -3., 3., 30, index=i)
      for i in range(4)),
    Observable("min_delta_r", r"$\min\,\Delta R_{\ell\ell}$", "min_delta_r", 0., 4., 40),
    Observable("z1_delta_r", r"$\Delta R_{\ell\ell}(Z_1)$", "z1_delta_r", 0., 6., 40),
    Observable("z2_delta_r", r"$\Delta R_{\ell\ell}(Z_2)$", "z2_delta_r", 0., 6., 40),
)


def histogram(rows, observable, scale_pb):
    """Sum hypotheses before squaring; include covariance with the shape total.

    Flow bins are explicit and excluded only from unit-area normalization.
    A source contributing to two bins correlates both with the same total.
    """
    if scale_pb is None or not math.isfinite(scale_pb):
        return {"status": "normalization_unavailable"}
    edges = observable.edges
    grouped = defaultdict(float)
    try:
        for row in rows:
            value = observable.value(row)
            weight = float(row["event_weight"]) * scale_pb
            if not math.isfinite(weight):
                raise ValueError("nonfinite physical event weight")
            index = int(np.searchsorted(edges, value, side="right"))
            grouped[(int(row["source_index"]), index)] += weight
    except (KeyError, IndexError) as error:
        return {"status": "branch_unavailable", "branch": observable.branch, "detail": str(error)}
    total_by_source = defaultdict(float)
    for (source, index), weight in grouped.items():
        if 0 < index < len(edges):
            total_by_source[source] += weight
    sums = np.zeros(len(edges) + 1)
    variances = np.zeros_like(sums)
    covariance = np.zeros_like(sums)
    entries = np.zeros_like(sums, dtype=int)
    for (source, index), weight in grouped.items():
        sums[index] += weight
        variances[index] += weight * weight
        covariance[index] += weight * total_by_source[source]
        entries[index] += 1
    return dict(status="ok", edges=edges.tolist(), cross_section_pb=sums.tolist(),
        variance_pb2=variances.tolist(), covariance_with_regular_total_pb2=covariance.tolist(),
        source_events=entries.tolist(), regular_total_pb=float(np.sum(sums[1:-1])),
        regular_total_variance_pb2=sum(w * w for w in total_by_source.values()))


def sample_observables(rows, normalization, signal_region):
    from hfourlepton.make_fourlepton_report import normalize_channel
    result = {}
    low, high = signal_region
    for region in ("selected", "signal_region"):
        selected = rows if region == "selected" else [r for r in rows if low < r["m4l"] < high]
        channels = {"inclusive": selected}
        for channel in ("4e", "4mu", "2e2mu"):
            channels[channel] = [r for r in selected if normalize_channel(r.get("channel")) == channel]
        result[region] = {channel: {spec.name: histogram(values, spec,
            normalization.get("scale_pb_per_weight")) for spec in OBSERVABLES}
            for channel, values in channels.items()}
    return result


def combine(histograms):
    if not histograms or any(h["status"] != "ok" for h in histograms):
        return {"status": "component_unavailable"}
    first = histograms[0]
    if any(h["edges"] != first["edges"] for h in histograms):
        raise ValueError("observable bin edges differ")
    result = dict(status="ok", edges=first["edges"])
    for name in ("cross_section_pb", "variance_pb2", "covariance_with_regular_total_pb2", "source_events"):
        result[name] = np.sum([h[name] for h in histograms], axis=0).tolist()
    for name in ("regular_total_pb", "regular_total_variance_pb2"):
        result[name] = sum(h[name] for h in histograms)
    return result


def density(hist):
    if hist["status"] != "ok" or hist["regular_total_pb"] <= 0.:
        return None
    values = np.asarray(hist["cross_section_pb"])[1:-1]
    variance = np.asarray(hist["variance_pb2"])[1:-1]
    covariance = np.asarray(hist["covariance_with_regular_total_pb2"])[1:-1]
    total = hist["regular_total_pb"]
    widths = np.diff(hist["edges"])
    errors2 = (variance / total**2 + values**2 * hist["regular_total_variance_pb2"] / total**4
               - 2. * values * covariance / total**3) / widths**2
    return values / total / widths, np.sqrt(np.maximum(0., errors2))


def component(sample):
    if sample["category"] == "signal":
        for prefix, label in (("signal_gg_h_", "ggF"), ("signal_vbf_", "VBF"),
                              ("signal_wp_", "WH"), ("signal_wm_", "WH"),
                              ("signal_z_", "ZH"), ("signal_gg_z_", "ZH"), ("signal_tt_", "ttH")):
            if sample["name"].startswith(prefix):
                return label
        return "Other signal"
    return "qq continuum" if sample["name"].startswith("bkg_qq_") else (
        "gg continuum" if sample["name"].startswith("bkg_gg_") else "Other background")


COLORS = {"qq continuum": "#d59527", "gg continuum": "#9f5b22", "ggF": "#4385ba",
          "VBF": "#65a7cf", "WH": "#83bd9a", "ZH": "#677eb1", "ttH": "#ac87b3",
          "Other background": "#999999", "Other signal": "#8888dd", "Signal": "#205783"}


def render(report, output):
    """Write paired plots plus per-sample/channel data, without rereading ROOT."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output.mkdir(parents=True, exist_ok=True)
    (output / "plots").mkdir(exist_ok=True)
    (output / "data").mkdir(exist_ok=True)
    if report["audit"]["errors"]:
        (output / "index.html").write_text("<h1>Observables unavailable</h1><p>Resolve the report audit errors first.</p>")
        return {"status": "audit_failed", "plots": 0}
    samples = [s for s in report["samples"] if s["included_in_physics_totals"]]
    profiles = sorted({s["detector_profile"] for s in samples})
    luminosity = report["metadata"]["luminosity_fb"]
    energy = report["metadata"]["sqrt_s_tev"]
    sr = report["metadata"]["signal_region_gev"]
    plots = []
    data = {"luminosity_fb": luminosity, "sqrt_s_tev": energy, "signal_region_gev": sr,
            "uncertainties": "MC statistical only; signed source weights; correlated detector hypotheses",
            "samples": {f"{s['detector_profile']}/{s['name']}": s["observables"] for s in samples}}
    (output / "histograms.json").write_text(json.dumps(data, allow_nan=False) + "\n")
    for spec in OBSERVABLES:
        with (output / "data" / f"{spec.name}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["profile", "sample", "component", "region", "channel", "bin_type", "low", "high",
                             "cross_section_pb", "mc_error_pb", "events", "mc_error_events", "source_events"])
            for sample in samples:
                for region, channels in sample["observables"].items():
                    for channel, values in channels.items():
                        h = values[spec.name]
                        if h["status"] != "ok":
                            continue
                        for i, (value, variance, n) in enumerate(zip(h["cross_section_pb"], h["variance_pb2"], h["source_events"])):
                            flow = "underflow" if i == 0 else "overflow" if i == len(h["edges"]) else "regular"
                            writer.writerow([sample["detector_profile"], sample["name"], component(sample), region, channel,
                                flow, h["edges"][i-1] if i else "", h["edges"][i] if i < len(h["edges"]) else "",
                                value, math.sqrt(variance), value * luminosity * 1000., math.sqrt(variance) * luminosity * 1000., n])
        for profile in profiles:
            profile_samples = [s for s in samples if s["detector_profile"] == profile]
            for region in ("selected", "signal_region"):
                groups = defaultdict(list)
                signals = []
                for sample in profile_samples:
                    h = sample["observables"][region]["inclusive"][spec.name]
                    groups[component(sample)].append(h)
                    if sample["category"] == "signal":
                        signals.append(h)
                combined = {key: combine(rows) for key, rows in groups.items()}
                title = (f"{energy:g} TeV" if energy is not None else "Collider energy unspecified") + f" · {profile} response · " + (
                    "fiducial selection" if region == "selected" else f"{sr[0]:g} < m4l < {sr[1]:g} GeV")
                fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
                shape_groups = {"Signal": combine(signals), **{k: v for k, v in combined.items() if "continuum" in k or k == "Other background"}}
                unavailable = []
                for label, h in shape_groups.items():
                    values = density(h)
                    if values is None:
                        unavailable.append(label)
                        continue
                    y, error = values
                    axes[0].stairs(y, spec.edges, label=label, color=COLORS[label], linewidth=1.5)
                    axes[0].fill_between(spec.edges, np.r_[y-error, (y-error)[-1]],
                        np.r_[y+error, (y+error)[-1]], step="post", color=COLORS[label], alpha=.13)
                if unavailable:
                    axes[0].text(.02, .02, "Shape unavailable: " + ", ".join(unavailable),
                                 transform=axes[0].transAxes, fontsize=7)
                base = np.zeros(spec.bins)
                total = combine(list(combined.values()))
                ordered = [key for key in COLORS if key in combined] if total["status"] == "ok" else []
                for label in ordered:
                    y = np.asarray(combined[label]["cross_section_pb"])[1:-1] * luminosity * 1000.
                    axes[1].stairs(base+y, spec.edges, baseline=base, fill=True, color=COLORS[label], alpha=.85, label=label)
                    base += y
                if total["status"] == "ok":
                    error = np.sqrt(total["variance_pb2"])[1:-1] * luminosity * 1000.
                    axes[1].stairs(base, spec.edges, color="black", linewidth=.8)
                    axes[1].fill_between(spec.edges, np.r_[base-error, (base-error)[-1]],
                        np.r_[base+error, (base+error)[-1]], step="post", color="black", alpha=.17, label="MC stat.")
                else:
                    axes[1].text(.05, .5, "Yield plot unavailable: a component\nlacks the observable or normalization.",
                                 transform=axes[1].transAxes, fontsize=9)
                for ax in axes:
                    ax.set_xlabel(spec.label)
                    ax.set_xlim(spec.low, spec.high)
                    ax.tick_params(which="both", direction="in", top=True, right=True)
                    ax.axhline(0., color="black", linewidth=.5)
                    ax.legend(frameon=False, fontsize=8)
                    ax.margins(y=.2)
                axes[0].set_ylabel("Unit-area density" + (f" [1/{spec.unit}]" if spec.unit else ""))
                axes[1].set_ylabel(rf"Events / bin at {luminosity:g} fb$^{{-1}}$")
                fig.suptitle(title, fontsize=11)
                fig.tight_layout()
                name = f"{profile}_{region}_{spec.name}"
                for extension in ("png", "svg", "pdf"):
                    fig.savefig(output / "plots" / f"{name}.{extension}", dpi=170)
                plt.close(fig)
                plots.append(dict(name=name, observable=spec.name, profile=profile, region=region, title=title))
    cards = "\n".join(f'<article data-profile="{p["profile"]}" data-region="{p["region"]}"><h2>{html.escape(p["observable"])}</h2>'
        f'<img src="plots/{p["name"]}.png" alt="{html.escape(p["title"])} {p["observable"]}">'
        f'<p><a href="plots/{p["name"]}.pdf">PDF</a> · <a href="plots/{p["name"]}.svg">SVG</a> · '
        f'<a href="data/{p["observable"]}.csv">CSV (all channels)</a></p></article>' for p in plots)
    options = "".join(f'<option>{p}</option>' for p in profiles)
    page = '''<!doctype html><html><meta charset="utf-8"><title>Four-lepton observables</title>
<style>body{font:16px system-ui;color:#172b40;background:#f4f6f9;max-width:1200px;margin:35px auto;padding:0 20px}
article{background:white;padding:20px;margin:22px 0;border-radius:9px}img{width:100%}a{color:#196a96}select{padding:8px;margin:8px}</style>
<h1>Four-lepton observables</h1><p><a href="../index.html">Analysis audit and cutflows</a> · <a href="histograms.json">All histogram data</a></p>
<p>Unit-area shapes and physical event yields use signed weights, with no extra signal scaling.
Bands show Monte Carlo statistical uncertainties, grouping correlated detector hypotheses by source event.
Shape errors include their normalization covariance. Scale, PDF and detector systematics are not included.</p>
<p>Leptons are ordered by transverse momentum; Z1 is the opposite-sign, same-flavour pair closest to the Z mass.
The first view uses the final fiducial selection before a Higgs-mass cut; the signal-region view adds the mass window.
The available continuum uses the campaign's hard-generation mass window. Underflow and overflow are exported;
unit-area shapes are normalized over the displayed range. Negative bins are not clipped. All three decay channels
are combined in the figures and supplied separately in CSV/JSON.</p>
<label>Response <select id="profile">OPTIONS</select></label>
<label>Selection <select id="region"><option value="selected">Fiducial selection</option><option value="signal_region">Higgs mass window</option></select></label>
CARDS<script>function update(){document.querySelectorAll('article').forEach(x=>x.hidden=x.dataset.profile!==document.querySelector('#profile').value||x.dataset.region!==document.querySelector('#region').value)}
document.querySelectorAll('select').forEach(x=>x.onchange=update);update();</script></html>'''
    (output / "index.html").write_text(page.replace("OPTIONS", options).replace("CARDS", cards))
    return {"status": "complete", "plots": len(plots), "variables": [s.name for s in OBSERVABLES],
            "path": "observables/index.html"}
