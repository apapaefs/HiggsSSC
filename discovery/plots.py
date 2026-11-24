"""Figures for the fixed-mass SSC analysis. Spectra are in fb per 1 GeV bin."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from input_data import BACKGROUNDS, SIGNAL
from likelihood import EDGES, MASS

INK, BLUE, ORANGE = "#20384b", "#136a88", "#b15027"
COLORS = ("#69a9bd", "#b3dcd2", "#e9a65b", "#b687a7", "#64707b")
NAMES = (r"Prompt $\gamma\gamma$", r"$gg\to\gamma\gamma$ box",
         r"$\gamma+$jet", "Dijet double fake", r"Drell–Yan $e^+e^-$")


def save(fig, out, name):
    for ext in ("png", "pdf"):
        metadata = {"CreationDate": None, "ModDate": None} if ext == "pdf" else None
        fig.savefig(out / f"{name}.{ext}", dpi=165, bbox_inches="tight", metadata=metadata)
    plt.close(fig)


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(color="#91a0aa", labelcolor=INK, labelsize=10)
    ax.spines[["bottom", "left"]].set_color("#91a0aa")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e5e8", linewidth=.7)


def plot_spectrum(data, out):
    fig, (top, low) = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True,
                                   height_ratios=(2.2, 1), layout="constrained")
    for ax in (top, low):
        style(ax)
    fig.suptitle(r"Diphoton mass spectrum at the 40 TeV SSC", x=.02, ha="left",
                 weight="bold", color=INK, fontsize=18)
    bottom = np.zeros(len(MASS))
    for sample, color in zip(BACKGROUNDS, COLORS):
        heights = bottom + data[sample]
        top.stairs(heights, EDGES, baseline=bottom, fill=True, color=color, linewidth=0)
        bottom = heights
    top.stairs(data["sideband_fit_fb"], EDGES, color=INK, linewidth=1.7)
    top.errorbar(MASS, data["background_fb"], yerr=data["background_mc_sigma_fb"],
                 ls="none", color=INK, marker=".", markersize=2.5, elinewidth=.7)
    handles = [Patch(facecolor=c, label=n) for c, n in zip(COLORS, NAMES)]
    handles += [Line2D([], [], color=INK, label="Sideband fit"),
                Line2D([], [], color=INK, marker=".", ls="none", label="Total background MC")]
    top.legend(handles=handles, frameon=False, ncol=2, loc="upper right", fontsize=9)
    top.set_ylabel("Cross section [fb / GeV]")
    top.set_ylim(0, max(bottom) * 1.14)
    for ax in (top, low):
        for edge in (120, 130):
            ax.axvline(edge, color=INK, lw=.85, ls=(0, (3, 3)), alpha=.65)
    low.axvspan(120, 130, color="#e5e9ec", alpha=.7, linewidth=0)
    low.axhline(0, color="#89969e", lw=.9)
    low.errorbar(MASS, data["residual_fb"], yerr=data["residual_mc_sigma_fb"],
                 color=INK, ms=3, fmt="o", elinewidth=.8, label="Background MC minus sideband fit")
    low.stairs(data[SIGNAL], EDGES, color=ORANGE, lw=2, label="SM Higgs signal")
    low.set_ylabel("Excess [fb / GeV]")
    low.set_xlabel(r"Reconstructed $m_{\gamma\gamma}$ [GeV]")
    low.set_xlim(105, 160)
    low.legend(loc="upper right", frameon=False, fontsize=8.7)
    save(fig, out, "diphoton-spectrum")


def plot_reach(rows, gem, out):
    fig, ax = plt.subplots(figsize=(9.2, 4.3), layout="constrained")
    values = {(r["channel"], r["heavy_flavour"]): r["luminosity_fb_inverse"] for r in rows
              if r["model"] in ("exponential2", "exact_poisson") and r["treatment"] != "smooth"}
    labels = [r"$\gamma\gamma+4\ell$", r"$\gamma\gamma$", r"$4\ell$"]
    for y, channel in enumerate(("combined", "diphoton", "fourlepton")):
        if channel == "diphoton":
            central = values[channel, "both"]
            ax.scatter(central, y, s=70, color=BLUE, zorder=3)
            ax.text(central + .20, y, f"{central:.2f}", va="center", weight="bold", color=INK)
            continue
        central, upper = values[channel, "central"], values[channel, "upper"]
        ax.plot([central, upper], [y, y], color="#c6cfd4", lw=2, zorder=1)
        ax.scatter(central, y, s=70, color=BLUE, zorder=3)
        ax.scatter(upper, y, s=74, marker="D", color=ORANGE, zorder=3)
        ax.text(central - .17, y, f"{central:.2f}", ha="right", va="center", weight="bold", color=INK)
        ax.text(upper + .20, y, f"{upper:.2f}", va="center", weight="bold", color=INK)
    ax.set_yticks(range(3), labels, fontsize=13)
    ax.set_ylim(-.45, 2.5)
    ax.set_xlim(0, 11.6)
    ax.set_xlabel(r"Integrated luminosity [fb$^{-1}$]")
    ax.secondary_xaxis("top", functions=(lambda l: 1.2*l, lambda m: m/1.2)).set_xlabel(
        r"Design months at 10 fb$^{-1}$ per year")
    ax.grid(axis="x", color="#e0e6ea", linewidth=.8)
    ax.spines[["left", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.suptitle(r"Median local $5\sigma$ discovery at $m_H=125$ GeV", x=.02,
                 ha="left", fontsize=17, color=INK, weight="bold")
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", color=BLUE,
                              label="Central background"),
                       Line2D([], [], ls="none", marker="D", color=ORANGE,
                              label=r"$4\ell$ heavy-flavour upper bound")],
              loc="lower right", frameon=False, fontsize=9)
    save(fig, out, "discovery-reach")

    fig, ax = plt.subplots(figsize=(7.8, 3.8), layout="constrained")
    style(ax)
    points = [(r["mass_gev"], r["estimated_five_sigma_fb_inverse"]) for r in gem]
    ax.scatter(*zip(*points), marker="s", s=60, color=ORANGE, label="GEM Table 18, scaled from 10 fb⁻¹")
    modern = values["diphoton", "both"]
    ax.scatter(125, modern, s=75, color=BLUE, label="This analysis, calibrated sideband fit")
    for mass, lum in [*points, (125, modern)]:
        ax.annotate(f"{lum:.2f} fb⁻¹", (mass, lum), xytext=(7, 4), textcoords="offset points", color=INK)
    ax.set(xticks=(120, 125, 130, 135, 140), xlim=(117, 144), ylim=(2.7, 4.55),
           xlabel="Higgs mass [GeV]", ylabel=r"Luminosity for $5\sigma$ [fb$^{-1}$]")
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    fig.suptitle("Inclusive diphoton discovery projections", fontsize=16, x=.02, ha="left", color=INK, weight="bold")
    save(fig, out, "gem-diphoton")
