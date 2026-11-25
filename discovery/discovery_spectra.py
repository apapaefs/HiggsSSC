"""Expected mass spectra for one 40 TeV SSC experiment at the combined discovery exposure."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, StrMethodFormatter
import numpy as np

from input_data import BACKGROUNDS, SIGNAL
from plots import COLORS, NAMES, save


SIGNAL_COLOR = "#a52a31"
BACKGROUND_COLOR = "#203444"


def _axes(ax):
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=10)
    ax.tick_params(which="major", length=5)
    ax.tick_params(which="minor", length=2.8)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.margins(x=0)


def _header(left, right, luminosity):
    left.set_title("SSC projection", loc="left", weight="bold", fontsize=16, pad=9)
    right.set_title(rf"{luminosity:.2f} fb$^{{-1}}$  (40 TeV)", loc="right", fontsize=13, pad=9)


def _diphoton(data, luminosity, spurious_mu, out):
    edges = np.r_[data["low_gev"], data["high_gev"][-1]]
    signal = luminosity * data[SIGNAL]
    # Use the same nominal background and S+B truth as likelihood.Calibration.
    background = luminosity * data["sideband_fit_fb"] + spurious_mu * signal
    fig, (top, excess) = plt.subplots(2, 1, figsize=(7.7, 6.5), sharex=True,
                                     height_ratios=(2.15, 1), layout="constrained")
    for ax in (top, excess):
        _axes(ax)
    _header(top, top, luminosity)
    top.text(.028, .94, r"$H\rightarrow\gamma\gamma$,  $m_H=125$ GeV",
             transform=top.transAxes, va="top", fontsize=12.5)
    stack = np.zeros_like(signal)
    components = []
    for key, name, color in zip(BACKGROUNDS, NAMES, COLORS):
        above = stack + luminosity * data[key]
        components.append(top.stairs(above, edges, baseline=stack, fill=True, facecolor=color,
                                     linewidth=0, label=name))
        stack = above
    total_line = top.stairs(background + signal, edges, baseline=None, color=SIGNAL_COLOR,
                            lw=2.1, label="Signal + background")
    background_line = top.stairs(background, edges, baseline=None, color=BACKGROUND_COLOR,
                                 ls="--", lw=1.45, label="Calibrated background")
    top.legend(handles=[total_line, background_line, *components], loc="upper right",
               ncols=2, frameon=False, fontsize=9, columnspacing=1.4, handlelength=2)
    top.set(ylabel="Expected events / 1 GeV", ylim=(0, max(stack.max(), (background + signal).max()) * 1.23))
    top.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))

    # The counting spread is shown around a fixed background mean, before nuisance profiling.
    poisson_sd = np.sqrt(background)
    band = excess.stairs(poisson_sd, edges, baseline=-poisson_sd, fill=True, facecolor="#dce2e8",
                         linewidth=0, label=r"Background counting spread ($\pm\sqrt{B}$)")
    excess.axhline(0, color=BACKGROUND_COLOR, lw=.8)
    peak = excess.stairs(signal, edges, baseline=None, color=SIGNAL_COLOR, lw=2.1,
                         label="Expected signal")
    excess.set(ylim=(-1.26 * poisson_sd.max(), max(1.23 * signal.max(), 1.65 * poisson_sd.max())),
               xlim=(edges[0], edges[-1]), xlabel=r"$m_{\gamma\gamma}$ [GeV]",
               ylabel="Excess / 1 GeV")
    excess.legend(handles=[peak, band], loc="upper right", frameon=False, fontsize=9)
    save(fig, out, "discovery-diphoton")


def _fourlepton_panel(ax, data, heavy_keys, luminosity, limits, stage):
    edges = np.r_[data["low_gev"], data["high_gev"][-1]]
    width = np.diff(edges)
    if not np.allclose(width, width[0]):
        raise ValueError("Expected uniform four-lepton mass bins")
    signal = luminosity * data["signal_fb"]
    qq = luminosity * data["qq_continuum_fb"]
    gg = luminosity * data["gg_continuum_fb"]
    heavy = luminosity * sum(data[key] for key in heavy_keys)
    background = qq + gg + heavy
    visible = (edges[:-1] >= limits[0]) & (edges[1:] <= limits[1])
    _axes(ax)
    qq_patch = ax.stairs(qq, edges, fill=True, facecolor=COLORS[0], linewidth=0,
                         label=r"$q\bar{q}\rightarrow4\ell$")
    gg_patch = ax.stairs(qq + gg, edges, baseline=qq, fill=True, facecolor=COLORS[1], linewidth=0,
                         label=r"$gg\rightarrow4\ell$")
    heavy_label = r"$Zb\bar{b}$, $t\bar{t}$"
    if not np.any(heavy[visible]):
        heavy_label += " (0 selected)"
    heavy_patch = ax.stairs(background, edges, baseline=qq + gg, fill=True, facecolor=COLORS[3],
                            linewidth=0, label=heavy_label)
    ax.stairs(background + signal, edges, baseline=background, facecolor="#f5dfdc", linewidth=0, fill=True)
    total = ax.stairs(background + signal, edges, color=SIGNAL_COLOR, lw=2.1, baseline=None,
                      label="Signal + background")
    ax.stairs(background, edges, color=BACKGROUND_COLOR, lw=.85, baseline=None)
    ax.set(xlim=limits, ylim=(0, (background + signal)[visible].max() * 1.43),
           xlabel=r"$m_{4\ell}$ [GeV]", ylabel=f"Expected events / {width[0]:g} GeV")
    ax.text(.035, .95, stage, transform=ax.transAxes, va="top", fontsize=12.5)
    return [total, qq_patch, gg_patch, heavy_patch]


def _fourlepton(fiducial, selected, luminosity, out):
    fig, (broad, narrow) = plt.subplots(1, 2, figsize=(10.1, 4.7), width_ratios=(1.55, 1),
                                       layout="constrained")
    _header(broad, narrow, luminosity)
    broad_keys = ("ttbar_fb", "zbb_ee_fb", "zbb_mumu_fb")
    # The continuum was generated over 90–180 GeV. Its reconstructed tails
    # outside this interval do not model its full rate at those masses.
    handles = _fourlepton_panel(broad, fiducial, broad_keys, luminosity, (90, 180), "Fiducial selection")
    broad.text(.035, .87, r"$H\rightarrow ZZ^*\rightarrow4\ell$,  $m_H=125$ GeV",
               transform=broad.transAxes, va="top", fontsize=10.5)
    broad.legend(handles=handles, loc="upper right", frameon=False, fontsize=9.3, handlelength=2)
    broad.xaxis.set_major_locator(MultipleLocator(15))
    broad.xaxis.set_minor_locator(MultipleLocator(5))
    narrow_keys = ("bkg_ttbar_fb", "bkg_zbb_ee_fb", "bkg_zbb_mumu_fb")
    _fourlepton_panel(narrow, selected, narrow_keys, luminosity, (120, 130), "Discovery selection")
    signal_sum = luminosity * np.sum(selected["signal_fb"])
    total_bkg = sum(np.sum(selected[k]) for k in ("qq_continuum_fb", "gg_continuum_fb", *narrow_keys))
    narrow.text(.035, .87, f"{signal_sum:.2f} signal, {luminosity * total_bkg:.2f} background",
                transform=narrow.transAxes, va="top", fontsize=10.5)
    narrow.xaxis.set_major_locator(MultipleLocator(2))
    narrow.xaxis.set_minor_locator(MultipleLocator(.5))
    save(fig, out, "discovery-fourlepton")


def plot_discovery_spectra(spectrum, fourlepton, luminosity, out, *, spurious_mu, fiducial_fourlepton):
    """Write expected counts at a common exposure using the central likelihood inputs."""
    _diphoton(spectrum, luminosity, spurious_mu, out)
    _fourlepton(fiducial_fourlepton, fourlepton, luminosity, out)
