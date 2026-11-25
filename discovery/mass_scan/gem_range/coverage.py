"""Plot the Higgs hypotheses tabulated for each channel by the GEM Collaboration.

The data in benchmarks.csv transcribe the original Technical Design Report,
Table 2-3, printed page 2-12. A blank cell means the table gives no estimate.
https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-262.pdf#page=52

The June 1993 update uses the same masses and revises the 800 GeV results
(Table 1, printed page 3). All plotted benchmarks refer to the original TDR.
https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=12
"""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator


HERE = Path(__file__).resolve().parent
CHANNELS = [
    ("inclusive_diphoton_sigma", r"Inclusive $\gamma\gamma$"),
    ("lepton_tagged_diphoton_sigma", r"Lepton-tagged $\gamma\gamma$"),
    ("four_lepton_sigma", r"Four leptons"),
    ("dilepton_dijet_sigma", r"$\ell^+\ell^- jj$"),
    ("dilepton_dineutrino_sigma", r"$\ell^+\ell^- \nu\bar{\nu}$"),
]


def main():
    with (HERE / "benchmarks.csv").open(newline="") as source:
        rows = list(csv.DictReader(source))

    fig, ax = plt.subplots(figsize=(9, 3.9))
    fig.subplots_adjust(left=0.22, right=0.965, top=0.75, bottom=0.19)
    fig.text(0.045, 0.93, "Higgs mass benchmarks in the GEM Technical Design Report",
             fontsize=15, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xlim(75, 875)
    ax.set_ylim(len(CHANNELS) - 0.6, -0.4)
    ticks = [80, 100, 125, 160, 200, 300, 400, 600, 800]
    ax.set_xticks(ticks, labels=[str(tick) for tick in ticks])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_yticks(range(len(CHANNELS)), labels=[label for _, label in CHANNELS])
    ax.tick_params(axis="y", length=0, pad=9)
    ax.grid(axis="y", color="#dedede", linewidth=0.7)
    ax.spines[["top", "left", "right"]].set_visible(False)
    ax.set_xlabel(r"Higgs mass hypothesis $m_H$ [GeV]", labelpad=6)

    signal_mass = ax.axvline(125, color="#a74217", linestyle="--", linewidth=1.5,
                             label="125 GeV Higgs")
    styles = {
        10: {"marker": "o", "markerfacecolor": "#244c6b", "markeredgecolor": "#244c6b"},
        30: {"marker": "D", "markerfacecolor": "white", "markeredgecolor": "#244c6b"},
    }
    for y, (column, _) in enumerate(CHANNELS):
        for exposure, style in styles.items():
            masses = [int(row["mass_gev"]) for row in rows
                      if row[column] and int(row["reference_luminosity_fb_inverse"]) == exposure]
            ax.plot(masses, [y] * len(masses), linestyle="none", markersize=7,
                    markeredgewidth=1.4, zorder=3, **style)

    legend = [Line2D([], [], linestyle="none", markersize=7, markeredgewidth=1.4,
                     label=rf"{exposure} fb$^{{-1}}$", **style)
              for exposure, style in styles.items()]
    fig.legend(handles=[*legend, signal_mass], loc="center left", bbox_to_anchor=(0.035, 0.83),
               ncol=3, frameon=False, handlelength=1.8, columnspacing=1.6)
    fig.savefig(HERE / "coverage.png", dpi=180)
    fig.savefig(HERE / "coverage.pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


if __name__ == "__main__":
    main()
