"""Translate the same-event variations into a factorization of the GEM window-rate difference."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text())


def compare():
    gem = read("legacy.json")
    matrix = read("hard_projection.json")
    vintage = read("vintage.json")
    history = gem["window"]["gem_interpolated_fb"]
    names = matrix["observable_order"]
    cov = np.asarray(matrix["covariance_fb2"])
    values = np.array([matrix["combined"][key]["window_fb"] for key in names])
    full = names.index("full")
    pole = names.index("production_z_pt10")
    weighted = names.index("production_z_pt10_ehlq")
    pdf_only = names.index("full_ehlq")
    factor_gg = gem["gem_1993"]["enhancement_on_qq_to_z_gamma"]
    tracker = gem["gem_1993"]["electron_fake_probability_in_tag"]
    tracker_scale = tracker / (1./33.)

    def ratio(i,j):
        gradient = np.zeros(len(values))
        gradient[i] = 1/values[j]
        gradient[j] = -values[i]/values[j]**2
        return {"value": values[i]/values[j], "paired_se": np.sqrt(gradient@cov@gradient)}

    def rate(i,scale=1.):
        return {"fb":scale*values[i], "mc_se_fb":scale*np.sqrt(cov[i,i])}

    # The order is explicit. Differences between adjacent rows are paired.
    cases = [("Published simulation",full,1.), ("Production Zγ, hard pT > 10 GeV",pole,1.),
             ("Also use EHLQ1",weighted,1.), ("Also apply GEM’s 20% increase",weighted,factor_gg),
             ("Also use GEM’s 3% tracker fake",weighted,factor_gg*tracker_scale)]
    series = []
    for key,i,scale in cases:
        item = {"name":key,**rate(i,scale)}
        if series:
            _,prev,prev_scale = cases[len(series)-1]
            grad = np.zeros(len(values))
            grad[prev] -= prev_scale
            grad[i] += scale
            item.update(change_from_previous_fb=grad@values,
                        change_paired_mc_se_fb=np.sqrt(max(0.,grad@cov@grad)))
        series.append(item)
    native = vintage["rates_with_object_probabilities_fb"]["123p3_126p7"]["pair_pt"]
    native_scale = factor_gg * tracker_scale
    residual = series[-1]["fb"]/history
    factors = {
        "full_over_resonant_pt10":ratio(full,pole),
        "nnpdf_over_ehlq_at_resonant_pt10":ratio(pole,weighted),
        "nominal_qq_over_gem_qq_plus_gg_prescription":1/factor_gg,
        "nominal_tracker_over_gem_tracker":1/tracker_scale,
        "historical_like_rate_over_table19":residual,
    }
    joined = (factors["full_over_resonant_pt10"]["value"] * factors["nnpdf_over_ehlq_at_resonant_pt10"]["value"]
              /factor_gg/tracker_scale*residual)
    np.testing.assert_allclose(joined,values[full]/history)
    d = values[full]-history
    pdf_contrast = rate(pdf_only)
    pdf_contrast["fraction_of_excess_removed_if_changed_first"]=(values[full]-values[pdf_only])/d
    last_scale = factor_gg * tracker_scale
    gradient = np.zeros(len(values))
    gradient[full] = (last_scale*values[weighted]-history)/d**2
    gradient[weighted] = -last_scale/d
    return dict(mass_window_gev=[123.3,126.7], direct_gem_fb=history,
                full_over_gem=values[full]/history, total_excess_fb=d, same_events=series,
                ordered_factors=factors, pdf_change_alone=pdf_contrast,
                fraction_of_excess_removed_in_order=(values[full]-series[-1]["fb"])/d,
                fraction_of_excess_removed_paired_mc_se=np.sqrt(gradient@cov@gradient),
                residual_over_gem_mc_se=series[-1]["mc_se_fb"]/history,
                native_vintage_qq_fb=native["mean"],
                native_vintage_with_gem_conventions_fb=native_scale*native["mean"],
                native_vintage_with_gem_conventions_mc_se_fb=native_scale*native["se"],
                tracker_probability_needed_with_full_modern_sample=history/values[full]/33.,
                pre_tracker_rate_fb={"modern":33.*values[full], "gem_implied":history/tracker})


def plot(d, output):
    labels = ["Full eeγ, modern PDFs", r"Production Zγ, hard photon $p_T > 10$ GeV", "Same events with EHLQ1",
              "EHLQ1, GEM gluon and fake prescriptions", "Independent PYTHIA 5.6 + our detector", "GEM Table 19"]
    rows = [d["same_events"][i] for i in (0,1,2,4)]
    rates = [r["fb"] for r in rows]+[d["native_vintage_with_gem_conventions_fb"],d["direct_gem_fb"]]
    errors = [r["mc_se_fb"] for r in rows]+[d["native_vintage_with_gem_conventions_mc_se_fb"],0.]
    fig,ax=plt.subplots(figsize=(8.7,3.7),layout="constrained")
    ax.set_axisbelow(True)
    ax.grid(axis="x",color="#dce4e8",linewidth=.8)
    for y,(x,e) in enumerate(zip(rates,errors)):
        ax.errorbar(x,y,xerr=e,fmt=("s" if y==5 else "o"),ms=6,capsize=3,
                    color=("#9d4b2a" if y==5 else "#486d57" if y==4 else "#176782"))
    ax.set(yticks=range(len(labels)),yticklabels=labels,xlim=(0,.31),
           xlabel=r"Selected $ee\gamma$ background in $123.3 < m_{\gamma\gamma} < 126.7$ GeV [fb]")
    ax.invert_yaxis()
    ax.set_title("Electron fakes in the 40 TeV lepton-tagged diphoton search",loc="left",pad=12)
    ax.spines[["top","right","left"]].set_visible(False)
    ax.tick_params(axis="y",length=0)
    for ext in ("pdf","png"):
        fig.savefig(output/f"ee-background-comparison.{ext}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE, help="Directory for the comparison table and plots")
    output = parser.parse_args().output
    output.mkdir(parents=True, exist_ok=True)
    report=compare()
    (output/"comparison.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    plot(report, output)
