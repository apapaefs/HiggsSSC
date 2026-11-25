"""Update the tagged GEM cut-flow normalization for a 125 GeV Higgs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

from scipy.special import erf

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from gem_rates import historical_generated_signal, selected_signal

SOURCE_SHA256 = "2e58f441ec25b5fb6f82dd83ff786d6e28337029215fd648a2128d55317465c6"
PRODUCTION_PB = {"wp_h": 3.2705023, "wm_h": 2.3891032, "tt_h": 5.4687732}
BR_H_GAMMAGAMMA = .00227
BR_W_LIGHT = 2 * .1086


def transfer(manifest=None):
    if manifest is not None:
        data = manifest.read_bytes()
        if hashlib.sha256(data).hexdigest() != SOURCE_SHA256:
            raise ValueError("Unexpected 40 TeV production manifest")
        for mode, cross_section in PRODUCTION_PB.items():
            samples = [s for s in json.loads(data)["samples"] if s["name"].startswith("signal_" + mode + "_")]
            if len(samples) != 3 or any(s["normalization"]["production_cross_section_pb"] != cross_section
                                        or s["sqrt_s_gev"] != 40000 or s["higgs_mass_gev"] != 125 for s in samples):
                raise ValueError("Production rate or beam energy mismatch for " + mode)
    acceptance = selected_signal(125.) / historical_generated_signal(125.)
    by_mode = {}
    for mode, pb in PRODUCTION_PB.items():
        tag_probability = 1 - (1 - BR_W_LIGHT)**2 if mode == "tt_h" else BR_W_LIGHT
        generated = 1000 * pb * BR_H_GAMMAGAMMA * tag_probability
        by_mode[mode] = dict(production_pb=pb, w_decay_probability=tag_probability,
                             leptonic_diphoton_fb=generated, selected_fb=generated*acceptance)
    total = sum(m["selected_fb"] for m in by_mode.values())
    return dict(production_order="MC@NLO", collider_energy_tev=40, true_mass_gev=125.,
                production_manifest="/home/shared/HiggsSSC/fourlepton/ssc40-20261118/manifest_ssc.json",
                production_manifest_sha256=SOURCE_SHA256, by_mode=by_mode, br_h_gamma_gamma=BR_H_GAMMAGAMMA,
                br_w_to_direct_e_or_mu=BR_W_LIGHT,
                branching_source="https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S043.7",
                historical_conditional_acceptance=acceptance, historical_selected_fb=selected_signal(125.),
                historical_generated_tagged_fb=historical_generated_signal(125.),
                generated_tagged_fb=sum(m["leptonic_diphoton_fb"] for m in by_mode.values()),
                selected_fb=total, selected_pm2sigma_fb=total*erf(2**.5))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, help="Optionally verify the frozen SSC production manifest on Timur")
    print(json.dumps(transfer(parser.parse_args().manifest), indent=2, allow_nan=False))
