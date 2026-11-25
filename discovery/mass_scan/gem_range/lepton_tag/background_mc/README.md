# Simulating the lepton-tagged diphoton continuum at 40 TeV

GEM quotes total selected background densities of 0.13 and 0.09 fb/GeV at 120 and 140 GeV. Linear interpolation and integration over 123.3–126.7 GeV give 0.408 fb. This is the reconstructed diphoton mass window for the fixed 125 GeV comparison. The discovery scan also needs the continuum throughout the tagged range of 80–140 GeV. [GEM Higgs study, sections 5.1–5.2 and Tables 19–20](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=60)

## Definition of the calculation

The proton beam energy is 20 TeV per beam. The background calculation uses leading-order MadGraph5_aMC@NLO matrix elements and the NNPDF4.0 NLO QED set (LHAPDF 335900). Herwig 7.3 supplies the parton shower, decays and hadronization. The production cross section in each LHE file supplies the absolute normalization. Each hard collision contributes detector outcomes whose probabilities include object identification and the 98.8% event trigger. The central value uses the project's isolated-photon response and the GEM electron rejection of 33. The jet-to-photon response includes isolation and photon identification. These detector efficiencies enter the estimate of the fake background. The generated collision samples provide an independent check of the corresponding selected cross sections conditional on those efficiencies. [GEM electron identification and jet response](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=43)

| Generated process | Contribution to GEM Table 19 |
| --- | --- |
| Charged-current ℓνγγ | Wγγ together with Wγ and lepton radiation, including interference at tree level |
| Neutral-current ℓ⁺ℓ⁻γγ | Radiative neutral-current contribution with Z/γ* interference |
| e⁺e⁻γ | One electron reconstructed as a photon, a genuine photon and a distinct tagged electron |
| tt̄, tt̄γ and tt̄γγ | Production photons, photons radiated in top or W decays, and electrons or jets reconstructed as photons |

The genuine electroweak samples include all three charged-lepton flavours. The detector selects a final electron or muon. The mass-window count requires two photons and a tagged lepton with transverse momentum above 20 GeV and absolute pseudorapidity below 2.5. The reconstructed pair must have transverse momentum above 40 GeV. The photon and lepton isolation cones have radii 0.45 and 0.30 with a 5 GeV calorimeter-energy threshold. The reconstruction uses the stable particles available in HwSim and the project’s detector smearing and fiducial cracks. Ordinary QCD jet-to-lepton fakes and a dedicated bb̄γγ sample are outside this calculation. GEM judged the latter negligible after isolation. The existing study quantified the available controls on the former. [Selection and background estimates](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=60) · [Inclusive fake-lepton controls](../background-control.md)

GEM enhanced its quark-initiated Zγ estimate by 20% to approximate gluon fusion. The nominal modern inputs use their generated cross sections. A separate refit applies GEM's 20% enhancement to the two Z samples. All background cross sections here are at leading order in the hard process. The signal normalization and the inclusive diphoton continuum come from the [existing discovery calculation](../README.md).

## Statistical interpretation

The statistical target for the direct count is a combined MC standard error of 0.0408 fb at 125 GeV. The covariance sums contributions from independent hard collisions, including correlations between different possible photon pairs in the same collision. Showers of one hard integration share its integration error. A failed shower counts in the generated denominator with zero selection weight. Real collision counts in the discovery simulation fluctuate according to the Poisson distribution in the tagged category.

The combined likelihood uses the same 125 GeV signal and untagged score as the [published GEM-background calculation](../tagged_scan.py). The tagged Poisson means are computed from the simulated continuum. The fit integrates bin contents across all overlapping trial windows, with 123.3 and 126.7 GeV explicit bin boundaries. It calibrates the global background tail and computes the fraction of signal-injected experiments crossing that threshold. Time assumes a single experiment recording 10 fb⁻¹ per year from the start of design operation.

The individual GEM entries imply 0.0918 fb for the charged-current rows together, 0.1403 fb for neutral-current radiation, 0.0374 fb for the electron fake and 0.1445 fb for the two top rows together. Their sum is 0.4140 fb. The 0.408 fb reference uses the separately rounded total row. The associated study groups the top contributions before comparing them because the allocation between production, decay radiation and misidentification depends on the generator. [GEM Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

The production cards and detailed reproduction instructions cover the [electroweak samples](ew_README.md), [top samples and radiative decays](top_backgrounds.md), and [detector response and normalization](response.md).

## Repeating the fit

Run from the repository root with the [discovery dependencies](../../../../requirements.txt). Set `campaign` to the detector-packet directory and `out` to the output directory. The following commands produce the direct spectrum, sideband cross-check, global discovery probabilities and figures.

```bash
code=discovery/mass_scan/gem_range/lepton_tag/background_mc
mkdir -p "$out"
python3 "$code/response_summary.py" "$campaign" "$out/background.json" --top-dileptonic
python3 "$code/response_sidebands.py" "$out/background.json" "$out/sidebands.json"
python3 "$code/refit.py" "$out/background.json" --output "$out"
```

The spectrum contains bin-integrated rates, their covariance and the cross section, seed and input identity for each independent integration. The fit reports a median and a 90% discovery probability. The MC variations change the tagged continuum by one window normalization error in both the background-only and signal-injected experiments. Shape and misidentification model changes are reported separately. The fit also records a predetermined 125 GeV result and validates the global-tail calculation against independent null experiments.
