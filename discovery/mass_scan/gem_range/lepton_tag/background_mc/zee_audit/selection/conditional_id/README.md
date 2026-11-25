# How often does a fake electron pass photon shower identification?

The direct historical test requires two counts from the same sample of isolated electrons. Among those satisfying the lepton-tagged analysis's photon tracking cut, how many also satisfy its photon calorimeter identification? We examined the GEM Higgs analysis and the separate central-tracker design study, then made the corresponding check on our simulated candidates. [GEM Higgs analysis, §4.4](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=45) · [Central tracker, §2.3.2 and Fig. 2-12](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-408.pdf#page=35)

## What GEM published

The tagged Higgs analysis uses a rejection factor of 33 for electrons at about 96% photon acceptance. Its supporting tracking study used isolated electromagnetic clusters, a full simulation of the central tracker and a fast calorimeter simulation. Table 15 gives electron and photon efficiencies for different cuts on the number of tracker hits. GEM studied the detailed calorimeter likelihood separately. Table 19 gives the final tagged electron-fake background. None of these tables gives the two counts for the same electrons at the tagged working point. [GEM tracker efficiencies](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=46) · [Calorimeter identification](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=18) · [Tagged results](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

The central-tracker design study generated Z decays to electron pairs, selected showers with transverse energy of at least 20 GeV, an electromagnetic energy fraction of at least 95% and isolation energy below 10 GeV in a radius-0.3 cone. It *then* measured electron rejection while varying the tracker hit requirement. That sequence specifies the denominator of its efficiency curve. It does not give the fraction passing the tagged photon shower likelihood among electrons already passing the tagged tracking requirement. The Higgs analysis cites a Wang and Sumarokov internal note as being in preparation. [Central-tracker study](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-408.pdf#page=35) · [Reference 16 in the Higgs study](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=91)

## Test on the simulated electrons

We used the 640,000-collision Herwig sample in the 123.3–126.7 GeV diphoton window. The denominator includes the 3% photon-tracking probability for the fake electron, the isolation of all three candidates, the additional-track vetoes, and the calorimeter selection of the genuine photon and the tag electron. We remove just the fake electron's shower identification, then apply its concentration and hadronic-energy cuts and the assumed 90% identification for a clean shower. There is zero added primary-shower leakage in this comparison. [Detector selection and inputs](../README.md) · [Calculation](probe.py)

| Stage | Selected background [fb] | Fraction of denominator |
| --- | ---: | ---: |
| Before fake-electron shower identification | 0.13930 ± 0.01473 | 1 |
| After cuts on surrounding deposits | 0.13923 ± 0.01473 | 99.95% |
| Including the assumed clean-shower efficiency | 0.12531 ± 0.01325 | 89.96% |

The explicit cuts on nearby deposits change the acceptance by about 0.05 percentage points. Using GEM's rate with the Herwig probe denominator after the surrounding-deposit cuts would require a clean-shower efficiency of 0.0374 / 0.13923 = 26.9%. Our 90% clean-shower efficiency is an input, and tracking is represented by a probability. These simulated events can test the effects of nearby particles but cannot measure a correlation between missing tracker hits and the primary shower's own shape. [Rates and paired errors](probe.json)

There is a separate numerical check with the 960,000 native PYTHIA 5.6 collisions. If all the remaining difference from GEM is attributed to the fake electron's clean-shower efficiency and the tracking probability is held at 3%, the two mass regions independently require:

| Reconstructed diphoton mass | Effective clean-shower efficiency |
| --- | ---: |
| 80–100 GeV | (25.74 ± 0.22)% |
| 123.3–126.7 GeV | (25.65 ± 1.04)% |

The difference is −0.09 ± 1.09 percentage points with the same twelve seeds used in both regions. Uncertainties cover our simulation statistics. These inferred efficiencies describe the rate suppression under the stated assumptions. Establishing its origin would require the joint GEM tracking and shower information. [Mass calibration and sensitivity to primary leakage](../mass_control/README.md)

## Reproduction

From the repository root, with NumPy and SciPy installed, run

```bash
python3 discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit/selection/conditional_id/probe.py
```
