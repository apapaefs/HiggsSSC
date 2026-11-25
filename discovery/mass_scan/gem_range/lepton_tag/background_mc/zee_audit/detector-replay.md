# Replaying the electron-fake cutflow

The calculation uses the 352,000 core and 288,000 complementary eeγ hard events from the final tagged-background fit. Each shower uses its original hard-event interval, Herwig card, plugin and random seed. The counter reuses the production particle identification, smearing, isolation, electron ancestry and two-leading-candidate probability functions. It checks the saved event count and compares every accepted 75–165 GeV pair, by event, mass and probability, against its original detector packet. The electron supplying the lepton tag must differ from the electron supplying the fake photon. A true photon must have the generated hard photon in its ancestry.

The results give successive cut rates for the reconstructed 123.3–126.7 and 75–165 GeV mass intervals. Each pre-identification row counts kinematically possible assignments of a genuine photon and one electron treated as a photon. After identification, the observable is an event rate with exclusive detector outcomes. The full covariance propagates repeated assignments in the same collision, each source's independent hard integration error and the failed showers in the generated denominator. The two disjoint generator regions add in cross section. The JSON includes their separate cutflows and the combined ones. [Detector cutflow and covariance matrices](detector-cutflow.json)

The nominal cutflow uses the original SSC calorimeter smearing, dead regions and fixed 5 GeV isolation thresholds. The unsmeared variation removes energy smearing. A further selection variation uses unsmeared energies, full angular coverage inside |η| < 2.5, and an isolation threshold of 5 GeV + 0.10 times the candidate transverse momentum, with cone radii 0.45 for photons and 0.30 for the tagged electron. This last variation approximates the acceptance model in the 1992 predecessor of the GEM study. That publication describes its uncorrected ideal calorimeter, an angular shower-shape cut and a 3% electron fake probability. The variation uses the same modern shower and detector identification probabilities for the before/after pair-transverse-momentum ratio. A direct absolute prediction for the 1992 apparatus would also require its shower-shape response and tagging implementation. [GEM TN-92-126, Sec. 7.1 and Table 15](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf#page=33)

## Measured selection at 125 GeV

These are cumulative rates for the 3.4 GeV mass window. The conditional efficiencies follow the order in the table. Plateau ID probabilities are applied after all kinematic cuts so their contributions can be read directly.

| Selection | Rate [fb] |
| --- | ---: |
| Reconstructed eγ mass in 123.3–126.7 GeV, distinct truth-matched electron available to tag | 172.30 ± 5.79 |
| Genuine photon pT ≥ 20 GeV and absolute η < 2.5 | 98.76 ± 3.84 |
| Electron supplying fake photon, same kinematics | 73.39 ± 3.00 |
| Tagged electron, same kinematics | 52.79 ± 2.31 |
| Fiducial cracks and central electromagnetic calorimeter gap for all three objects | 37.77 ± 2.02 |
| Genuine photon isolation | 27.86 ± 1.66 |
| Fake-electron isolation in the photon cone | 21.30 ± 1.49 |
| Tag-electron isolation | 19.81 ± 1.45 |
| Reconstructed pair pT > 40 GeV | 12.287 ± 1.248 |
| Genuine photon ID, 0.90 × 0.96 | 10.616 ± 1.079 |
| Tag ID, 0.90 | 9.555 ± 0.971 |
| Electron-to-photon transfer, 0.90/33 | 0.2606 ± 0.0265 |
| Event trigger, 0.988 | **0.25746 ± 0.02616** |

The fake electron has charge −1 in 0.14561 ± 0.01763 fb and +1 in 0.11185 ± 0.01708 fb. Both charges give allowed reconstructions in a single Higgs-window collision. The event-level normalization includes mutually exclusive ID assignments. The final rate agrees with every published eligible event and its individual mass and probability, in the 352,000 nominal and 288,000 wider-acceptance collisions. The wider region contributes 0.04476 ± 0.01918 fb. All 4,387 collisions contributing to the 75–165 GeV comparison pass this independent replay check.

| Reconstructed mass and selection | Before pair pT cut [fb] | After cut [fb] | Passing fraction |
| --- | ---: | ---: | ---: |
| 123.3–126.7 GeV, nominal response | 0.4151 ± 0.0304 | 0.2575 ± 0.0262 | 62.0 ± 3.3% |
| 75–165 GeV, nominal response | 15.380 ± 0.191 | 8.124 ± 0.141 | 52.82 ± 0.55% |
| 75–165 GeV, nominal geometry and isolation, unsmeared energies | 15.357 ± 0.190 | 8.095 ± 0.140 | 52.71 ± 0.55% |
| 75–165 GeV, ideal geometry, relative isolation and unsmeared energies, same ID probabilities | 28.368 ± 0.258 | 14.886 ± 0.187 | 52.47 ± 0.40% |
| 75–165 GeV, GEM 1992 Table 15, photon cone 0.45 | 5.1 | 3.2 | 62.75% |

The precursor gave 51 and 32 expected events for 10 fb⁻¹. The ratio of the full modern rate to that early GEM result decreases from 3.02 before the pair pT requirement to 2.54 after it. In the same modern simulation the earlier acceptance and isolation prescription approximately doubles the absolute broad-window rate while giving a similar pT passing fraction. The 1992 study used its own detector working points, and the 1993 result at 125 GeV has a different mass interval. [1992 GEM Table 15](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf#page=34) · [1993 GEM Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

## Reproduction on Timur

The raw student release includes the input hard events, original shower receipts and archived plugin binaries. With the discovery dependencies, ROOT and Herwig 7.3 available, run from the repository root. The `work` path should have space for the shower trees, approximately 25 GB for this run.

```bash
release=/home/shared/HiggsSSC/student-results/ssc40-tag-background-20261125
code=discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit
work=$PWD/zee-detector-replay
mkdir -p "$work/core" "$work/wide" "$work/showers"
tar -xzf "$release/raw/core-inputs.tar.gz" -C "$work/core"
tar -xzf "$release/raw/zee-confirmation-inputs.tar.gz" -C "$work/wide"
g++ -O2 -I"$PWD" $(root-config --cflags) "$code/detector_cutflow.cc" \
    $(root-config --libs) -o "$work/detector_cutflow"
python3 "$code/replay_showers.py" "$release/background.json.gz" "$work/core" "$work/core/ew" \
    "$release/zee-confirmation-spectrum.json.gz" "$work/wide" "$work/wide" "$work/showers" --workers 50
python3 "$code/detector_cutflow.py" "$release/background.json.gz" "$work/core" "$work/core/ew" \
    "$release/zee-confirmation-spectrum.json.gz" "$work/wide" "$work/wide" "$work/showers" \
    "$work/detector_cutflow" "$work/detector-cutflow.json" --workers 12
```
