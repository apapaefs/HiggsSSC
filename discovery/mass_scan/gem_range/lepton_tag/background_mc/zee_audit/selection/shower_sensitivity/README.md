# Radiation and the tagged electron-fake background at 40 TeV

We varied one generator setting at a time in native PYTHIA 5.6 with JETSET 7.3. Every run uses the same Born process, EHLQ1 hard PDFs, electron-only Z decay, 10 GeV hard transverse-momentum cut, GEM 1.20 gluon prescription and 3% tracker fake probability. All cases use the same lepton-tag, photon selection, 40 GeV diphoton transverse-momentum cut and detector response with zero added leakage from primary showers. The reconstructed diphoton mass window is 123.3–126.7 GeV. [Detector model](../README.md) · [GEM selection and background definitions](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=60)

Each setting has twelve independent seeds of 80,000 eeγ collisions. The nominal run repeats and checks all twelve previously published selected rates seed by seed. The other five settings use disjoint seeds. The quoted errors describe simulation statistics from the dispersion across seeds. The common nominal denominator correlates the five rate ratios. GEM published no component uncertainty on its interpolated 0.0374 fb.

| Generator setting | Pre-isolation rate [fb] | Selected rate [fb] | Selected / baseline | Mean bare eeγ recoil [GeV] |
| --- | ---: | ---: | ---: | ---: |
| GEM-era baseline | 0.1761 ± 0.0066 | 0.1312 ± 0.0054 | 1 | 36.8 ± 1.0 |
| Half the ISR scale | 0.1877 ± 0.0059 | 0.1386 ± 0.0045 | 1.056 ± 0.055 | 44.6 ± 1.0 |
| Double the ISR scale | 0.1705 ± 0.0034 | 0.1288 ± 0.0027 | 0.981 ± 0.045 | 36.5 ± 0.8 |
| Initial-state radiation disabled | 0.1164 ± 0.0036 | 0.0961 ± 0.0032 | 0.732 ± 0.039 | 1.4 ± 0.2 |
| Time-like shower photon radiation disabled | 0.1803 ± 0.0053 | 0.1393 ± 0.0038 | 1.062 ± 0.052 | 37.8 ± 0.7 |
| Additional parton interactions disabled | 0.1723 ± 0.0075 | 0.1314 ± 0.0058 | 1.001 ± 0.060 | 36.7 ± 1.0 |

Halving and doubling the ISR scale give shifts of **+5.6 ± 5.5%** and **−1.9 ± 4.5%** in the selected rate. Completely removing initial-state radiation lowers the rate by **26.8 ± 3.9%**, leaving **2.57 times GEM**. The hard cross sections for the six runs span 5.979–5.986 pb before the common 1.20 weight. Removing ISR chiefly reduces the pre-isolation population, while the fraction passing isolation and the remaining cluster and track cuts increases from 74.5% to 82.6%. No setting in this set supplies the approximately 72% reduction needed to reach GEM.

The ISR scale changes set `PARP(64)` to 0.25 or 4 with the default `MSTP(64)=2`. This halves or doubles the transverse-momentum scale entering the model's initial-state evolution. Disabling ISR, disabling time-like shower photon radiation off quarks and leptons, and removing extra parton interactions are separate diagnostic extremes. Hard-process photon production stays active in every case. As a check on the response to ISR, we measure the transverse momentum of the bare three-object system using the final hard-process photon and final Z-decay electrons. This quantity can also respond to final-state photon radiation and intrinsic parton transverse momentum. A genuine photon and a misidentified electron can already form a pair with transverse momentum above 40 GeV at Born level by recoiling against the tag electron. [Original PYTHIA and JETSET manual](https://cds.cern.ch/record/2296393/files/pythia.pdf)

The pre-isolation column includes the same kinematic, geometry, mass and pair momentum constraints and fixed object weights as the final rate. The acceptance between the two columns covers calorimeter isolation and noise, neighboring deposits, and extra tracks. Event-level sums over the two mutually exclusive electron assignments precede the independent seed averages. The mean three-object recoil and the mean thresholded transverse energy from nearby particles in the fake candidate's isolation cone use the pre-isolation population. The energy measurement excludes the primary shower and electronics noise. [Per-seed rates, recoil and cone observables](rates.json)

From the repository root, Python with NumPy and SciPy, C++17, ROOT and gfortran reproduce the calculation. The build fetches the historical sources and checks their hashes. The per-seed file records occasional JETSET fragmentation diagnostics, including 14 messages in the doubled-scale run. All 960,000 generated collisions in each setting had the expected hard photon and electron ancestry.

```bash
audit=discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit
work=$(mktemp -d)
bash "$audit/build_vintage.sh" "$work" --selection
python3 "$audit/selection/shower_sensitivity/compare.py" "$work/vintage-selection" "$work/rates.json" --workers 8
```
