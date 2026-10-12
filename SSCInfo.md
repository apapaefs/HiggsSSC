# Historical SSC detector inputs

The numerical baseline is the 1993 GEM Higgs study, *Higgs Searches with the
GEM Detector* (GEM TN-93-373, CALT-68-1856), which updated the GEM TDR using
detailed GEANT detector and object-ID studies. SDC is a comparator, not the
source of the numbers below. These are projected capabilities of unbuilt
detectors, not measured performance.

The executable approximations are documented separately for
[diphotons](hgammagamma/SSC_DETECTOR_RESPONSE.md) and
[four leptons](hfourlepton/README.md). Historical projections should not be
confused with outputs of this repository's simplified simulation.

## Detector baseline

| Object | Historical acceptance/efficiency | Resolution and fake assumptions |
| --- | --- | --- |
| Photons/electrons | EMC `0.1 < \|eta\| < 3`; high-resolution analyses veto `1.01 < \|eta\| < 1.16`; shower ID about 90%, combined geometric/cut efficiency 80–85% | EM sampling 6% barrel, 8–8.5% endcap, constant 0.4%; detailed parameterization uses 8.5%. Position `4.4 mm/sqrt(E)`. |
| Muons | `0.1 < \|eta\| < 2.5`; 85% geometric × 95% reconstruction | At 500 GeV, momentum resolution 5% at eta=0, 12% at |eta|=2.5. No universal jet-to-muon fake probability was specified. |
| Jets | HCAL to `\|eta\|=5.5` | Central energy resolution `60%/sqrt(E) + 4%` in quadrature; the heavy-Higgs forward-tag study used `200%/sqrt(E) + 6%`. |

Energies are in GeV. GEM EM noise was parameterized in transverse energy:
thermal 0.100 GeV (barrel), 0.175 GeV (endcap), plus pileup 0.120 GeV below
`|eta|=1.4`, increasing by `0.366*(|eta|-1.4)` GeV above it, combined in
quadrature. This noise model does not simulate event overlay or occupancy.

## Diphoton channel

GEM required isolated photons with `pT>20 GeV`, `|eta|<2.5`, shower-shape ID,
HCAL veto, tracker discrimination and a `|cos(theta*)|` selection. Final
accepted signal rates were approximately 18–27% of production × BR across the
studied mass points; post-selection trigger efficiency was 98.8%.

| Higgs mass [GeV] | 80 | 100 | 120 | 140 | 160 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline mass resolution [GeV] | 0.66 | 0.77 | 0.84 | 0.91 | 1.0 |
| Degraded EMC [GeV] | 1.3 | 1.5 | 1.7 | 1.8 | 2.0 |

The degraded case used sampling terms 14%/17% and a 1% constant term.
Quark-jet fake rates after isolation and photon ID ranged from
`4.1e-4` to `3.4e-4`; gluon-jet rates from `1.2e-4` to `0.83e-4` over
80–160 GeV. Electron-to-photon tracker-veto probabilities were 0.15% near
`mZ` and 2% elsewhere, conditional on an isolated EM cluster. Apply the
separate EM shower-ID efficiency, not an additional electron-reconstruction
factor. Jet-to-electron rejection was 123,750 (about `8e-6`).

The historical diphoton significances at 10 fb⁻¹ were 4.6, 7.8, 9.0 and 7.3
sigma at 100, 120, 140 and 150 GeV respectively. These forecasts are not a
125 GeV sensitivity result for this code. Prompt diphotons, gamma+jet,
dijet fakes and Drell–Yan fakes all matter; the default implementation does
not yet include every contribution.

## Four-lepton channel

For `140 < mH < 2mZ`, GEM required four isolated leptons with `pT>10 GeV`,
`|eta|<2.5`, the electron crack veto, `R=0.35` isolation below 5 GeV, and
pair windows 70–100 GeV and 10–100 GeV. Trigger efficiency after cuts was
above 98% for 4e and above 99% for the other channels. Backgrounds included
ZZ/ZZ*, heavy flavour and top production; lowering lepton pT to 5 GeV
increased jet background by a factor of 3–4 for a 33% signal gain.

| Higgs mass [GeV] | 140 | 150 | 160 | 170 | 180 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 4e resolution [GeV] | 1.05 | 1.06 | 1.13 | 1.23 | 1.33 |
| 4mu resolution [GeV] | 1.59 | 1.62 | 1.73 | 1.84 | 2.22 |
| 2e2mu resolution [GeV] | 1.36 | 1.46 | 1.56 | 1.71 | 1.77 |
| Combined significance at 10 fb⁻¹ | 11 | 13 | 8.1 | 5.7 | 10 |

These are GEM Table 22 forecasts. Extrapolation to 125 GeV requires simulation
and new acceptance/background estimates. GEM also studied heavier Higgs masses
and forward tagging; that regime is outside this repository's nominal study.

## Historical versus modern analysis

A historically motivated baseline uses cuts, object-ID/fake models and
mass-peak/sideband reasoning. LEP-era likelihood-ratio and CLs combinations,
neural networks and VBF forward tagging provide later comparators. Modern
XGBoost should be labeled a modern exploratory extension, not an original
GEM technique. The 1999 ATLAS TDRs provide useful independent detector/search
comparisons; their assumptions must not be silently substituted for GEM's.

## References

### SSC-Era Detector And Higgs Studies

- GEM Collaboration, *GEM Technical Design Report*, SSC-GEM-TN-93-262,
  SSCL-SR-1219, April 1993. Used here as the detector-performance baseline.
  Fermilab Library Server:
  <https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-262.pdf>.
- S. Mrenna, S. Shevchenko, X. R. Shi, H. Yamamoto, and R. Y. Zhu, *Higgs
  Searches with the GEM Detector*, GEM TN-93-373, CALT-68-1856, June 20, 1993.
  This update to the GEM TDR is the numerical source for the GEM
  $H\to\gamma\gamma$ and $H\to ZZ^{(\ast)}\to 4\ell$ tables summarized above.
  Fermilab Library Server:
  <https://lss.fnal.gov/archive/other/calt-68-1856.pdf>.
- SDC Collaboration, *Solenoidal Detector Collaboration Technical Design
  Report*, SSC-era detector TDR. Used here as a comparator rather than the
  numerical source summarized above; it is listed in the Fermilab SSC report
  index.

### LHC TDR And Higgs Search Methodology

- ATLAS Collaboration, *ATLAS detector and physics performance: Technical
  Design Report, 1*, CERN-LHCC-99-014, ATLAS-TDR-14, CERN Document Server:
  <https://cds.cern.ch/record/391176>.
- ATLAS Collaboration, *ATLAS detector and physics performance: Technical
  Design Report, 2*, CERN-LHCC-99-015, ATLAS-TDR-15, CERN Document Server:
  <https://cds.cern.ch/record/391177>.
- LEP Higgs Working Group for Higgs Boson Searches and the ALEPH, DELPHI, L3,
  and OPAL Collaborations, "Search for the Standard Model Higgs Boson at LEP,"
  arXiv:hep-ex/0306033: <https://arxiv.org/abs/hep-ex/0306033>.
- L3 Collaboration, "Search for the Standard Model Higgs boson in $e^+e^-$
  collisions at $\sqrt{s}$ up to 202 GeV," arXiv:hep-ex/0012019:
  <https://arxiv.org/abs/hep-ex/0012019>.
- A. L. Read, "Presentation of search results: the $CL_s$ technique,"
  *J. Phys. G* **28** (2002) 2693-2704,
  DOI: `10.1088/0954-3899/28/10/313`.

### Weak-Boson-Fusion Phenomenology

- D. Rainwater and D. Zeppenfeld, "Searching for $H\to\gamma\gamma$ in weak
  boson fusion at the LHC," arXiv:hep-ph/9712271:
  <https://arxiv.org/abs/hep-ph/9712271>.
- D. Rainwater, D. Zeppenfeld, and K. Hagiwara, "Searching for
  $H\to\tau\tau$ in weak boson fusion at the LHC," arXiv:hep-ph/9808468:
  <https://arxiv.org/abs/hep-ph/9808468>.
