# Validation record

A successful pilot establishes executable interfaces, event accounting and
output generation. It does not establish precision, detector realism,
matching/isolation stability or discovery reach. Preserve raw logs, manifests,
input cards and generated products alongside any quoted numerical result.

## Common workflow audit — October 2026

Fresh validation uses an isolated instructor-owned Timur checkout and a new
private installation; no generated LHE, grids or normalization records were
borrowed. The shared Herwig/HwSim runtime is used read-only.

- LO setup installed the model source and all 101 LO PDF members. HO setup
  subsequently added its dependencies to the same installation.
- The final-source `lo_pilot_02` completed all four samples with 30 LHE,
  30 HwSim and 30 analyzed source events each. All source weights were positive.
  Response-row counts were 81 (signal), 40 (prompt diphoton), 38 (gamma+jet)
  and 63 (Drell–Yan); rows represent weighted detector outcomes, not extra
  independent events. Tree/summary weight closure was within `1.5e-14`.
- All 41 manifest artifact hashes/sizes and 16 pilot hashes matched. Repeating
  the command reused all 12 generation/shower/analysis stages. Report-only
  regeneration at 200 fb⁻¹ preserved pilot evidence. The report contained
  12 PNGs, 12 SVGs and 13 CSVs; a representative plot was inspected.
- Baseline cuts at 100 fb⁻¹ used all four samples. No background hypotheses
  passed the mass window, so significance remained undefined; this pilot
  cannot establish background rejection or sensitivity.
- 441 local regression tests passed, with seven optional ML tests skipped.
  Documentation links, shell examples and dependency-free CLI help were checked.

Fresh full POWHEG, all 109 ihixs points, the four-sample HO pilot and both
four-lepton responses remain in progress. Failed development attempts are
preserved. Older checks below do not certify these unfinished acceptance tests.

## Independent HO entry point — October 12, 2026

The preceding HO-only entry point was exercised in an isolated Timur checkout:

- Fresh private PDF installation, Cuba 4.2.2, MG5 3.5.15 configuration, POWHEG
  compilation, analysis compilation and Herwig/HwSim loading completed.
- The fresh 13 TeV, `mH=125 GeV`, pure-HEFT upstream benchmark with
  `PDF4LHC15_nnlo_100` returned `45.1815532 +/- 0.0109436459 pb`, compatible
  with the 45.1816 pb reference. This is **not** the 40 TeV production rate.
- A tiny POWHEG stage-1 initialization completed; it checked runtime and PDF
  loading, not production integration precision.
- `pipeline_smoke_release` generated, showered and analyzed 30 events in each
  of the four samples and produced an HTML report with 12 PNGs. It reused a
  previously validated signal LHE and ihixs record. Negative source weights
  were retained: signal 3, prompt diphoton 4, gamma+jet 8, Drell–Yan 4.
  Tree/summary weight closure was exact for backgrounds and below `1e-18`
  for signal. A representative plot was inspected; sparse negative bins do
  not support shape conclusions.
- Fourteen entry-point tests passed locally and on Timur; 111 existing ihixs,
  HO campaign and normalization tests passed locally. macOS used
  `TMPDIR=/private/tmp` to avoid equivalent-path fixture mismatches.

Fresh full POWHEG integration and the 109-point calculation were not repeated
for that earlier downstream test. Failed setup/pilot attempts were preserved.

## Existing ihixs production — October 5, 2026

A separate refined Timur calculation completed 109 points, both benchmarks,
PDF/coupling checks and repaired central comparison. At 40 TeV pp,
`mH=125 GeV`, `muR=muF=62.5 GeV`, pure HEFT N3LO QCD with the aN3LO QED PDF,
it obtained `226.82244161097535 +/- 0.06711554812896317 pb`.
The quoted error is numerical only; scale/PDF uncertainties are separate.
This is evidence for that saved calculation, not a substitute student record.

## Earlier diphoton checks

On September 17, 2026, tests used MG5 3.5.15, Herwig 7.3.0, ThePEG 2.3.0,
LHAPDF 6.5.6, ROOT 6.40.02 and HwSim. Fifty-three regression tests passed;
all three NLO background processes exported. Thirty-event signal and
Drell–Yan chains preserved signed weights and completed response analysis;
Drell–Yan response closure differed by `1.4e-14`. Resume and reports were
checked. Photon-background export alone did not establish their integration.
Herwig signal warnings included an event-level momentum mismatch of about
5 MeV; their presence remains a simulation-quality observation.

On July 16, 2026, both diphoton analyzers compiled on Timur and 40 Python
campaign/analysis tests passed locally and there. A 50,000-event signal
response check gave mean mass 124.997 GeV and RMS 0.877 GeV in 121–129 GeV,
with weight closure at floating-point precision. A near-Z dielectron check
reproduced fake-pair probability `1.8225e-6`; gamma+jet and synthetic tests
exercised flavour assignment, fake-rate interpolation and deterministic
smearing. These checks validate bookkeeping, not projected detector performance.

## Four-lepton reference checks

The [validation cards](hfourlepton/validation/) contain direct SM decay surveys
at `mH=125.09 GeV` using MG5 3.5.15. Widths `Gamma(4e)=Gamma(4mu)=1.353e-7 GeV`
and `Gamma(2e2mu)=2.436e-7 GeV` give fractions 0.263127188, 0.263127188,
0.473745624 (percent-level survey uncertainty), normalized to the external
`BR(H -> 4l)=1.251e-4`. They do not replace PROPHECY4F.

Earlier Herwig 7.3/HwSim forced-decay checks used ten stable-top events.
The prompt-only two-event factor was `0.0453183201610`; virtual W lines from
heavy-hadron decays received no extra prompt factor. A one-event full-support
semileptonic proposal gave finite weight `0.0286571568`, propagated to all
145 optional LHE weights. These establish handler placement and propagation,
**not** accelerated-sample physics closure. Such strata remain
`closure_validated: false` and excluded from nominal totals.

## Repeatable software checks

```bash
TMPDIR=/private/tmp python3 -m unittest discover -s tests -p 'test_ihixs*.py'
python3 -m unittest discover -s tests -p 'test_ho_*.py'
python3 -m unittest discover -s tests -p 'test_gammagamma*.py'
```

Use an activated runtime for tests requiring PyROOT or compiled analyzers.
The [four-lepton analyzer reference](hfourlepton/LOAnalysis/Code/README.md)
describes its core tests and synthetic ROOT fixture. Unit tests mock some
expensive boundaries; actual pilots must also check saved counts, signed
normalization, response closure and representative reports.
