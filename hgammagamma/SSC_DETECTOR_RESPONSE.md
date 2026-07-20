# SSC/GEM Detector Response For The Gamma-Gamma Analysis

This document is the implementation and operating reference for the optional
SSC/GEM detector response in the leading-order `H -> gamma gamma` workflow.
The historical inputs and their broader physics context are summarized in
[`SSCInfo.md`](../SSCInfo.md). The executable implementation is
[`LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC.cc`](LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC.cc).

## Quick Start And Option Semantics

The campaign has two detector-response choices:

```bash
# Recommended physics workflow: GEM smearing, efficiencies, and photon fakes.
python3 hgammagamma/run_gammagamma_campaign.py \
  --detector-response ssc \
  --run-tag run_02_ssc

# Original unsmeared HwSim photon analysis, for genuine-photon samples only.
python3 hgammagamma/run_gammagamma_campaign.py \
  --detector-response none \
  --run-samples signal_gg_h_aa,bkg_prompt_aa \
  --run-tag run_02_legacy
```

`ssc` is the default and can also be selected with
`DETECTOR_RESPONSE=ssc`. It runs `HwSimPostAnalysis_gammagamma_SSC`.
`none` runs the original `HwSimPostAnalysis_gammagamma` executable.

Important: `none` is a legacy baseline, not an otherwise identical SSC
analysis with only the Gaussian randomization disabled. The legacy executable
uses `pT > 10 GeV` and `|eta| < 6`, with no GEM transition veto, isolation, or
efficiency weights. The SSC executable uses the selection documented below.
Consequently, an `ssc` versus `none` comparison changes acceptance and object
selection as well as detector response. Use separate run tags because both
executables otherwise use the same output filename pattern.

The campaign rejects `bkg_gamma_j` and `bkg_dy_ee` under `none`: the legacy
executable has no jet-to-photon or electron-to-photon transfer model. It is
better to stop than to produce physically misleading zero-fake samples.

The cut and XGBoost analysis cards support a provenance guard:

```yaml
analysis:
  detector_response: ssc
```

The analysis then checks every campaign `.dat` file and stops if its recorded
detector profile does not match. Use `none` in a card intended for legacy
outputs.

## Simulation Chain And Samples

The automated chain is:

```text
MG5_aMC -> LHE -> Herwig shower/hadronization -> HwSim ROOT
        -> selected post-analysis executable -> plots/cuts/XGBoost
```

The campaign defaults are 20 TeV per proton beam, NNPDF40 LO through LHAPDF
ID 331900 in both MG5 and Herwig, and deterministic sample seeds
`seed_base + sample_index`. `HW-template.in` has `SavePartons Yes`, allowing
the gamma+jet response to distinguish quark and gluon jets using the hard
record.

The default samples are:

| Sample | MG5 process | Response topology |
| --- | --- | --- |
| `signal_gg_h_aa` | `g g > h [noborn=QCD]`, MadSpin `h > a a` | `genuine` |
| `bkg_prompt_aa` | `p p > a a` | `genuine` |
| `bkg_gamma_j` | `p p > a j` | `gammajet` |
| `bkg_dy_ee` | `p p > e+ e-` | `dielectron` |

The full `p p > e+ e-` sample contains on-shell Z and virtual-photon
Drell--Yan. The tree-level `p p > a a` sample does not contain the separate
loop-induced `g g > a a` continuum, which remains a future sample. Signal and
continuum interference is also not included.

## Electromagnetic Response

Photons, electrons promoted to photons, and jets promoted to photons are made
massless and smeared along their original direction. The implementation uses

\[
\left(\frac{\sigma_E}{E}\right)^2 =
\frac{a^2}{E} + b^2 +
\left(\frac{\sigma_{E_T}^{\rm noise}}{p_T}\right)^2,
\]

with energies and transverse momenta in GeV. The constants are:

| Region | Sampling term `a` | Constant term `b` | Thermal `ET` noise |
| --- | ---: | ---: | ---: |
| barrel, `|eta| < 1.01` | 0.060 | 0.004 | 0.100 GeV |
| outside the barrel | 0.085 | 0.004 | 0.175 GeV |

The pileup-noise term is 0.120 GeV for `|eta| < 1.4` and
`0.120 + 0.366 (|eta| - 1.4)` GeV outside it. Thermal and pileup terms are
combined in quadrature. `--no-pileup-noise` removes only the pileup term; the
thermal term remains. A Gaussian draw has a small positive energy floor. No
angular, pointing, vertex, or conversion smearing is applied.

Each object's random seed is a deterministic hash of the analysis seed,
source-event index, object collection, and object index. Results are therefore
reproducible and independent of conditional control flow through other object
collections.

## Jet Response

Ordinary reconstructed jets use

\[
\frac{\sigma_E}{E} =
\begin{cases}
\sqrt{0.60^2/E + 0.04^2}, & |\eta| \leq 3,\\
\sqrt{2.00^2/E + 0.06^2}, & |\eta| > 3.
\end{cases}
\]

The full jet four-vector is multiplied by the smeared-energy ratio. This
preserves its direction and `m/E`; it scales, rather than preserves, the jet
mass. A jet promoted to a photon instead becomes a massless EM candidate and
uses the electromagnetic response above. No jet angular smearing is applied.

## Acceptance And Isolation

All thresholds are applied after energy smearing, except the unchanged muon
diagnostic objects.

| Object or use | Minimum `pT` | Accepted `|eta|` | Additional condition |
| --- | ---: | ---: | --- |
| photon or fake-photon candidate | 20 GeV | 0.1 to 2.5 | veto `1.01 < |eta| < 1.16` |
| reconstructed jet diagnostic | 20 GeV | 0.1 to 5.5 | none |
| jet used for photon isolation | 10 GeV | 0.1 to 5.5 | veto candidate for `DeltaR < 0.4` |
| electron diagnostic | 10 GeV | 0.1 to 2.5 | same crack veto |
| electron promoted to photon | 20 GeV | 0.1 to 2.5 | crack veto and jet isolation |
| muon diagnostic | 10 GeV | 0.1 to 2.5 | no momentum smearing |

Both light jets and b-jets participate in the isolation collection. The
isolation is a deliberately coarse jet-based proxy, not the complete GEM
calorimeter and tracker isolation algorithm.

## Efficiencies

The genuine-photon probability is the product of shower identification and
the electron-veto acceptance:

| Quantity | Default |
| --- | ---: |
| photon shower ID | 0.90 |
| photon electron veto, continuum | 0.96 |
| photon electron veto, near Z | 0.86 |
| effective genuine photon, continuum | `0.90 * 0.96 = 0.864` |
| effective genuine photon, near Z | `0.90 * 0.86 = 0.774` |

“Near Z” means that the relevant pair mass satisfies
`|m - 91.1876 GeV| < 10 GeV` by default.

The electron, muon, and jet efficiency settings are diagnostic in this
gamma-gamma executable:

| Object | Default | Effect |
| --- | ---: | --- |
| electron | 0.90 | weighted before/after count in `.dat` only |
| muon | `0.85 * 0.95 = 0.8075` | weighted before/after count in `.dat` only |
| jet | 1.00 | weighted before/after count in `.dat` only |

They do not thin objects and do not alter diphoton hypotheses. GEM does not
give a universal offline jet efficiency for this analysis, so 1.00 is an
explicit baseline assumption. Muons are not smeared in the present code.

## Jet-To-Photon Fakes

The `gammajet` topology uses one genuine photon and one selected hard-jet fake.
The rates are linearly interpolated in reconstructed diphoton mass and clamped
to the endpoint values:

| `m(gamma gamma)` [GeV] | 80 | 100 | 120 | 140 | 160 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| quark jet | 4.1e-4 | 3.9e-4 | 3.7e-4 | 3.6e-4 | 3.4e-4 |
| gluon jet | 1.2e-4 | 1.1e-4 | 1.0e-4 | 0.91e-4 | 0.83e-4 |

These GEM rates are already post-isolation and post-photon-identification; no
extra 0.90 shower-ID factor is applied to a fake jet. `--jet-fake-scale`
multiplies the tabulated probability.

The selected fake candidate is matched to a hard quark or gluon within
`DeltaR < 0.4` using the unsmeared jet direction. If no geometric match is
found but the hard gamma+jet record has exactly one colored outgoing parton,
that parton supplies the flavor. If flavor is still unknown, the probability
is a mixture controlled by `--unmatched-jet-quark-fraction`, whose conservative
default is 1.0. The `.dat` counters distinguish geometric matching, hard-record
fallback, and unknown flavor.

## Electron-To-Photon Fakes

The `dielectron` topology promotes an isolated reconstructed electron and
positron to photon candidates. The tracker-veto probability per electron is:

| Source-pair region | Tracker-veto probability |
| --- | ---: |
| `|m(ee) - 91.1876 GeV| < 10 GeV` | 0.0015 |
| continuum | 0.020 |

Each probability is multiplied by the independent 0.90 EM shower-ID
efficiency and by `--electron-fake-scale`. Thus the default near-Z two-fake
probability for an accepted pair is
`(0.90 * 0.0015)^2 = 1.8225e-6`. The mass-region choice uses the unsmeared
source `m(ee)`. The model assumes independent fakes. It does not additionally
apply electron reconstruction efficiency or genuine-photon track-veto
acceptance.

## Weighted Exclusive Outcomes

ID and fake probabilities are represented deterministically as exclusive
weighted hypotheses, rather than rare Bernoulli throws. For two accepted
sources with probabilities `p1` and `p2`, the tree contains outcomes with

\[
P_0=(1-p_1)(1-p_2),\quad
P_{1a}=p_1(1-p_2),\quad
P_{1b}=(1-p_1)p_2,\quad
P_2=p_1p_2.
\]

If only one source passes kinematic and isolation acceptance, the zero- and
one-candidate outcomes are retained. For every source event,

\[
\sum_{\rm output\ rows} w_{\rm event}
= w_{\rm HwSim}\,w_{\rm scale}.
\]

One source event can therefore create several ROOT rows. Counts such as
`events_with_two_selected_photons` and `valid_diphoton_hypotheses` count valid
two-candidate rows, not expected event yields. Physical yields must use the
sum of `eventweight`. Cut mode reads the complete tree and the baseline card
explicitly requires `n_selected_photons >= 2`. XGBoost filters to those
two-candidate rows before training and scoring. Both workflows retain the
complete-tree weight sum as the normalization denominator.

## Configuration Reference

Campaign-level controls are:

| CLI | Environment | Meaning |
| --- | --- | --- |
| `--detector-response ssc|none` | `DETECTOR_RESPONSE` | choose SSC/GEM or legacy post-analysis |
| `--seed-base N` | `SEED_BASE` | base seed; sample index is added |
| `--run-samples LIST` | `RUN_SAMPLES` | select campaign samples |
| `--run-tag TAG` | `RUN_TAG` | output/generation tag |

The SSC executable accepts the following direct controls:

| Option | Default | Meaning |
| --- | ---: | --- |
| `-t TAG` | empty | output tag |
| `-n N`, `-nmin N`, `-nmax N` | all events | event range |
| `-w X`, `--weight-scale X` | 1 | sample-wide event weight |
| `--response-mode genuine|gammajet|dielectron` | `genuine` | source topology |
| `--seed N` | 14101983 | deterministic smearing seed |
| `--photon-id-efficiency X` | 0.90 | photon shower ID |
| `--photon-electron-veto-efficiency X` | 0.96 | continuum track veto |
| `--photon-electron-veto-efficiency-near-z X` | 0.86 | near-Z track veto |
| `--electron-efficiency X` | 0.90 | electron diagnostic efficiency |
| `--muon-efficiency X` | 0.8075 | muon diagnostic efficiency |
| `--jet-efficiency X` | 1.00 | jet diagnostic efficiency assumption |
| `--electron-fake-rate-near-z X` | 0.0015 | tracker-veto probability |
| `--electron-fake-rate-continuum X` | 0.020 | continuum tracker-veto probability |
| `--near-z-half-width X` | 10 GeV | mass half-window |
| `--jet-fake-scale X` | 1 | jet-fake systematic multiplier |
| `--electron-fake-scale X` | 1 | electron-fake systematic multiplier |
| `--unmatched-jet-quark-fraction X` | 1 | unknown-flavor mixture |
| `--no-pileup-noise` | off | omit pileup, retaining thermal noise |

The campaign assigns the response mode, seed, tag, and weight scale. To vary
the other response parameters, rerun the SSC executable directly on an
existing `<sample>_hwsim_roots.input`, using a new tag, or add the desired
arguments to a controlled campaign variation.

## Output Contract And Provenance

For an input list named `<sample>_hwsim_roots.input` and tag `<tag>`, the
post-analysis writes:

| Output | Contents |
| --- | --- |
| `*-<tag>.top` | weighted TopDrawer histograms |
| `*-<tag>.dat` | analysis/profile identity, settings, counters, diagnostics, and closure |
| `*-<tag>.evp` | source indices that have a valid two-candidate hypothesis |
| `*-<tag>_var.root` | `Data2` analysis tree |

The SSC tree preserves the legacy ten-element `variables` array and adds:

```text
generatorweight[0] = HwSim event weight * sample weight_scale
responseweight[0]  = exclusive detector-response probability
eventweight[0]     = generatorweight[0] * responseweight[0]
photonorigin[0:2]  = 0 none, 1 genuine photon, 2 jet, 3 electron
sourceevent[0]     = original input-chain entry
```

Despite its name, `generatorweight` includes the sample-wide `weight_scale`.
The `.dat` file records `analysis`, `detector_response`, `response_mode`,
`weighted_hypotheses`, seed, efficiency/fake settings, q/g counters, object
diagnostics, `sum_weight`, `sum_tree_weight`, and
`tree_weight_closure_difference`.

Reports expose this provenance in HTML and CSV. Historical outputs without an
explicit `detector_response` are supported: a known SSC analysis name is
inferred as `ssc`, a metadata-free legacy file as `none`, and an unrecognized
future analysis as `unknown`. A plot report warns if one run tag mixes detector
profiles. Analysis cards with `detector_response` reject such a mismatch.

For physical rate plots, use:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --run-tag run_02_ssc \
  --no-density \
  --normalization event_xsec
```

Density plots are unit-area shapes and intentionally hide rate differences.

## Validation Record

The following development checks were run on 2026-07-16. They are technical
validation of implementation and bookkeeping, not a detector-performance or
physics validation of the phenomenological model.

- Both post-analysis executables compiled against the Herwig/HwSim ROOT stack
  on `timur.kennesaw.edu`.
- All 40 gamma-gamma Python campaign and analysis unit tests passed locally
  and on Timur, including response selection, provenance guards, and
  weighted-row reading. Executable smoke checks and synthetic ROOT inputs
  collectively cover q/g rates, electron fakes, deterministic smearing, and
  weight closure.
- A 50,000-event signal response test gave mean `m(gamma gamma) = 124.997 GeV`
  and RMS `0.877 GeV` in the 121--129 GeV window. The input-weight sum was
  50,000 and the output-tree sum agreed within floating-point precision.
- A near-Z Drell--Yan smoke test reproduced the expected pair probability
  `1.8225e-6` with zero reported closure difference.
- A gamma+jet smoke test exercised hard-parton flavor assignment and its
  geometric/fallback counters with no unknown-flavor hypotheses.
- Separate synthetic-input executable checks exercise both quark- and
  gluon-rate interpolation and endpoint behavior.

Recommended repeatable checks after any response change are:

```bash
python3 -m unittest tests.test_gammagamma_campaign tests.test_gammagamma_analysis

python3 hgammagamma/run_gammagamma_campaign.py \
  --dry-run --nevents 10 --detector-response ssc

python3 hgammagamma/run_gammagamma_campaign.py \
  --dry-run --nevents 10 --detector-response none \
  --run-samples signal_gg_h_aa,bkg_prompt_aa
```

On the current Timur checkout, add
`--mg5-dir /home/apapaefs/Projects/MG5_aMC_v3_5_15` if the clone-local MG5
directory is absent.

For SSC output, require
`abs(tree_weight_closure_difference)` to be consistent with accumulated
floating-point roundoff relative to `sum_weight`.

## Known Limitations

- No angular/pointing, vertex, conversion, or detector-material model is
  included.
- Pileup is represented only by an EM noise term; there are no overlaid events,
  pileup jets, occupancy effects, or pileup-dependent isolation.
- Muon momenta are not smeared.
- Jet, electron, and muon efficiencies are diagnostic counts only; only
  photon efficiencies and photon-fake probabilities enter diphoton weights.
- The offline jet efficiency is fixed to an assumed 100% by default.
- `gammajet` selects one hard fake candidate rather than summing over every
  shower jet as an independent fake source.
- There is no `dijet` double-fake response mode yet.
- Isolation is a simple `DeltaR` veto against reconstructed jets, not the full
  GEM shower-shape, calorimeter, tracker, and underlying-event procedure.
- Object-ID probabilities are independent; correlated efficiencies and
  systematic nuisance models are not implemented.
- The input arrays are capped at 100 objects per collection. Truncation is
  counted in `.dat`.
- The default samples are leading order. Higher-order background corrections,
  the loop-induced prompt continuum, and signal-background interference are
  not included by default.
- No standard detector-response uncertainty envelope has yet been defined;
  the fake-rate scale options are variation hooks, not a complete uncertainty
  prescription.
- The legacy `none` option is not an apples-to-apples no-smearing ablation of
  the SSC profile. Such a study would require a new SSC option that preserves
  the SSC acceptance, isolation, efficiencies, and weighted bookkeeping while
  disabling only energy randomization.

## Numerical Sources

The implemented numerical baseline comes from:

- GEM Collaboration, *GEM Technical Design Report*, SSC-GEM-TN-93-262,
  SSCL-SR-1219 (1993):
  <https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-262.pdf>.
- S. Mrenna et al., *Higgs Searches with the GEM Detector*, GEM TN-93-373,
  CALT-68-1856 (1993):
  <https://lss.fnal.gov/archive/other/calt-68-1856.pdf>.

See [`SSCInfo.md`](../SSCInfo.md) for the provenance of individual forecast
numbers and for comparison with other SSC- and LHC-era studies.
