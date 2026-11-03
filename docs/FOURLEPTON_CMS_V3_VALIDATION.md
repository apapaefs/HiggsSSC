# CMS four-lepton ancestry validation — 3 November 2026

The HwSim ancestry correction passes the existing approximate CMS acceptance
benchmark for all five Higgs production modes. The selection, detector response
and normalization are distinct: this is a particle-level comparison with
perfect response, not a CMS detector simulation or an SSC resolution validation.

## Object definition and implementation checks

The private `HwSimAncestry.so` plugin records the actual parent graph and ThePEG
particle-copy links. CMS v3 uses direct hard/H/Z/W electrons and muons, excludes
hadron and tau daughters, dresses with photons from the emitting lepton's
lineage within delta-R 0.3, and excludes only that lepton's FSR from isolation.
Hadronic photons contribute to isolation. The shared installed HwSim library,
the original completed campaign, and the original analyzer executable are
unchanged. Implementation and student commands are in the
[ancestry workflow guide](../hfourlepton/hwsim/README.md).

The smoke test used 200 events each from ggF→4mu, qqZH→4e, ttH→2e2mu and
qq→2e2mu. All 800 events reproduced the old stable-particle momenta and generator
weights exactly. The analyzer rejects missing or malformed ancestry, including
mixed legacy/new ROOT inputs. Core tests cover photon lineage, hadronic photons,
tau daughters, continuum leptons without an explicit Z, and isolation overlap.

The full `fourl13_minimal_nlofix/revisions/cms_v3` pass reuses the original hard
and MadSpin events. All 27 showers and analyses completed with unchanged source
counts and signed-weight moments. ihixs normalization for ggF, native rates for
other modes, physical branching fractions, and failed showers' zero response
are unchanged. The run consumed 231,000 hard events and saved 230,966 showers.
The strict report has zero errors and zero warnings.

## Signal acceptance

The completed 21 signal samples give the following fiducial acceptances in the
105–140 GeV four-lepton mass range. Errors in the simulation column are the
existing analysis's Monte Carlo errors. Reference values are from
[CMS HIG-19-001, Table 8](https://arxiv.org/pdf/2103.04956).

| Mode | CMS (%) | Legacy v2 (%) | Ancestry v3 (%) | Relative difference from CMS | Selected effective events |
|---|---:|---:|---:|---:|---:|
| ggF | 40.20 | 42.53 | 40.48 ± 0.37 | +0.69% | 12,210 |
| VBF | 44.50 | 46.60 | 44.98 ± 0.52 | +1.07% | 7,459 |
| WH | 32.90 | 35.64 | 33.98 ± 0.42 | +3.29% | 6,611 |
| ZH | 34.00 | 37.91 | 36.20 ± 0.56 | +6.48% | 4,152 |
| ttH | 31.50 | 35.92 | 33.56 ± 0.85 | +6.54% | 1,541 |

Every mode passes the unchanged requirements of 10% relative agreement and
at least 1,000 selected effective events. ZH and ttH remain above the CMS central
values by about 6.5%; this is a broad implementation benchmark, not precision
agreement with the experimental uncertainties. Associated-production
combinatorial signal is not separately fitted, and no ATLAS likelihood or
matrix-element pairing is reproduced. No threshold or normalization was tuned
to obtain this result.

## Background statistics extension

With CMS v3, the 120–130 GeV region at 100 fb^-1 contains 247.24 ± 1.97 signal
and 86.35 ± 8.45 irreducible-background events (Monte Carlo errors only). The
background effective count remains about 104, below the required 400.
In the preceding audit qq→4l supplied 99.94% of its variance. An independent campaign,
`fourl13_qq200k_truth`, is generating 200,000 events for each qq final state
(4e, 4mu, 2e2mu), using CMS v3 and the same NLO matching, PDFs and 90–180 GeV
generation window. The command specifies seed base 20261103; the minimal preset
derives and records the campaign-specific effective seed. Generation uses
32 cores while up to eight ancestry-validation workers run alongside it.

The new qq samples replace the previous three channels in the combined
manifest. They are not added as extra physical cross sections. The assembly
helper verifies both input campaigns and physics settings before running the
strict report, acceptance check, cut analysis and XGBoost. The effective-count
requirements remain 400 for the cut analysis and 100 for the ML working point.
The final statistical support must be read from those outputs after generation
finishes; the projected improvement is not a measured result.

The intended combined output is
`hfourlepton/runs/fourl13_minimal_nlofix/revisions/cms_v3_qq200k/` on Timur.
Generation and validation launch logs are under `inputs/pipeline/logs/`.

## Reproduction identifiers

- Implementation commit: `4ba92dac78fd83d010f0acfe53dfa57dbbf13ca4`.
- HwSim base revision: `51ba33aea55ddbee0df8309d425b6fbdde1ceee3` plus the tracked ancestry patch.
- Private plugin SHA-256: `94f7b191523cba59d774bf8711afdc3caec32afa2d023f5a52a4a2819134417a`.
- CMS v3 manifest SHA-256: `4659859ccc6d8b72659dfd6a9d2c87a43c245d2f627441aa2ec12a44cc491725`.
- Acceptance result: `pass_approximate_acceptance_check`.
