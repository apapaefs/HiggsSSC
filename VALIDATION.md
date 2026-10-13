# Validation record

Pilots check execution, event accounting and outputs. They do not establish
precision, detector realism, matching/isolation stability or discovery reach.
Keep raw logs, manifests, cards and products with any numerical result.

## Common workflow audit — 12 October 2026

Validation used two isolated instructor-owned Timur checkouts with fresh private
installations. No generated LHE, grids or normalization records were borrowed;
the shared Herwig/HwSim runtime was read-only. Failed attempts are preserved.

| Acceptance check | Status | Evidence |
| --- | --- | --- |
| Private setup and additive HO setup | Passed | Own model source, all 101 LO PDF members, then HO dependencies |
| Four-sample LO pilot, resume, report and cuts | Passed | Counts, signed weights and hashes checked below |
| Fresh ihixs build, benchmark and normalization | Passed | 109 points and 443 raw artifact hashes checked |
| Fresh full POWHEG integration and signal | Passed | Four seeds; all integration stages and consistency checks passed |
| Four-sample HO pilot and report | Passed | 30 source events per sample; signed-weight closure and 489 hashes checked |
| Nine four-lepton samples, both responses, strict report | Passed | 18 sample/profile combinations; independent weight audit passed |

### Completed measurements

- `lo_pilot_02` completed all four samples with 30 LHE, 30 HwSim and 30
  analyzed source events each. Source weights were positive. Response rows
  numbered 81 (signal), 40 (prompt diphoton), 38 (gamma+jet), 63 (Drell–Yan);
  these are weighted detector outcomes, not independent events. Tree/summary
  closure was within `1.5e-14`.
- All 41 manifest artifact hashes/sizes and 16 pilot hashes matched. Resume
  reused all 12 generation/shower/analysis stages. Report-only regeneration
  at 200 fb⁻¹ preserved pilot evidence. The report had 12 PNGs, 12 SVGs and
  13 CSVs; a representative plot was inspected.
- Baseline cuts at 100 fb⁻¹ used all four samples. No background hypotheses
  passed the mass window; significance correctly remained undefined. This
  pilot cannot establish background rejection or sensitivity.
- Fresh ihixs completed central + six scale variations + 100 PDF replicas +
  two native-order comparisons. At 40 TeV pp, `mH=125 GeV`,
  `muR=muF=62.5 GeV`, pure-HEFT N3LO with
  `NNPDF40_an3lo_as_01180_qed_mhou`, the rate was
  `226.82098 +/- 0.07033 pb` (numerical error only). Independent reparsing
  matched all 109 raw results/coupling checks; all 443 raw artifact hashes
  matched. Recomputed uncertainties were scale `-9.24092/+1.57075 pb` and
  PDF replica standard deviation `1.92285 pb`.
- `scale_2` and `scale_4` each required one precision retry: `epsrel` changed
  from `1e-4` to `1e-6`, retaining `maxeval=50000000`. Original attempts remain
  saved. Every final relative numerical error was below the 0.05% gate.
  The separate 13 TeV upstream benchmark with `PDF4LHC15_nnlo_100` returned
  `45.1815532 +/- 0.0109436459 pb`, compatible with reference 45.1816 pb.
- Fresh POWHEG used four seeds with 8, 8, 7 and 7 events, `ncall2=500000`
  and `check_bad_st2=1`. `ho_pilot_01` completed all four samples with 30 LHE,
  HwSim and analyzed source events each. Negative source-event counts were
  2 (signal), 4 (prompt diphoton), 8 (gamma+jet) and 4 (Drell–Yan).
  Response rows numbered 71, 45, 36 and 68 respectively; aggregate closure
  was within `1.1e-14`. Signal `BR=0.00227` was applied once.
- The HO audit verified 489 unique file hashes, including all normalization
  artifacts and 20 pilot outputs. Its report contained 12 PNGs, 12 SVGs and
  13 CSVs; a representative plot was inspected. Sparse signed samples can
  produce negative bins; this is not a precision prediction.
- `fourl_pilot_03` completed nine samples with both `perfect` and `ssc`
  responses: 30 source events per sample/profile, 540 analyzed entries in
  total. The nine profile pairs reuse identical LHE inputs and normalization;
  these are not 540 independent generated events. Seven hard-process groups
  were generated. All source weights in this LO pilot were positive.
- Independent XML/ROOT checks matched all 76,680 finite optional weights
  exactly through LHE, HwSim and analysis. Per-source response closure was
  within `2.23e-16`; retained MadSpin weight ratios agreed within `6.8e-8`
  serialization rounding. Nine undefined upstream signal alternatives were
  explicitly recorded. Selected cross sections matched the strict report,
  which had zero audit errors. Representative plots were inspected.
- Four-lepton resume preserved all 38 analysis hashes and the pilot record.
  All 188 manifest artifact records matched across 116 unique files. On
  source `4560999`, its regression suite ran 460 tests with three skips for
  ihixs source-anchor checks, outside that four-lepton-only installation.

### Recovery and final checks

Real HO resume exposed a byte-identical POWHEG merge rewrite that changed its
timestamp and correctly failed shower provenance validation. The merge now
preserves identical files and refuses conflicting output; failed attempts are
retained. The affected LHE's recorded timestamp was restored only after its
original SHA-256 was independently verified, with a backup and recovery log.
No event content or scientific manifest was changed.

On final code `d6f4c3a`, the real four-seed merge and identical-output retry
passed. HO resume reused all completed stages; report-only regeneration at
200 fb⁻¹ preserved all 16 LO and 20 HO physics hashes and the LHE's identity.
Additive four-lepton/ML setup also passed in the LO→HO installation, including
private repair of its incomplete inherited NumPy at the same version.
The full Timur suite ran **466 tests, all passing with no skips**. Documentation
checks covered 18 tables, 71 local links, 22 shell examples, the embedded Python
example and dependency-free help. External references were checked directly or
against primary-source citations where CERN's site blocked automated access.

## Reference checks and scientific limits

Numerical tests of smearing, fake-rate
interpolation, near-Z dielectron probability `1.8225e-6`, flavour assignment
and response closure are described by the
[detector model](hgammagamma/SSC_DETECTOR_RESPONSE.md). They are not an
experimental validation of GEM performance. Herwig event losses and momentum
warnings remain quality diagnostics even when output completion is accepted.

The [four-lepton validation cards](hfourlepton/validation/) use direct SM
matrix elements at `mH=125.09 GeV` in MG5 3.5.15. Widths
`Gamma(4e)=Gamma(4mu)=1.353e-7 GeV`, `Gamma(2e2mu)=2.436e-7 GeV` give channel
fractions 0.263127188, 0.263127188, 0.473745624 with percent-level survey
uncertainty. Their sum is normalized to external `BR(H -> 4l)=1.251e-4`;
the survey does not replace PROPHECY4F.

Forced-decay smoke checks established handler placement and propagation to
optional LHE weights, not accelerated heavy-flavour physics closure. Those
strata remain `closure_validated: false` and excluded from nominal totals.
Fresh stable-H generation also exposed undefined alternative scale weights;
see the [four-lepton variation policy](hfourlepton/README.md#optional-weight-provenance).
Preserving valid weights does not establish a complete uncertainty ensemble.

## Repeatable checks

With the configured Python/ROOT runtime:

```bash
python3 -m unittest discover -s tests
```

On macOS, prepend `TMPDIR=/private/tmp` to avoid equivalent-path fixture
mismatches. The [four-lepton analyzer reference](hfourlepton/LOAnalysis/Code/README.md)
provides core tests and synthetic ROOT fixtures. Unit tests mock expensive
boundaries; actual pilots must also verify saved counts, signed normalization,
response closure and representative reports.
