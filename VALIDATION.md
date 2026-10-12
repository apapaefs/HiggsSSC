# Validation record

Pilots check execution, event accounting and outputs. They do not establish
precision, detector realism, matching/isolation stability or discovery reach.
Keep raw logs, manifests, cards and products with any numerical result.

## Common workflow audit — October 2026

Validation uses an isolated instructor-owned Timur checkout and a fresh private
installation. No generated LHE, grids or normalization records were borrowed;
the shared Herwig/HwSim runtime is read-only. Failed attempts are preserved.

| Acceptance check | Status | Evidence |
| --- | --- | --- |
| Private setup and additive HO setup | Passed | Own model source, all 101 LO PDF members, then HO dependencies |
| Four-sample LO pilot, resume, report and cuts | Passed | Counts, signed weights and hashes checked below |
| Fresh ihixs build, benchmark and normalization | Passed | 109 points and 443 raw artifact hashes checked |
| Fresh full POWHEG integration and signal | Pending | In progress; not certified by earlier initialization tests |
| Four-sample HO pilot and report | Pending | Requires completed fresh signal input |
| Nine four-lepton samples, both responses, strict report | Pending | Full campaign acceptance remains in progress |

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
- 441 local regression tests passed; seven optional ML tests were skipped.
  Documentation links, shell examples and dependency-free help were checked.

## Reference checks and scientific limits

Earlier downstream HO pilots exercised signed weights in all four samples
and closed response weights to floating-point precision, using existing signal
LHE/normalization inputs. Those checks support bookkeeping but do not certify
the pending fresh pipeline above. Numerical tests of smearing, fake-rate
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
