# Recover a campaign

First run `python3 run_pipeline.py status --workflow WORKFLOW --run-tag TAG`.
Read the named stage log. Keep the failed directory and manifests: they contain
the evidence needed to resume safely.

| Symptom | Meaning and next step |
| --- | --- |
| No `.top`/`.dat` files when making a report | Generation/showering may be complete but analysis is missing. Repeat the original `run` command; it checks and resumes stages. |
| Missing normalization record | A clone contains source, not a calculated rate. HO `run` builds, benchmarks and calculates ihixs before analysis. |
| Configuration changed | Use a new tag for changed event counts, cards or scientific settings. Do not remove the manifest to force reuse. |
| Unsupported option | Check `git log -1 --oneline`, `git remote -v` and the exact script's `--help`; update the checkout actually being run. |
| MG5 error with exit status 0 | MG5 can report failure inside its log. Inspect `mg5-process.log` or `mg5-generate.log`; the exit status alone is insufficient. |
| Ninja download/build failure | Use common `setup`, which configures the private MG5 and available loop libraries. Do not replace a completed process directory. |
| PDF missing or shared PDF directory unwritable | Repeat common `setup` for the required workflow. It installs PDFs privately; do not write into `/home/shared`. |
| Long silence | Integration/build stages can be quiet. Check the log timestamp and process CPU activity; silence alone does not establish a stalled job. |
| Missing run tag | Check spelling and the workflow: each channel has its own run directory. |

## POWHEG consistency failure

`The fraction of inconsistent runs is too large` is a failed integration
consistency check. A traceback naming stage 3 does not mean stage 3 caused it:
inspect all stage-2 logs and the seed identified by POWHEG. Common HO setup uses
`ncall2=500000` and retains `check_bad_st2=1`.

Preserve the failed integration. A controlled retry may reuse validated
stage-1 grids in a new run directory and increase stage-2 statistics. Do not
reuse suspect stage-2 grids, disable the check, or treat a new random seed as
proof that the original inconsistency was resolved. The advanced wrapper is
`hgammagamma/run_powheg_hjminnlo.py`; inspect its `--help` and record the retry
settings.

## ihixs integration failure

The normal HO command resumes completed points and refines failed points.
Keep the record and complete calculation directory together. A missing or
failed Cuba result is not a usable cross section. The raw point logs distinguish
precision failures from configuration/build errors; see the
[ihixs reference](../hgammagamma/HOAnalysis/IHIXS_STUDENT_GUIDE.md).

## Herwig event errors or exhausted input

Finite LHE input must not be replayed to fill a requested count. HO completion
checks require finalized ROOT output, source-consumption evidence and signed
weight accounting. Discarded source events remain in the normalization
population with zero simulated response. A terminal error is recoverable only
when those checks succeed; warnings and event losses remain quality diagnostics.
For a changed shower configuration, use a new tag and retain the old evidence.

## Existing campaigns

New common runs are isolated by tag. Existing HO and four-lepton manifests
can resume through their original entry points and saved settings. Legacy LO
outputs without manifests remain readable by reports and analysis; start new
generation under a new tag. Copying outputs does not transfer provenance.
Reports and cuts can be directed to existing saved analysis through
`--analysis-root` and `--output-dir`; see [analysis](ANALYSIS.md).
