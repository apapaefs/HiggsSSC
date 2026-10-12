# Full HO Diphoton Pipeline On Timur

Each student can install the additional dependencies, generate the signal
and backgrounds, calculate the ihixs normalization, shower, analyze and
make the report in their own account. No files from the instructor's home
directory are required. Use
[`run_gammagamma_ho_pipeline.py`](../run_gammagamma_ho_pipeline.py) for this
workflow. The existing campaign runner remains available for individual
stages and recovery of older campaigns.

This workflow targets **pp collisions at 40 TeV**, HJMiNNLO NNLO+PS Higgs
event shapes, the inclusive **pure-HEFT N3LO ggF** normalization from ihixs,
and the three MC@NLO backgrounds. The analysis applies the Higgs branching
fraction and SSC/GEM detector response separately. A small pilot checks the
workflow; it does not establish statistical precision or physics readiness.
See the [physics conventions](README.md#processes-and-conventions).

## 1. Use Your Own Checkout

Log into Timur, then clone into a directory owned by your account:

```bash
mkdir -p "$HOME/Projects"
cd "$HOME/Projects"
git clone https://github.com/apapaefs/HiggsSSC.git
cd HiggsSSC
```

If you already have a checkout, use that repository directory and update
the source with `git pull --ff-only origin main`. Do not clone over an
existing directory. The pipeline initializes only the public ihixs
submodule; access to the paper is unnecessary.

Keep a long-running terminal alive with `screen -S higgs-ho` or `tmux`.
Run all following commands from the `HiggsSSC` repository root on Timur.

## 2. Install And Check Once

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py setup --jobs 4
```

The setup loads `herwig/stable` automatically and checks its Herwig,
HwSim, ROOT/PyROOT, LHAPDF, FastJet and compiler environment. These are
Timur's shared runtime prerequisites; the pipeline does not install or
modify the system Herwig/ROOT stack. On another Linux host, provide a
compatible stack and select `--herwig-env /path/to/Herwig` during setup,
with `--cc`, `--cxx` and `--fc` if needed. That host still needs its own
validation.

Setup installs under `hgammagamma/HOAnalysis/inputs/pipeline/`:

- A Python virtual environment using the activated PyROOT installation.
- The full NLO, NNLO and aN3LO QED PDF sets, plus the ihixs benchmark PDF.
  Their index and all 101 members per set go into a private PDF directory,
  which is placed first in the search path.
- MadGraph5_aMC@NLO 3.5.15 from a pinned upstream archive. On Timur it uses
  the available shared Ninja/OneLOop libraries read-only; otherwise MG5's
  bundled reduction libraries are used. The optional Ninja/Collier
  auto-download is disabled. An NLO process export checks this setup.
  A recorded startup repair imports `logging` before MG5's Python-version
  warning, which otherwise crashes the pinned release on Timur's Python 3.9.
  The archive also follows upstream release packaging by retaining
  `create_release.py` outside `bin/`; leaving it there would insert ANSI
  development-banner escapes into LHE headers. Automatic updates are
  disabled for this pinned installation.
- Cuba 4.2.2, compiled with double precision and installed privately.
- Pinned POWHEG-BOX-V2/HJ source and the compiled HJMiNNLO executable.

It also initializes the pinned ihixs source and builds the repository's
analysis executable. Cuba and MG5 archives are checked against pinned
SHA-256 hashes; POWHEG and HJ use fixed Git commits. Setup records these
versions, Python packages and binary hashes in `installed.json`.
Installation and generated results are
ignored by Git. Allow several GB of disk space and time for source builds.

Setup can be repeated after an interrupted download or build. It retains
previous logs and requires the same saved runtime/compiler settings. A
different installation can use a new `--work-dir`; pass that same option
to subsequent commands.

To check an installed environment without generating events:

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py check
```

## 3. Run The Complete Pilot

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py run \
  --run-tag ho_pilot_01 --nevents 30 --jobs 4
```

This command performs the entire chain:

1. Build the benchmark and production ihixs variants, run the benchmark,
   and calculate all 109 required normalization points. Failed precision
   points use the existing bounded refinement procedure; failed checks
   prevent installing a normalization record.
2. Generate a POWHEG HJMiNNLO signal LHE. The pipeline uses `ncall2=500000`
   for stage 2 and retains `check_bad_st2=1`.
3. Generate the three MG5 backgrounds, shower all four samples with Herwig,
   and run the response analysis, including signed-weight and normalization
   checks performed by the campaign runner.
4. Produce the report using `event_xsec` normalization and 100 fb^-1.

**The first run includes the full ihixs calculation and POWHEG grid
integration even for 30 events.** Event count limits the event sample,
not the normalization calculation. Later runs reuse the validated ihixs
record in `inputs/pipeline/normalization/ggf-ssc40-n3lo.json`. Keep it
together with its `ihixs/` calculation directory.

The report is:

```text
hgammagamma/HOAnalysis/plots/ho_pilot_01/index.html
```

Each command prints its log path. For a read-only summary:

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py status --run-tag ho_pilot_01
```

Integrations can be quiet between completed points. Follow the log path
printed by the command in a second terminal, for example:

```bash
tail -f hgammagamma/HOAnalysis/inputs/pipeline/logs/ihixs-calculate.log
```

Use `Ctrl-C` in the log-viewing terminal to stop following the log; leave
the calculation terminal running.

To inspect the plan without installing dependencies or running anything,
append `--dry-run` to `setup` or `run`.

## 4. Run A Larger Sample

After a successful pilot, choose a new tag:

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py run \
  --run-tag ho_run_01 --nevents 10000 --jobs 4
```

The runner requires the recorded pilot outputs to remain intact before
starting a run larger than 30 events. Review its event counts, signed
weights, normalization and a representative plot before interpreting a
larger sample. Repeating the exact same run command resumes completed
stages. Changing event counts, seed jobs, input paths or the POWHEG card
requires a new tag; the old results remain available. A lock prevents two
commands from using the same dependency/normalization directory at once.

## Optional: Reuse Inputs You Already Calculated

An existing validated ihixs record and/or compatible signal LHE can be
supplied explicitly. Use actual paths; both options are optional:

```bash
python3 hgammagamma/run_gammagamma_ho_pipeline.py run \
  --run-tag ho_reuse_pilot_01 --nevents 30 --jobs 4 \
  --signal-lhe "$PWD/hgammagamma/HOAnalysis/inputs/powheg-hjminnlo-merged.lhe" \
  --signal-normalization "$PWD/hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json"
```

The same strict normalization and LHE checks still apply. Omitting these
options computes the inputs in your own account; copying the instructor's
normalization is not a prerequisite.

## Existing Campaigns And Failures

For an older campaign that already completed `generate` and `shower`, use
the existing campaign runner's `--stage analyze --resume` with its original
options, followed by the report command. An empty report search means the
analysis outputs are missing. Installing Cuba alone does not create them.
Calculate the missing normalization using the
[ihixs guide](IHIXS_STUDENT_GUIDE.md) if necessary. This pipeline refuses
to take over a nonempty campaign that it did not create.

If POWHEG reports **“The fraction of inconsistent runs is too large”**,
inspect its stage-2 logs and grids. Higher sampling reduces this problem
but does not guarantee statistical consistency. The pipeline stops and
preserves the run; it never disables the consistency check or silently
drops a suspicious seed. Repeating `--resume` reuses completed stage-2
grids and therefore cannot by itself repair inconsistent grids. Use a
separate integration attempt with adequate stage-2 sampling, then supply
its validated LHE to a new pipeline tag. Preserve the failed grids and
logs so the inconsistent seed can be investigated.

For a build or MG5 error, inspect the printed log before rerunning. MG5
errors are detected even if it returns exit status zero. Failed setup
attempts keep their logs and extracted sources for inspection.

## Validation Of This Entry Point

On October 12, 2026, setup was exercised in a separate Timur checkout with
Python 3.9.25 and the `herwig/stable` runtime. Private Cuba, MG5, all four
PDF sets, POWHEG and the analysis executable were installed and checked.
The newly compiled POWHEG executable completed a technical stage-1
initialization with 100 calls; those grids are not physics results.

Fresh upstream and production ihixs builds linked against the private
Cuba installation. The 13 TeV upstream pure-HEFT benchmark with
`PDF4LHC15_nnlo_100` returned `45.1815532 +/- 0.0109436459 pb`, agreeing
with the expected `45.1816 pb` within its required tolerance. This is the
benchmark, not the 40 TeV production normalization.

The `pipeline_smoke_release` pilot completed all four samples through
generation, showering, analysis and the HTML report, producing 12 PNG
plots. Each sample saved and analyzed 30 events. The consumed source
populations retained 3 negative-weight signal events, 4 prompt-diphoton
events, 8 gamma+jet events and 4 Drell–Yan events. Tree/summary signed-weight
closure was exact for the backgrounds and differed by less than `1e-18`
for the signal. A representative plot was checked visually; this sparse
pilot has negative histogram bins and does not validate shape precision.

This downstream smoke test reused the existing validated signal LHE and
109-point ihixs record. Fresh POWHEG production integration and the complete
109-point calculation were not repeated as part of this entry-point test;
the normal command invokes their existing validated wrappers when those
inputs are not supplied. Earlier failed setup/smoke attempts were retained.

The workflow tests mock expensive external program boundaries and cover
stage ordering, invalid-normalization rejection, preservation of existing
campaigns, pilot checks, private PDF installation and MG5 packaging.
The 14 new tests passed on both the development machine and Timur. The
existing 111 ihixs, HO campaign and normalization tests also passed using
a canonical temporary path on macOS (`TMPDIR=/private/tmp`).
