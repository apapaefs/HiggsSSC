# Timur: One HiggsSSC Checkout

Use `/home/apapaefs/Projects/HiggsSSC` for the LO and HO analyses,
ihixs calculations, inputs and reports. All commands below run from
this repository root.

The canonical layout is:

| Item | Path relative to HiggsSSC |
| --- | --- |
| Pinned ihixs source | `external/ihixs/` |
| MG5 3.5.15 runtime | `MG5_aMC_v3_5_15/` |
| Signal LHE | `hgammagamma/HOAnalysis/inputs/powheg-hjminnlo-merged.lhe` |
| Local generation PDFs | `hgammagamma/HOAnalysis/inputs/lhapdf/` |
| ihixs builds and integrations | `hgammagamma/HOAnalysis/normalization/ihixs-ssc40/` |
| Validated N3LO rate | `hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json` |
| New production campaign | `hgammagamma/HOAnalysis/runs/ho_100k_02/` |
| Production report | `hgammagamma/HOAnalysis/plots/ho_100k_02/` |

The complete 100,000-event signal LHE and the full
`NNPDF40_nlo_as_01180_qed` and `NNPDF40_nnlo_as_01180_qed` sets were copied
into `HOAnalysis/inputs/` and verified against their source files. The
missing MG5 runtime files were copied into the existing repository-local
installation, preserving its custom model files and updating `mg5_path`.
The ihixs work directory and existing LO source files were preserved.
Inputs, installed generators, builds and generated results are ignored
by Git; the scripts, cards and documentation share one GitHub repository.

The signal LHE has SHA-256:

```text
01807025f5e6feb21d7d8f836f3161c4aeadce5c44af81edc305b3256ce51187
```

The previous prepared production campaign `ho_100k_01` and completed pilot
`ho_timur_smoke` were copied intact to
`HOAnalysis/runs/archive-before-consolidation/`. Their source directories
remain as backups. These archives retain historical absolute paths in
cards, manifests and generator products, so do not resume them at the
copied location. Use the fresh tags below. The prepared production
backgrounds had no completed stages and no generated events; the full
signal LHE is reused. Copy inventories and the MG5 configuration change
are recorded in the archive's `CONSOLIDATION.json`.

Earlier checks used Python 3.9.25, GCC 11.5, Herwig 7.3.0, ThePEG 2.3.0,
LHAPDF 6.5.3, ROOT 6.40.04 and the installed HwSim plugin. The old
30-event pilot preserved three negative-weight events and response-weight
closure to `8.2e-20`; the signal LHE passed the 100,000-event validation.
Those checks predate ihixs normalization and consolidation. The new
workflow, generator copy and regression tests are left for you to run.

## After ihixs Finishes

First complete the build, benchmark and all 109 integrations described in
[`README.md`](README.md#calculate-the-ho-signal-normalization). The
validated N3LO record above is mandatory for HO signal analysis. The
existing NNLO signal input and NLO backgrounds retain their generation PDFs.

Log in to Timur, preferably inside `screen` or `tmux`, and set up:

```bash
cd /home/apapaefs/Projects/HiggsSSC
git pull --ff-only origin main
source /etc/profile.d/modules.sh
module load herwig/stable

export LHAPDF_DATA_PATH="$PWD/hgammagamma/HOAnalysis/inputs/lhapdf:/home/shared/Herwig/share/LHAPDF"
SIGNAL_NORM="$PWD/hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json"
SIGNAL_LHE="$PWD/hgammagamma/HOAnalysis/inputs/powheg-hjminnlo-merged.lhe"
LUMINOSITY_FB=100  # set your intended luminosity in inverse femtobarns

HO_OPTIONS=(
  --nb-core 8
  --herwig-module herwig/stable --herwig Herwig
  --mg5-dir "$PWD/MG5_aMC_v3_5_15" --mg5-python python3
  --mg5-fortran /usr/bin/gfortran --mg5-cxx /usr/bin/g++
  --analysis-cxx /usr/bin/g++
  --signal-normalization "$SIGNAL_NORM"
  --signal-lhe "$SIGNAL_LHE"
)
```

Optionally run the Python regression suites yourself before the campaign:

```bash
python3 -m unittest discover -s tests -p 'test_ihixs*.py' -v
python3 -m unittest discover -s tests -p 'test_ho_normalization_*.py' -v
python3 -m unittest discover -s tests -p 'test_gammagamma*.py' -v
```

Run a fresh 30-event pilot, then inspect its report:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage all --run-tag ho_n3lo_pilot_01 --nevents 30 \
  "${HO_OPTIONS[@]}" --resume

python3 hgammagamma/make_gammagamma_report.py \
  --analysis-root hgammagamma/HOAnalysis/runs/ho_n3lo_pilot_01 \
  --run-tag ho_n3lo_pilot_01 \
  --output-dir hgammagamma/HOAnalysis/plots/ho_n3lo_pilot_01 \
  --no-density --normalization event_xsec --luminosity-fb "$LUMINOSITY_FB"
```

This includes the first integration/shower validation for the backgrounds.
NLO integration still takes time when requesting only 30 events. Check
the event counts, signed-weight diagnostics, branching-ratio normalization
and report before launching production.

Run 100,000 events per sample with the same options:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage generate --run-tag ho_100k_02 --nevents 100000 \
  "${HO_OPTIONS[@]}" --resume

python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage shower --run-tag ho_100k_02 --nevents 100000 \
  "${HO_OPTIONS[@]}" --resume

python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage analyze --run-tag ho_100k_02 --nevents 100000 \
  "${HO_OPTIONS[@]}" --resume

python3 hgammagamma/make_gammagamma_report.py \
  --analysis-root hgammagamma/HOAnalysis/runs/ho_100k_02 \
  --run-tag ho_100k_02 \
  --output-dir hgammagamma/HOAnalysis/plots/ho_100k_02 \
  --no-density --normalization event_xsec --luminosity-fb "$LUMINOSITY_FB"
```

This showers 100,000 existing signal events and generates 100,000 MC@NLO
events for each of prompt gamma-gamma, gamma+jet and Drell--Yan before
showering and applying the SSC/GEM response. The event count is before
detector selection. Samples run sequentially; MG5 uses up to eight cores.
Run the commands in order, waiting for each to finish successfully.
The `generate` stage builds background processes as needed; `analyze`
builds the analysis executable. The signal LHE already exists, so
`generate` validates it rather than generating another signal sample.
The runner initializes the module system in its subprocesses.
The explicit C++ compiler removes the module's `-std=c++14` suffix from
MG5's compiler executable setting.

Repeat the campaign command with identical options and `--resume` after
an interruption. Completed stages with matching configuration are reused;
partial Herwig ROOT products are preserved for inspection. Changed
configuration requires a new tag. The report is
`HOAnalysis/plots/ho_100k_02/index.html`; physical rates and expected yields
are in `data/sample_summary.csv` within that report directory. Only the
HO Higgs signal uses the ihixs N3LO rate, with the physical diphoton BR
applied once. Yields use your specified luminosity.

The physics limitations and earlier validation status are in
[`README.md`](README.md) and [`VALIDATION.md`](VALIDATION.md). No builds,
tests, integrations or analyses were run while consolidating the checkout.

## Recover From The Ninja Installer Error

The first gamma+jet export on Timur failed while downloading
`HEPToolsInstaller_V168.tar.gz`: the incomplete archive left no
`HEPToolInstaller.py`, so the automatic Ninja installation failed.
There was no exported gamma+jet `mg5_process` directory. The signal and
prompt diphoton samples had already completed `generate` and were preserved.

Timur already has a Ninja installation with MG5's required library and
Fortran-module layout in the Herwig stack. The repository-local
`MG5_aMC_v3_5_15/input/mg5_configuration.txt` was backed up and configured
with these settings:

```text
ninja = /home/shared/Herwig/opt/MG5_aMC_v3_5_1/HEPTools/ninja/lib
collier = None
```

The Ninja library directory contains `libninja.a`, the OneLOop library,
and the expected Fortran module/header files in `../include/`. Selecting
this existing library avoids the failed automatic download. The bundled
CutTools and IREGI sources remain available. This is a site-local MG5
runtime setting; the generator installation and configuration backups
are ignored by Git. No campaign cards, fingerprints or completed LHE
files were changed by the repair.

After applying the common shell setup above, retry the same production
generation command:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage generate --run-tag ho_100k_02 --nevents 100000 \
  "${HO_OPTIONS[@]}" --resume
```

This reuses the completed signal and prompt diphoton stages, then retries
gamma+jet and continues to Drell--Yan. Once it succeeds, run `shower`,
`analyze` and the report command above with the same options. There is
no need for a new run tag or a repeat ihixs calculation. Library files
and configuration were inspected; no dependency build, generator,
analysis or tests were run while preparing this repair.

## Recover A Finalized Shower At End Of The LHE File

The signal shower in `ho_100k_02` consumed all 100,000 LHE records and
saved 99,994 events. Six `D_s0(2590)+/-` decays failed, so Herwig discarded
those events and reached the end of the finite file while trying to fill
its target of 100,000 successful events. It reported
`More events requested than available in LesHouchesReader` after HwSim
had finalized the ROOT file. Metadata inspection confirmed that the
latest `Data` tree has 99,994 entries and the file is neither corrupt nor
marked as recovered by ROOT.

The runner now audits this specific termination with `--resume`. It
requires the unchanged complete LHE, final Herwig/HwSim statistics, the
single expected EOF run error, matching discard counts and readable
closed ROOT files. It saves separate source/attempted/generated/ROOT
counts and hashes before marking the shower stage complete. An arbitrary
truncated file or unrelated failure is still refused. This preserves the
existing event file, original requested count and generation fingerprint.
The card audit counts the two `set ...:Cuts` assignments. The separate
`create ThePEG::Cuts /Herwig/Cuts/NoCuts` declaration creates the object
and is not a third cut assignment.

In the same shell with the common `HO_OPTIONS` above, run:

```bash
cd /home/apapaefs/Projects/HiggsSSC
git pull --ff-only origin main

python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage shower --run-tag ho_100k_02 --nevents 100000 \
  "${HO_OPTIONS[@]}" --resume
```

The signal is audited and adopted without rerunning it, then the three
background showers proceed. Do not delete the signal ROOT file or change
`--nevents` to the saved count. After this command succeeds, use the
`analyze` and report commands above. The ihixs normalization record is
required before analysis; no new ihixs calculation is needed.

For new showers, yields use the signed sum of all consumed LHE weights
in Herwig's weight convention, including records discarded during
showering. Their detector response is unavailable and is recorded as
zero simulated response; the surviving sample is not rescaled by an
unweighted event-count fraction. The report separates saved weights,
source normalization and shower-quality diagnostics. The existing signal
log also records momentum-consistency warnings, with a maximum violation
of about 4.12 TeV, independently of the six failed decays. Recovery
preserves these diagnostics; it does not correct those simulation issues.

The recovery and signed-weight tests are provided for you to run:

```bash
python3 -m unittest discover -s tests -p 'test_ho_shower_completion.py' -v
python3 -m unittest discover -s tests -p 'test_gammagamma_ho_campaign.py' -v
python3 -m unittest discover -s tests -p 'test_ho_normalization_*.py' -v
```

No showers, analysis, builds or tests were run while preparing this repair.

## Resume ihixs After A Production Precision Failure

The `ihixs-ssc40-refined01` build and benchmark passed, but `scale_2`
returned raw `eftn3lo = 226.524792 +/- 0.121277162 pb`: a 0.053538%
numerical error, above the required 0.05%. About 97% of its numerical
variance came from the `q1q2` contribution. The pinned ihixs source
relaxes that N3LO term's relative target by a factor of 10,000. Its
Cuhre integrator also uses a fixed minimum of 1,000 evaluations, so the
increased Vegas statistics did not force more evaluations of that term.

Use the same refined settings, workspace and installed builds. The new
`--refine-failed` option reuses passing points and retries only inaccurate
ones with tighter relative targets in separate directories. For this
profile, the first retry uses `epsrel = 1e-6`; at most two further retries
use `1e-7` and `1e-8`. The 0.05% final-error requirement remains in force.
The settings file, plan, original results and ROOT event samples are kept.

In the shell where `IHIXS_REFINED` is defined, run:

```bash
cd /home/apapaefs/Projects/HiggsSSC
git pull --ff-only origin main
python3 hgammagamma/run_ihixs_normalization.py \
  --stage calculate "${IHIXS_REFINED[@]}" --resume --refine-failed
```

If using a fresh shell, reconstruct the options explicitly:

```bash
source /etc/profile.d/modules.sh
module load herwig/stable
export LHAPDF_DATA_PATH="$PWD/hgammagamma/HOAnalysis/inputs/lhapdf:/home/shared/Herwig/share/LHAPDF"
python3 hgammagamma/run_ihixs_normalization.py \
  --stage calculate --resume --refine-failed \
  --settings hgammagamma/HOAnalysis/normalization/ihixs-ssc40-refined01-settings.json \
  --work-dir hgammagamma/HOAnalysis/normalization/ihixs-ssc40-refined01 \
  --herwig-module herwig/stable --lhapdf-dir /home/shared/Herwig \
  --cuba-dir "$HOME/.local/cuba-4.2.2-gcc11" --cc cc --cxx c++ --jobs 8
```

There is no need to repeat `build` or `benchmark`. A successful calculation
must still finish all 109 points before creating `ggf-ssc40-n3lo.json`.
Then use the analysis and report commands above. The mocked regression
tests are in `tests/test_ihixs_normalization.py` for you to run; no
calculations, builds or tests were run while preparing this update.

## Recover ihixs After A Parser Heap Abort

The resumed `ihixs-ssc40-refined01` calculation completed all seven scale
points and 100 PDF replicas, including two precision refinements. Its
`nnlo_native` attempt printed the NNLO result and wrote `ihixs.out`, then
aborted with return code `-6` (`SIGABRT`) and
`corrupted size vs. prev_size`. The two native-PDF comparison points are
still needed. The aborted output cannot be adopted as a completed run.

Static inspection found two bounds errors in the pinned option parser:
its option-name allocation omits space for the terminating NUL byte, and
its terminating option entry is written past the allocated array. These
are a possible cause of the heap abort. `--repair-parser` applies only
those two corrections in separate source copies. The physics cards,
PDFs, hard alpha_s and integration settings stay the same. This option
recovers an existing work directory with an original completed build
manifest; the default build and pinned submodule remain unchanged.

Keep the original refined settings and work directory. The repair adds
`parser-repair/` inside that directory, with separate `source-upstream`,
`source-lhapdf`, `build-upstream`, `build-lhapdf`, `build-manifest.json`,
`benchmark.json` and `runs/benchmark/`. The original builds, benchmark and
107 verified completed points remain intact. With `--resume`, the failed
`nnlo_native` directory is archived before retrying; unfinished and new
points use the repaired production executable. The final record retains
both benchmarks and the hashes for the original and repaired builds and
selected production attempts. `calculate` also runs a repaired 40 TeV
central check in `parser-repair/runs/central_check/` (or a separate
`central_check__precision_N/` directory), requiring a numerical
error of at most 0.05% and agreement with the cached central result within
that same relative tolerance. Both central results are recorded; a
mismatch stops publication. Thus this recovery reuses the 107 completed
points and runs one central verification, the two unfinished comparisons
and a fresh 13 TeV benchmark.

Run the following yourself on Timur. The loop stops at the first failure
without closing your shell. It builds and benchmarks the separate repair
before resuming the outstanding calculations:

```bash
cd /home/apapaefs/Projects/HiggsSSC
git pull --ff-only origin main
source /etc/profile.d/modules.sh
module load herwig/stable
export LHAPDF_DATA_PATH="$PWD/hgammagamma/HOAnalysis/inputs/lhapdf:/home/shared/Herwig/share/LHAPDF"

for ihixs_stage in build benchmark calculate; do
  IHIXS_REPAIR_STAGE=()
  if [ "$ihixs_stage" = calculate ]; then
    IHIXS_REPAIR_STAGE=(--refine-failed)
  fi
  python3 hgammagamma/run_ihixs_normalization.py \
    --stage "$ihixs_stage" --repair-parser --resume "${IHIXS_REPAIR_STAGE[@]}" \
    --settings hgammagamma/HOAnalysis/normalization/ihixs-ssc40-refined01-settings.json \
    --work-dir hgammagamma/HOAnalysis/normalization/ihixs-ssc40-refined01 \
    --herwig-module herwig/stable --lhapdf-dir /home/shared/Herwig \
    --cuba-dir "$HOME/.local/cuba-4.2.2-gcc11" --cc cc --cxx c++ --jobs 8 || break
done
```

The final `ggf-ssc40-n3lo.json` is written only after all 109 points,
both benchmarks, the central comparison and the unchanged 0.05% error
requirement pass. Preserve the original work directory and its
`parser-repair/` tree with that record.
This repair remains unverified until you run the build, benchmark and
calculation. No builds, tests, integrations or analyses were run while
preparing it.

The ihixs regression pattern includes the original normalization tests and
the new parser-repair and record-provenance tests. Run it yourself:

```bash
python3 -m unittest discover -s tests -p 'test_ihixs*.py' -v
```
