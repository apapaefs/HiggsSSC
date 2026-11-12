# HiggsSSC

Study a 125 GeV Higgs boson at a 40 TeV proton–proton collider using
MadGraph, Herwig/HwSim and a simplified SSC/GEM detector response.
Start with a small leading-order (LO) sample, inspect it, then move to
higher order (HO) or the four-lepton channel.

## Analyze the prepared four-lepton events on Timur

Students can use the completed **40 TeV HO four-lepton dataset** without
generating events, running MadSpin, or repeating the showers. The instructor
publishes the shared release once using the guide below; `attach` checks that
it is complete and readable before connecting it. From your own
Timur account, clone the repository (or `git pull --ff-only` in an existing
clone), then run:

```bash
git clone https://github.com/apapaefs/HiggsSSC.git
cd HiggsSSC
python3 student_events.py attach
python3 student_events.py setup
python3 student_events.py check
```

`attach` connects your workspace to
`/home/shared/HiggsSSC/fourlepton/ssc40-20261112`. It copies only small, editable
analysis cards into `hfourlepton/student_runs/ssc40/cards/`; the event files
stay in the shared directory. `setup` installs the Python analysis packages in
that workspace's `.venv`. This route needs neither a generator installation
nor access to the instructor's home directory. Do **not** run
`run_pipeline.py setup` or `run_pipeline.py run` for this exercise.

Run the analyses in a [screen session](docs/TERMINAL.md):

```bash
# Cutflow and all 17 observables, comparing perfect and SSC responses.
python3 student_events.py report

# Cut-based analysis and classifier training/evaluation with SSC response.
python3 student_events.py cuts
python3 student_events.py xgboost
```

Results go to `hfourlepton/student_runs/ssc40/results/baseline/`, in
`report_both/`, `cuts_ssc/` and `xgboost_ssc/`. Each contains `index.html` and
machine-readable results; its log and command record are alongside it.
These commands read existing analysis ROOT trees and preserve the signed
weights, source-event grouping, higher-order normalization and detector
response. They cannot launch generation or showering. For this release, a
full SSC analysis reads over 600 million detector-response records; allow
tens of minutes on Timur. `attach` itself does not run an analysis.

Edit your private `cards/cuts_ssc.json` or `cards/xgboost_ssc.json` to change
the cuts or classifier settings. Use a new result name for each variation:

```bash
python3 student_events.py cuts --name tighter_cuts
python3 student_events.py xgboost --name classifier_v2
python3 student_events.py cuts --profile perfect --name perfect_comparison
```

Existing result directories are never overwritten. `--luminosity-fb 100`
is the default. For another private workspace, pass
`--workspace /path/in/your/account` to **every** command. Repeating `attach`
preserves your edited cards. `check --full` additionally hashes every large
event file; normal analysis also applies the existing physics and file audits.

The shared release contains the following nominal samples in each of the
`4e`, `4mu` and `2e2mu` channels:

| Process | Hard events per channel |
| --- | ---: |
| ggF Higgs | 20,000 |
| VBF Higgs | 10,000 |
| W+H, W−H, qqZH, ggZH, ttH (each) | 5,000 |
| qq → 4ℓ continuum | 1,000,000 |
| gg → 4ℓ continuum | 2,000 |

There are **27 samples and 3,171,000 requested hard events**, with both
`perfect` and `ssc` analysis trees. `perfect` means ideal detector response;
`ssc` includes the implemented SSC/GEM resolutions and efficiencies. The
continuum samples cover the 90–180 GeV generation window. Reducible-background
control samples are not part of this nominal set. Shower losses and their
normalization denominators are recorded in the manifests.

The release also provides the LHE inputs and HwSim shower ROOT files under
`nominal/events/` for advanced reconstruction work, original provenance under
`provenance/`, and the instructor's completed plots, cuts and classifier under
`reference/`. The default student commands use the saved reconstruction;
changing an analysis card does not change the detector simulation or its
baseline object selection.

An additional 200,000 qq events per channel are supplied separately through
`supplemental/manifest_ssc.json`. They were used to extend the **frozen**
classifier's test, with results under `reference/extended_test/`. Do not add
them to the nominal samples as another background rate or mix them into
training when reproducing that test. The nominal classifier test alone had
88.8 effective background events; the independent extension reached 160.7.
New student training runs evaluate their own split and statistical adequacy;
they do not automatically reproduce the extended test.

For instructors publishing a new release, see the
[shared-event publication guide](hfourlepton/SHARED_EVENTS.md).

## 1. Clone and install on Timur

Log in to your own Timur account. Run all commands from the repository root:

```bash
git clone https://github.com/apapaefs/HiggsSSC.git
cd HiggsSSC
python3 run_pipeline.py setup --workflow gammagamma-lo
python3 run_pipeline.py check --workflow gammagamma-lo
```

Setup uses Timur's `herwig/stable` module and installs private dependencies
under `inputs/pipeline/`. It does not need an instructor's generated files or
access to the private paper submodule. The shared runtime must provide Herwig,
HwSim, ROOT/PyROOT, LHAPDF and FastJet; setup checks it before building.
Elsewhere, supply `--herwig-env /path/to/Herwig` and compatible compilers;
Timur is the primary tested environment.

## 2. Run a small LO campaign

Use a persistent [screen session](docs/TERMINAL.md) so disconnecting from SSH
does not end your shell.

```bash
python3 run_pipeline.py run --workflow gammagamma-lo --run-tag lo_pilot_01
```

This generates **30 events per sample** for the signal and three backgrounds,
showers them, analyzes the detector response and makes an HTML report.
Defaults are four workers and 100 fb⁻¹. A pilot checks execution and event
bookkeeping; it is too small for a physics conclusion.

```bash
python3 run_pipeline.py status --workflow gammagamma-lo --run-tag lo_pilot_01
```

Open `hgammagamma/LOAnalysis/plots/lo_pilot_01/index.html` using the
[remote viewing instructions](SSHFS_REMOTE_MOUNT.md). Check the event counts,
all four samples, signed weights and plot labels. The underlying files are in
`hgammagamma/LOAnalysis/runs/lo_pilot_01/`.

## 3. Continue or increase statistics

Repeat the **same run command** to resume completed stages. Use a new tag when
changing event counts or scientific settings:

```bash
python3 run_pipeline.py run --workflow gammagamma-lo \
  --run-tag lo_10k_01 --nevents 10000 --jobs 4
```

To remake a report from saved analysis outputs:

```bash
python3 run_pipeline.py report --workflow gammagamma-lo \
  --run-tag lo_10k_01 --luminosity-fb 100
```

Diphoton reports show unit-area shapes and physical cross-section plots side
by side for each observable, with separate PNG, SVG and CSV downloads.

A report cannot replace missing analysis. Read
[recovery instructions](docs/TROUBLESHOOTING.md) before deleting files or
starting another calculation.

## 4. Higher order or four leptons

HO adds POWHEG HJMiNNLO signal generation, NLO backgrounds and an independent
ihixs N3LO signal-rate calculation. Setup extends the existing installation:

```bash
python3 run_pipeline.py setup --workflow gammagamma-ho
python3 run_pipeline.py run --workflow gammagamma-ho --run-tag ho_pilot_01
```

The first HO run includes the full 109-point normalization calculation and
POWHEG integration, even for a 30-event pilot. It can take hours. See
[HO conventions](hgammagamma/HOAnalysis/README.md) and
[normalization](hgammagamma/HOAnalysis/IHIXS_STUDENT_GUIDE.md).

The four-lepton workflow runs nine irreducible samples, reuses their generated
LHE events for separate `perfect` and `ssc` showers/analyses, and compares them:

```bash
python3 run_pipeline.py setup --workflow fourlepton-lo
python3 run_pipeline.py run --workflow fourlepton-lo --run-tag fourl_pilot_01
```

Its report is `hfourlepton/runs/fourl_pilot_01/report/index.html`.
Each mass spectrum pairs a unit-area shape with physical event counts,
separately for each response profile, channel and process category.
See [four-lepton physics and advanced samples](hfourlepton/README.md).

The new `fourlepton-ho` workflow prepares 13 TeV LHC validation and 40 TeV SSC
campaigns, including ggF, VBF, VH and ttH, matched higher-order backgrounds,
ihixs normalization, cut analysis and XGBoost. Start with `prepare`, which
writes inputs without running generators. Follow the
[four-lepton HO run guide](hfourlepton/HIGHER_ORDER.md) for both energy-specific
pilots and the subsequent validation. Numerical LHC validation awaits those runs.

For a smaller complete student campaign, use `--preset minimal`. It includes
all 27 HO samples with process-dependent counts (231,000 events in total):

```bash
python3 run_pipeline.py prepare --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal --jobs 40
python3 run_pipeline.py run --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal --jobs 40
```

First complete setup and the energy-specific pilot described in the
[minimal-preset guide](hfourlepton/HIGHER_ORDER.md#minimal-student-preset).
Add `--event-multiplier 2` with a new tag to double every sample; use
`--energy-tev 40` and a separate tag for the SSC. The 1–2 day target on 40 cores
still needs a full timing run. You launch all simulations yourself.
For 200,000 events in each qq-background channel with the other minimal counts
unchanged, add `--qq-event-multiplier 10` (771,000 events total). The global
event multiplier can still scale this entire set uniformly.

For an existing 13 TeV HO campaign, `reprocess --revision NAME` rebuilds the
CMS selection from saved showers and writes separate report, validation,
cut and XGBoost results. It does not rerun event generation. See the
[CMS reprocessing instructions](hfourlepton/HIGHER_ORDER.md#reprocess-the-cms-selection-using-saved-showers).

## Next steps and references

- [Additional cuts and optional XGBoost](docs/ANALYSIS.md).
- [Diphoton physics](hgammagamma/README.md) and [detector model](hgammagamma/SSC_DETECTOR_RESPONSE.md).
- [Validation evidence and remaining limitations](VALIDATION.md).
- [Terminal basics](docs/TERMINAL.md) and [troubleshooting](docs/TROUBLESHOOTING.md).
- [Historical detector sources](SSCInfo.md) and [study questions](SSCNarrative.md).

`python3 run_pipeline.py --help` needs no simulation libraries. Add `--dry-run`
to inspect an action without installing or running anything. Existing channel
scripts remain available for advanced use; their `--help` lists stage-specific
options. Update source with `git pull --ff-only`; keep your run directories and
manifests. Generated events, dependencies and reports are not included in Git.
