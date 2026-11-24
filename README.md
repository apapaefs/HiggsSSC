# HiggsSSC

Study a 125 GeV Higgs boson at a 40 TeV proton–proton collider using
MadGraph, Herwig/HwSim and a simplified SSC/GEM detector response.
Start with a small leading-order (LO) sample, inspect it, then move to
higher order (HO) or the four-lepton channel.

## How much SSC luminosity would a Higgs discovery have needed?

At 40 TeV, the combined diphoton and four-lepton median local 5σ reach for a
125 GeV SM Higgs is **1.83 fb⁻¹ (2.20 design months)**, or **2.56 fb⁻¹
(3.07 months)** at the four-lepton heavy-flavour upper bound. Students can
refit the full background campaign and redraw the figures from a clean checkout.

```bash
python3 -m pip install -r discovery/requirements.txt
python3 student_events.py discovery
```

[Results, plots and the comparison with GEM](discovery/README.md) ·
[Inputs, checksums and replication on Timur](discovery/REPLICATION.md)

## Analyze the prepared 40 TeV events on Timur

In your Timur checkout, run `git pull --ff-only` to get the student analyses and
the compact background distributions for both channels. For a first checkout, use
`git clone https://github.com/apapaefs/HiggsSSC.git`, then `cd HiggsSSC`.
The following commands use the prepared events on the shared disk. You can
start at the channel you want.

### Four leptons

```bash
python3 student_events.py attach
python3 student_events.py setup
python3 student_events.py check
python3 student_events.py report
python3 student_events.py cuts
python3 student_events.py xgboost
```

The release has **30 samples**. In each of the `4e`, `4mu` and `2e2mu` channels
it supplies ggF, VBF, W⁺H, W⁻H, qqZH, ggZH and ttH signal and the qq and gg
continuum backgrounds. Those 27 samples comprise 3,171,000 requested hard events.
The three additional samples are Z(→ee) + bb, Z(→μμ) + bb and dileptonic tt̄,
with one million unbiased events each. Both response profiles use the same
generated events. The two control samples with a Z are distinguished by its
prompt decay, while reconstructed events are classified by all four leptons.

The report compares ideal (`perfect`) and SSC/GEM (`ssc`) detector responses and
includes cutflows and 17 observables. Cuts and XGBoost default to the SSC response.
The per-process cutflows show where the reducible backgrounds enter. In the SSC
sample, zero events from these three processes pass the final fiducial selection.
With ideal detector response one Z(→μμ) + bb event passes at 93.5 GeV. Both
profiles therefore give zero simulated heavy-flavour yield in the 120–130 GeV
window. The finite-sample upper bound for the SSC response is
reported separately in the [background study](docs/ANALYSIS.md).

### Diphotons

```bash
python3 student_events.py attach --channel diphoton
python3 student_events.py setup --channel diphoton
python3 student_events.py check --channel diphoton
python3 student_events.py report --channel diphoton
python3 student_events.py cuts --channel diphoton
python3 student_events.py xgboost --channel diphoton
```

The SSC defaults include the ggF Higgs signal and five backgrounds. Alongside
prompt diphotons from quark-initiated production, photon + jet and Drell–Yan,
the plots include the **loop-induced gg → γγ continuum at LO** and **dijets with
two jets misidentified as photons**. The dijet curve uses the upper GEM fake-rate
scenario, treating both jets as quark jets and their fake probabilities as
independent. It enters the total once. The all-gluon alternative and the
assumptions behind the bound are documented with the
[event tables](hgammagamma/HOAnalysis/student_backgrounds/README.md).
Rare prompt-photon processes are outside this sample set. The Higgs signal
uses its existing ihixs normalization.

### Workspaces and results

The shared releases are under `/home/shared/HiggsSSC/fourlepton/ssc40-20261118`
and `/home/shared/HiggsSSC/diphoton/ssc40-20261118`. The Git checkout also
supplies the loop and dijet event tables and both detector responses for the
three heavy-flavour samples. The standard 40 TeV higher-order generation
workflows use these distributions when they make reports or run the default
cuts and classifier. Each channel has a private `student_runs/ssc40-20261118/`
under `hfourlepton/` or `hgammagamma/`. The workspace contains editable cards,
its Python environment and the results. A generator installation is unnecessary.
Use a [screen session](docs/TERMINAL.md) for longer analyses. The four-lepton
SSC report scans over 600 million detector-response records.

Open `results/baseline/report_both/index.html` for the four-lepton report or
`results/baseline/report_ssc/index.html` for diphotons. Both channels write the
cut and classifier reports under `cuts_ssc/` and `xgboost_ssc/` alongside the
report directory. Logs and execution records identify the release and analysis
code for each run. All three commands default to 100 fb⁻¹.

Students can also browse the completed default plots on Timur under
`/home/shared/HiggsSSC/student-results/ssc40-20261118/`. Open the report,
cuts or classifier directory for your channel, or
[copy a report home for viewing](SSHFS_REMOTE_MOUNT.md).

Edit the corresponding files under `cards/` and give each variation a new name.
The diphoton classifier card's `baseline_cuts_config` selects the cut card used
in its comparison. Set it to the path of your private cut card to compare to
your own cuts.

```bash
python3 student_events.py cuts --name tighter_cuts
python3 student_events.py cuts --profile perfect --name perfect_comparison
python3 student_events.py xgboost --channel diphoton --name classifier_v2
```

On upgrade, four-lepton `attach` carries edits from the former default
`hfourlepton/student_runs/ssc40/` cards into the new workspace. The former
results remain accessible there. Repeating `attach` preserves edits in an
already connected workspace. Use `--workspace /path/in/your/account` on each
command to choose a private directory. Use `check --full` to verify every file
against its SHA-256 checksum.

The independent four-lepton qq extension (200,000 events per channel) and the
original instructor classifier results are historical references. Student
classifier runs train and choose a threshold on the default sample set.
The discovery study used frozen classifiers. For the underlying events and
publication procedure, see the [four-lepton release guide](hfourlepton/SHARED_EVENTS.md).

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
