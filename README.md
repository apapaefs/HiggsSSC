# HiggsSSC

Study a 125 GeV Higgs boson at a 40 TeV proton–proton collider using
MadGraph, Herwig/HwSim and a simplified SSC/GEM detector response.
Start with a small leading-order (LO) sample, inspect it, then move to
higher order (HO) or the four-lepton channel.

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
See [four-lepton physics and advanced samples](hfourlepton/README.md).

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
