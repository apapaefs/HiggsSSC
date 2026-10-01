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
python3 -m unittest discover -s tests -p 'test_ihixs_normalization.py' -v
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
  --stage all --run-tag ho_100k_02 --nevents 100000 \
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
The `all` stage builds the analysis executable and background processes
as needed. The runner initializes the module system in its subprocesses.
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
