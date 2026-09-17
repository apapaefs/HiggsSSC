# Timur: 100,000 events per sample

Use `/home/apapaefs/Projects/HiggsSSC-HO` for this higher-order campaign.
This clean checkout preserves the uncommitted work in the older
`/home/apapaefs/HiggsSSC` and `/home/apapaefs/Projects/HiggsSSC` checkouts.

The separate input directory `/home/apapaefs/Projects/HiggsSSC-HO-inputs`
contains the complete 100,000-event HJMiNNLO signal LHE and the full
`NNPDF40_nlo_as_01180_qed` and `NNPDF40_nnlo_as_01180_qed` PDF sets. They
were copied from the locally validated inputs and are not committed to Git.
The signal file has SHA-256:

```text
01807025f5e6feb21d7d8f836f3161c4aeadce5c44af81edc305b3256ce51187
```

The synchronized setup was checked on Timur with Python 3.9.25, GCC 11.5,
Herwig 7.3.0, ThePEG 2.3.0, LHAPDF 6.5.3, ROOT 6.40.04, and the installed
HwSim plugin. All 53 gamma-gamma regression tests passed. A separate
30-event signal pilot completed showering and detector analysis, preserving
three negative-weight events and response-weight closure to `8.2e-20`.
The full signal input passed the 100,000-event LHE validation, and cards
for all four production samples were prepared without generating events.

Run this after logging in to Timur, preferably inside `screen` or `tmux`:

```bash
cd /home/apapaefs/Projects/HiggsSSC-HO
export LHAPDF_DATA_PATH=/home/apapaefs/Projects/HiggsSSC-HO-inputs/lhapdf:/home/shared/Herwig/share/LHAPDF
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage all --run-tag ho_100k_01 --nevents 100000 --nb-core 8 \
  --herwig-module herwig/stable \
  --mg5-dir /home/apapaefs/HiggsSSC/MG5_aMC_v3_5_15 \
  --mg5-fortran /usr/bin/gfortran --mg5-cxx /usr/bin/g++ \
  --signal-lhe /home/apapaefs/Projects/HiggsSSC-HO-inputs/powheg-hjminnlo-merged.lhe \
  --resume
```

This showers 100,000 existing signal events and generates 100,000 MC@NLO
events for each of prompt gamma-gamma, gamma+jet, and Drell--Yan before
showering and applying the SSC/GEM response. The event count is before
detector selection. Samples run sequentially; MG5 uses up to eight cores.
The runner initializes the module system itself, including in a non-login
shell. The explicit C++ compiler removes the module's `-std=c++14` suffix
from MG5's compiler executable setting.

Outputs are under `hgammagamma/HOAnalysis/runs/ho_100k_01`. The `--resume`
option reuses complete stages with matching configuration. A partial Herwig
ROOT output requires inspection; it is never silently overwritten.

To create the report afterwards:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --analysis-root hgammagamma/HOAnalysis/runs/ho_100k_01 \
  --run-tag ho_100k_01 \
  --output-dir hgammagamma/HOAnalysis/plots/ho_100k_01 \
  --no-density --normalization event_xsec
```

The local validation status and physics limitations are in
[`VALIDATION.md`](VALIDATION.md) and [`README.md`](README.md). In particular,
the photon-background processes have been exported but still require their
first integration/shower validation. Synchronizing this checkout and its
inputs does not launch the 100,000-event campaign.
