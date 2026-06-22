# Higgs To Gamma Gamma Study

This directory is for the Higgs boson decay study

```text
h -> gamma gamma
```

using the `loop_sm_haa` MadGraph5_aMC@NLO model in this repository.

## Get The Repository

Start by cloning the repository onto the machine where you will run the
campaign.  If your GitHub SSH key is set up, use:

```bash
git clone git@github.com:apapaefs/HiggsSSC.git
cd HiggsSSC
```

If SSH is not set up, use the HTTPS URL instead:

```bash
git clone https://github.com/apapaefs/HiggsSSC.git
cd HiggsSSC
```

All commands below assume that you are running from this repository root.

## Physics Target

The goal is to generate and study Higgs boson decays to two photons, using an
effective Higgs-photon-photon coupling added to the Standard Model loop model.

In MG5 notation the process is

```text
h > a a
```

where `h` is the Higgs boson and `a` is the photon.

## Download MadGraph

From the repository root, download MadGraph5_aMC@NLO version `v3.5.15` from

```text
https://launchpad.net/mg5amcnlo
```

For example:

```bash
wget https://launchpad.net/mg5amcnlo/3.0/3.7.x/+download/LTS_MG5aMC_v3.5.15.tgz
```

This downloads the tarball

```text
HiggsSSC/LTS_MG5aMC_v3.5.15.tgz
```

Unpack it with:

```bash
tar -xzf LTS_MG5aMC_v3.5.15.tgz
```

You should now have:

```text
HiggsSSC/MG5_aMC_v3_5_15
```

The repository already contains the custom model at:

```text
MG5_aMC_v3_5_15/models/loop_sm_haa
```

If you use a separate MG5 installation outside this repository, copy the model
there and pass that MG5 path to the campaign runner:

```bash
cp -r MG5_aMC_v3_5_15/models/loop_sm_haa /path/to/MG5_aMC_v3_5_15/models/
python3 hgammagamma/run_gammagamma_campaign.py \
  --mg5-dir /path/to/MG5_aMC_v3_5_15 \
  --dry-run \
  --nevents 10
```

## Model

Use the model

```text
loop_sm_haa
```

from

```text
MG5_aMC_v3_5_15/models/loop_sm_haa
```

This model contains the effective `h a a` interaction copied from
`loop_sm_scalar`.

## Minimal MG5 Check

From the directory containing the repository clone:

```bash
cd HiggsSSC
module use "$PWD/modulefiles"
module load python/311
module load higgsssc/mg5
mg5_aMC
```

The `higgsssc/mg5` module is clone-local: it uses the `HiggsSSC` checkout from
which you ran `module use "$PWD/modulefiles"`.  This means each student can
clone the repository into their own home directory and get their own `MG5_DIR`.

If you want to run MG5 directly from the MG5 directory instead of using the
wrapper, use:

```bash
cd MG5_aMC_v3_5_15
python3.11 ./bin/mg5_aMC
```

Inside MG5:

```text
import model loop_sm_haa
display interactions h a a
generate h > a a
```

The process should generate one diagram.

## Starter Process Card

A minimal process card for decay-only generation is:

```text
import model loop_sm_haa
generate h > a a
output hgammagamma_decay
```

For a production-and-decay study, generate the Higgs production process first
and let MadSpin handle the decay. For gluon-fusion Higgs production use:

```text
import model loop_sm_haa
generate g g > h [noborn=QCD]
output gg_h_gammagamma
```

The `[noborn=QCD]` option tells MG5 to use the loop-induced gluon-fusion
process rather than looking for a QCD Born contribution.

If MG5 asks to install recommended libraries or helper packages while setting
up this loop-induced process, accept the recommended installs.  Let MG5 manage
those dependencies unless you already have a known working site installation.

## Beam Energy

For the SSC setup we are using 20 TeV proton beams colliding with 20 TeV proton
beams, for a total centre-of-mass energy of 40 TeV. MG5 stores beam energies in
GeV, so each beam should be set to `20000`.

After generating and outputting the process, launch it:

```text
launch gg_h_gammagamma
```

At the launch prompt, set the beam energies before running:

```text
set ebeam1 20000
set ebeam2 20000
```

You can also edit `gg_h_gammagamma/Cards/run_card.dat` directly and set:

```text
     1        = lpp1    ! beam 1 type
     1        = lpp2    ! beam 2 type
 20000.0      = ebeam1  ! beam 1 total energy in GeV
 20000.0      = ebeam2  ! beam 2 total energy in GeV
```

## MadSpin Decay

Use MadSpin to decay the produced Higgs boson:

```text
h > a a
```

After `output gg_h_gammagamma`, edit
`gg_h_gammagamma/Cards/madspin_card.dat`. If the file does not exist yet, copy
the default card first:

```bash
cp gg_h_gammagamma/Cards/madspin_card_default.dat gg_h_gammagamma/Cards/madspin_card.dat
```

Use this minimal MadSpin card:

```text
set spinmode none

decay h > a a

launch
```

To enable MadSpin when launching the run, type:

```text
launch gg_h_gammagamma
madspin=ON
set ebeam1 20000
set ebeam2 20000
```

Then continue the launch as usual. MG5 will generate the hard process
`g g > h [noborn=QCD]` and MadSpin will decay the Higgs according to the
`decay h > a a` line in `madspin_card.dat`.

## Notes To Record

When adding results to this directory, record:

- the MG5 version;
- the exact process card;
- the exact MadSpin card;
- the Higgs mass and width used in the param card;
- `ebeam1 = 20000` and `ebeam2 = 20000` for the 40 TeV setup;
- any generation cuts;
- the number of generated events;
- the random seed;
- whether events are parton-level only or showered/hadronized.

For final normalization, do not use this simple MG5+MadSpin setup as the only
source of the inclusive rate or branching ratio. We will calculate the total
cross section with other tools, for example `iHixs 2`, and take the
`h -> gamma gamma` branching ratio from the CERN Yellow Report branching-ratio
page:

```text
https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
```

We will also use higher-order Monte Carlo generators later as additional
inputs and cross-checks.

## Suggested Outputs

Useful first plots or checks include:

- photon transverse momentum;
- photon pseudorapidity;
- diphoton invariant mass;
- photon separation;

## LO Campaign Pipeline

This section describes the current automated leading-order workflow.  Before
starting, follow the repository and MG5 setup steps above.  All commands below
assume that you are running from the `HiggsSSC` repository root.

The main scripts are:

```text
hgammagamma/run_gammagamma_campaign.py
hgammagamma/make_gammagamma_report.py
```

The first script generates and analyzes samples.  The second script combines
the analysis outputs into stacked plots and a small HTML webpage.

### What The Campaign Does

For each enabled sample, the campaign script:

1. writes an MG5 process card;
2. generates an MG5 process directory;
3. sets the beam energies to 20 TeV per beam;
4. sets the number of events, random seed, and basic generation cuts;
5. runs MG5 event generation;
6. runs MadSpin if the sample has a decay, for example `h > a a`;
7. finds the generated LHE file;
8. writes a Herwig input file from `HW-template.in`;
9. runs `Herwig read` and `Herwig run`;
10. runs `HwSimPostAnalysis_gammagamma`;
11. writes `.top`, `.dat`, `.evp`, and `_var.root` analysis outputs.

At this stage there is no smearing and no internal fake-photon construction.
Reducible backgrounds can be generated as samples, but jet-to-photon and
lepton-to-photon fake shapes still need the later fake-photon analysis code.

### Setup Checklist

You need:

- MG5_aMC installed and available at `MG5_aMC_v3_5_15` in the repository, or
  another path supplied with `--mg5-dir`;
- the `loop_sm_haa` model inside the MG5 `models/` directory;
- Herwig and HwSim available through the Herwig environment;
- ROOT available for the post-analysis executable;
- Python 3 with `matplotlib` for the report script.

On `timur.kennesaw.edu`, which is the Red Hat development machine for this
workflow, Herwig is provided by the module:

```bash
module load herwig/stable
```

The Python runner will load this module automatically on Linux when no
`--herwig-env` activation script is supplied.  To be explicit, you can also
pass:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --herwig-module herwig/stable
```

On the local development machine used to build this workflow, the Herwig
environment is found automatically at:

```text
~/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate
```

On another machine with a Herwig activation script rather than modules, pass
the Herwig environment explicitly if needed:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --herwig-env /path/to/herwig/bin/activate
```

If MG5 is not inside this repository, pass its path:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --mg5-dir /path/to/MG5_aMC_v3_5_15
```

### First Dry Run

Before generating events, always do a dry run.  This prints the commands and
paths without running MG5 or Herwig:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --dry-run \
  --nevents 10 \
  --run-samples signal_gg_h_aa
```

To check both default samples:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --dry-run \
  --nevents 10
```

The default enabled samples are:

```text
signal_gg_h_aa  : g g > h [noborn=QCD], then h > a a
bkg_prompt_aa   : p p > a a
```

### Run A Tiny Smoke Test

Start with 10 events.  This is only a technical test, not a physics result:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --nevents 10 \
  --run-samples signal_gg_h_aa
```

Then test the prompt diphoton background:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --nevents 10 \
  --run-samples bkg_prompt_aa
```

If both commands finish, the full chain is working.

### Run The Default Campaign

To run both the Higgs signal and the prompt diphoton background:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --nevents 10000
```

The default beam energy is already:

```text
ebeam1 = 20000 GeV
ebeam2 = 20000 GeV
```

The default run tag is:

```text
run_01
```

You can change these:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --nevents 50000 \
  --run-tag run_02 \
  --ebeam 20000
```

### Where Outputs Go

Signal samples are written under:

```text
hgammagamma/LOAnalysis/Signal/events/<sample>/
```

Background samples are written under:

```text
hgammagamma/LOAnalysis/Backgrounds/events/<sample>/
```

For example, after running `signal_gg_h_aa`, look for:

```text
hgammagamma/LOAnalysis/Signal/events/signal_gg_h_aa/
```

Important outputs include:

```text
mg5_process/Events/<run>/unweighted_events.lhe.gz
herwig/events/<sample>.root
<sample>_hwsim_roots.input
<sample>_hwsim_roots-<run_tag>.top
<sample>_hwsim_roots-<run_tag>.dat
<sample>_hwsim_roots-<run_tag>.evp
<sample>_hwsim_roots-<run_tag>_var.root
```

The `.top` files contain the histograms.  The `.dat` files contain the analysis
summary, including the number of events read and the number of events with two
selected photons.  The `_var.root` files contain the analysis tree.

### Gamma Gamma Variables

The `_var.root` tree stores:

```text
variables[0] = m_gg
variables[1] = pt_gamma1
variables[2] = eta_gamma1
variables[3] = pt_gamma2
variables[4] = eta_gamma2
variables[5] = deltaR_gg
variables[6] = deltaPhi_gg
variables[7] = pt_gg
variables[8] = y_gg
variables[9] = n_selected_photons
eventweight[0] = evweight * weight_scale
```

All transverse momenta and invariant masses are in GeV.

### Editing Samples

Open:

```text
hgammagamma/run_gammagamma_campaign.py
```

Near the top there is a list called `SAMPLES`.  Each entry has the form:

```python
Sample(name, category, model, process, madspin_decay, weight_scale)
```

The default entries are:

```python
YR4_BR_H_TO_GAMMAGAMMA = 2.27e-3
SIGNAL_GGH_K_FACTOR = 2.0
SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT = SIGNAL_GGH_K_FACTOR * YR4_BR_H_TO_GAMMAGAMMA

Sample("signal_gg_h_aa", "Signal", "loop_sm_haa",
       "g g > h [noborn=QCD]", "h > a a",
       SIGNAL_GGH_TO_GAMMAGAMMA_WEIGHT)

Sample("bkg_prompt_aa", "Backgrounds", "sm",
       "p p > a a", "", 1.0)
```

For the signal, MadSpin decays every generated Higgs boson to photons.  The MG5
cross section is therefore treated as the inclusive LO `g g > h` production
rate, and the default signal weight applies:

```text
weight_scale = K_ggH * BR(H -> gamma gamma)
             = 2.0 * 2.27e-3
             = 4.54e-3
```

The branching ratio is the Standard Model value at `mH = 125.09 GeV` used in
the LHC Higgs Cross Section Working Group
[Yellow Report 4](https://arxiv.org/abs/1610.07922) tables.  The `K_ggH`
factor is a first-pass total cross-section correction; replace it later with a
more precise prediction if needed.

There are also commented entries for:

```text
g g > a a [noborn=QCD]   loop-induced prompt continuum
p p > a j                photon + jet reducible background
p p > j j                dijet reducible background
p p > e+ e-              electron fake background
```

To enable a commented sample, remove the leading `#` and rerun the campaign.

The last number, `weight_scale`, multiplies the event weights in the
post-analysis.  For now this is the place to apply simple total normalization
factors, fake-rate factors, or cross-section rescalings sample by sample.

If you already produced signal `.top` and `.dat` files with `weight_scale = 1`,
rerun the post-analysis or rerun the campaign so the `.dat` file records the
new signal weight.

### Making The HTML Plot Report

After the samples have been analyzed, build the report:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --run-tag run_01
```

The report is written to:

```text
hgammagamma/LOAnalysis/plots/gammagamma_run_01/index.html
```

Open this file in a browser.  If the report was produced on
`timur.kennesaw.edu` and you want to view it from a local Linux machine, see
[Viewing Remote Results With SSHFS](../SSHFS_REMOTE_MOUNT.md).

The report contains:

- stacked histograms;
- backgrounds first, Higgs signal stacked on top;
- one PNG and one SVG per plot;
- one CSV per plot;
- a `sample_summary.csv` file;
- a zip file containing the plot images.

The plots use a compact `SSCwf?` publication style inspired by the
`H -> gamma gamma` figures in the ATLAS discovery paper,
<https://arxiv.org/pdf/1207.7214>: white canvas, black axes with inward ticks,
stepped stacked histograms, and in-plot simulation labels.

By default, the report normalizes each sample to:

```text
cross section * weight_scale * analysis efficiency
```

where the cross section is read from the MG5 banner, and `weight_scale` and the
analysis efficiency are read from the `.dat` file:

```text
analysis efficiency = sum_diphoton_weight / sum_weight
```

The report divides by bin width by default, so a mass distribution such as
`m_gg` is plotted in:

```text
pb / GeV
```

To make non-density plots instead, use:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --run-tag run_01 \
  --no-density
```

To make a report from selected samples only:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --run-tag run_01 \
  --samples signal_gg_h_aa,bkg_prompt_aa
```

### Running Cut And XGBoost Analyses

The repo-root analysis CLI reads the `_var.root` files written by the LO
campaign and produces a small analysis report above the campaign outputs.  It
uses the same MG5 cross sections and `.dat` weight scales used by the plot
report.

Install the Python analysis dependencies in the environment where PyROOT is
available:

```bash
python3 -m pip install -r requirements-analysis.txt
```

PyROOT itself usually comes from the local ROOT installation rather than from
`pip`.

For a rectangular cut analysis, create a YAML card such as:

```yaml
analysis:
  name: baseline_cuts
  run_tag: run_01
  luminosity_fb: 100.0
  cuts:
    - variable: n_selected_photons
      min: 2
    - variable: m_gg
      min: 120.0
      max: 130.0
```

Then run, from the repository root:

```bash
python3 analyze_lo_varfiles.py cuts --config baseline_cuts.yaml
```

Cuts are inclusive and combined with logical AND.  The allowed variable names
are:

```text
m_gg, pt_gamma1, eta_gamma1, pt_gamma2, eta_gamma2,
deltaR_gg, deltaPhi_gg, pt_gg, y_gg, n_selected_photons
```

By default the output is written to:

```text
hgammagamma/LOAnalysis/analyses/<run_tag>/<analysis_name>/
```

The output directory contains:

- `summary.csv`;
- `summary.json`;
- `index.html`.

The summary includes the number of selected Monte Carlo events, the analysis
efficiency, the selected cross section, and the expected event yield at the
luminosity in the YAML card.  The expected event count is calculated as:

```text
expected events = selected cross section [pb] * luminosity [fb^-1] * 1000
```

To analyze only selected samples, add:

```yaml
  samples:
    - signal_gg_h_aa
    - bkg_prompt_aa
```

To write somewhere else, add:

```yaml
  output_dir: /path/to/output
```

The XGBoost mode uses the same sample discovery and normalization, but trains a
binary signal-versus-background classifier and chooses a score threshold that
maximizes the expected significance on the test split.  A minimal card is:

```yaml
analysis:
  name: xgboost_baseline
  run_tag: run_01
  luminosity_fb: 100.0
  xgboost:
    test_size: 0.35
    seed: 12345
    max_events: 1000
```

Run it with:

```bash
python3 analyze_lo_varfiles.py xgboost --config xgboost_baseline.yaml
```

This writes the same `summary.csv`, `summary.json`, and `index.html` files,
plus XGBoost outputs such as `metrics.json`, `scores.csv`, `roc.png`,
`feature_importance.png`, and the trained model JSON.  If `xgboost`,
`scikit-learn`, or `tqdm` are missing, this subcommand exits with a dependency
message; the cut analysis does not require those optional packages.

### Linux Notes

The Python runner is written so it can move to Linux.  On Linux the script uses
`LD_LIBRARY_PATH` for optional OpenLoops/COLLIER runtime libraries.  The
macOS-specific `install_name_tool` patching is only used on macOS.

On `timur.kennesaw.edu`, Herwig is available through:

```bash
module load herwig/stable
```

Load the repository MG5 helper module before running the campaign:

```bash
module use "$PWD/modulefiles"
module load python/311
module load higgsssc/mg5
```

This sets `MG5_DIR`, exposes the safe `mg5_aMC` wrapper, and patches the MG5
launcher if needed.  Then the usual Timur campaign command is:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --mg5-dir "$MG5_DIR" \
  --nevents 10000
```

The runner loads `herwig/stable` automatically on Linux when `--herwig-env` is
not set.  If you prefer to show the Herwig module choice explicitly, use:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --mg5-dir "$MG5_DIR" \
  --herwig-module herwig/stable \
  --nevents 10000
```

If OpenLoops/COLLIER is needed for a loop-induced sample and is not found
automatically, pass it explicitly:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --collier-library /path/to/libcollier.so \
  --nevents 10000
```

### Quick Troubleshooting

If MG5 cannot find the model, check that the model directory exists:

```text
MG5_aMC_v3_5_15/models/loop_sm_haa
```

If Herwig cannot find the PDF set, either install the PDF set or choose one
that is installed:

```bash
python3 hgammagamma/run_gammagamma_campaign.py \
  --herwig-pdf NNPDF31_nnlo_as_0118
```

If the report script says it found no samples, check that the campaign produced
matching `.top` and `.dat` files for the requested run tag.

If a reducible background has zero diphoton events, remember that fake-photon
construction is not implemented yet.  That is expected for `j j` and
`e+ e-` until the analysis code learns how to reinterpret jets or leptons as
photons.
