# Higher-Order Gamma-Gamma Campaign

The campaign runner is [`../run_gammagamma_ho_campaign.py`](../run_gammagamma_ho_campaign.py).
It showers the existing POWHEG HJMiNNLO signal with Herwig, forces
`H -> gamma gamma`, and generates the three backgrounds used in the LO
detector-response study with MadGraph5_aMC@NLO and Herwig. All four samples
use the same HwSim reconstruction and
[SSC/GEM response analysis](../SSC_DETECTOR_RESPONSE.md).
The checked local configuration, test results and remaining validation work
are recorded in [`VALIDATION.md`](VALIDATION.md).
The synchronized Timur checkout, input paths, and 100,000-event command
are documented in [`TIMUR.md`](TIMUR.md).
On Timur, use `/home/apapaefs/Projects/HiggsSSC` for the LO and HO analyses,
ihixs, inputs and reports. The HO LHE/PDF inputs are in
`hgammagamma/HOAnalysis/inputs/`, the normalization calculation is in
`hgammagamma/HOAnalysis/normalization/`, and campaigns and reports are in
`HOAnalysis/runs/` and `HOAnalysis/plots/`. These generated files and the
local MG5 installation are ignored by Git; source and settings share the
same GitHub repository. The ihixs source remains the pinned
`external/ihixs` submodule within this checkout.

`HJ/HJMiNNLO` is built from a Higgs-plus-jet process, but the MiNNLOPS
construction makes it appropriate as an inclusive `gg -> H` NNLO+PS signal
sample. Do not treat it as requiring an analysis jet.

The HO diphoton signal now requires an inclusive N3LO QCD normalization from
[ihixs](https://github.com/dulatf/ihixs), linked as the pinned
`external/ihixs` submodule. The event shapes remain HJMiNNLO NNLO+PS. This
normalization applies only to `signal_gg_h_aa` in an HO campaign; backgrounds,
LO analyses and the four-lepton analysis retain their existing rates.
The new ihixs workflow and regression tests have not yet been run. Earlier
validation results in `VALIDATION.md` and `TIMUR.md` predate this change.

## Processes And Conventions

| Sample | Hard process and matching | Hard/shower PDF | Response |
| --- | --- | --- | --- |
| `signal_gg_h_aa` | POWHEG HJMiNNLO, followed by Herwig `h -> gamma gamma` | NNPDF40_nnlo_as_01180_qed, 336100 | genuine |
| `bkg_prompt_aa` | `p p > a a [QCD]`, MC@NLO | NNPDF40_nlo_as_01180_qed, 335900 | genuine |
| `bkg_gamma_j` | `p p > a j [QCD]`, MC@NLO | same NLO PDF | gammajet |
| `bkg_dy_ee` | `p p > e+ e- [QCD]`, MC@NLO | same NLO PDF | dielectron |

Defaults are proton beams of 20 TeV each, `mH = 125 GeV`, the Herwig
angular-ordered shower, hadronization and the installed underlying-event
tune. The MPI PDF is retained from that tune. The backgrounds use
`loop_sm-no_b_mass` and five-flavour `p`/`j` definitions, including `b` and
`b~`, consistently with the NLO PDF. This explicitly differs from the old
LO cards' four-light-flavour parton lists. Drell--Yan includes both the
virtual photon and Z contributions and their interference.

MG5 uses `parton_shower = HERWIGPP`, `ickkw = 0`, and `event_norm = average`.
`--parton` postpones the shower but retains MC@NLO subtraction terms; this
is not a fixed-order event sample. The background template applies the
recoil and spin settings from Herwig 7.3's distributed `LHE-MCatNLO.in`.
The POWHEG template follows `LHE-POWHEG.in` and uses SCALUP as a shower veto.
Both preserve signed variable weights and forbid recycling the LHE file.
This is the external-LHE angular-ordered shower interface; no additional
truncated-shower implementation is supplied here.

Background central scales are explicitly `muR = muF = HT/2`
(`dynamical_scale_choice = 3`), with MG5 scale-reweight information retained
in the LHE. The current HwSim/response output analyzes the nominal weight;
it does not yet propagate the full scale/PDF variation ensemble.

NLO prompt-photon samples need an infrared-safe definition. The starting
generation settings are smooth-cone isolation (`gamma_is_j = False`,
`ptgmin = 10 GeV`, `R0gamma = 0.4`, `epsgamma = 1`, `xn = 1`, `isoEM = True`),
`|eta| < 6`, anti-kt jets with R = 0.4, and a 10 GeV Born-jet threshold
only for gamma+jet. There is no generation jet threshold for gamma-gamma
or Drell--Yan. Drell--Yan has `pT(l) > 10 GeV` and `m(ll) > 30 GeV`.
Generation cuts are configurable with `--gen-*` and `--isolation-*`.
Check cut and isolation dependence before using the results for physics:
the smooth-cone definition is distinct from the GEM jet-based isolation
proxy, and photon fragmentation contributions are not included.

The loop-induced `gg -> gamma gamma` box sample was disabled in the LO
campaign and remains separate from the NLO quark-initiated continuum.
Dijet double fakes are not enabled. Gamma+jet retains the existing
single-jet fake choice and quark/gluon matching model, even with an extra
NLO parton; its flavour and multiple-jet ambiguities remain detector-model
systematics. See the response reference for its limitations.

## Calculate The HO Signal Normalization

Run the commands below yourself from the repository root, before an HO
signal `analyze` or `all` stage. No numerical cross section is supplied as
a placeholder. The analysis refuses to substitute the native LHE rate
when the required ihixs record is absent or incompatible.

Initialize the pinned source:

```bash
git submodule update --init external/ihixs
```

The wrapper needs CMake, a C/C++ compiler, LHAPDF 6, Boost headers and
Cuba 4.2. Install Cuba in a prefix containing `include/cuba.h` and
`lib/libcuba.a`, with compilers compatible with the LHAPDF installation.
`--cuba-dir` selects that prefix; `--boost-dir` optionally selects the
directory containing the `boost/` headers. Cuba was not found in the
inspected Timur Herwig prefix, so provide your own installation. The
wrapper makes separate upstream and production source/build copies;
the pinned submodule is kept unchanged.

On Timur, install [Cuba 4.2.2](https://feynarts.de/cuba/) in your own
prefix before building ihixs. Run this in a separate shell, or return to
the HiggsSSC repository afterwards:

```bash
source /etc/profile.d/modules.sh
module load herwig/stable
CUBA_PREFIX="$HOME/.local/cuba-4.2.2-gcc11"
mkdir -p "$HOME/.local/src"
cd "$HOME/.local/src"
curl -fL https://feynarts.de/cuba/Cuba-4.2.2.tar.gz -o Cuba-4.2.2.tar.gz
tar -xzf Cuba-4.2.2.tar.gz
cd Cuba-4.2.2

CC=/usr/bin/gcc CFLAGS="-O2 -fPIC" ./configure \
  --prefix="$CUBA_PREFIX" \
  --libdir="$CUBA_PREFIX/lib" \
  --includedir="$CUBA_PREFIX/include" \
  --with-real=8
make -j1 lib
make -j1 install TOOLS_DEFAULT=
ls -l "$CUBA_PREFIX/include/cuba.h" "$CUBA_PREFIX/lib/libcuba.a"
```

The standard double-precision interface (`--with-real=8`) matches ihixs.
Build Cuba serially: its archive-member recipes can update the same
`libcuba.a` concurrently. `TOOLS_DEFAULT=` skips the optional Qt viewer.
These commands install a library dependency; they do not calculate a
Higgs cross section.

On Timur, activate the runtime and install the additional PDF sets:

```bash
source /etc/profile.d/modules.sh
module load herwig/stable
cd /home/apapaefs/Projects/HiggsSSC
mkdir -p hgammagamma/HOAnalysis/inputs/lhapdf
export LHAPDF_DATA_PATH="$PWD/hgammagamma/HOAnalysis/inputs/lhapdf:/home/shared/Herwig/share/LHAPDF"
lhapdf update
lhapdf install NNPDF40_an3lo_as_01180_qed_mhou
lhapdf install NNPDF40_nnlo_as_01180_qed
lhapdf install PDF4LHC15_nnlo_100
```

The first set is the approximate N3LO QCD, NLO QED, five-flavour NNPDF4.0
fit with `alpha_s(MZ) = 0.118`. It is the highest-order matched NNPDF4.0 QED
choice used here. Member 0 is central; members 1--100 supply the PDF
uncertainty. The NNLO QED set supplies comparison rates and remains the
event-generation/shower PDF. PDF4LHC15 is used only for the upstream
13 TeV benchmark. See the
[NNPDF QED documentation](https://nnpdf.mi.infn.it/nnpdf4-0-qed/) and
[aN3LO set metadata](https://lhapdfsets.web.cern.ch/current/NNPDF40_an3lo_as_01180_qed_mhou/NNPDF40_an3lo_as_01180_qed_mhou.info).

Choose your actual Cuba prefix below. On the laptop replace the module
option with `--herwig-env /path/to/Herwig/prefix`, and select compilers
compatible with that stack. `--powheg-input` can point to the actual
production `powheg.input`; its physics settings must match the tracked
SSC reference card. Event counts, random seeds and integration statistics
do not change the inclusive-rate profile.

```bash
CUBA_PREFIX="$HOME/.local/cuba-4.2.2-gcc11"
IHIXS_OPTIONS=(
  --herwig-module herwig/stable
  --lhapdf-dir /home/shared/Herwig
  --cuba-dir "$CUBA_PREFIX"
  --cc /usr/bin/gcc --cxx /usr/bin/g++ --jobs 8
  --powheg-input hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input
)

# Optional: inspect the 109-point plan without running or writing anything.
python3 hgammagamma/run_ihixs_normalization.py --dry-run "${IHIXS_OPTIONS[@]}"

python3 hgammagamma/run_ihixs_normalization.py --stage build "${IHIXS_OPTIONS[@]}"
python3 hgammagamma/run_ihixs_normalization.py --stage benchmark "${IHIXS_OPTIONS[@]}"
python3 hgammagamma/run_ihixs_normalization.py --stage calculate "${IHIXS_OPTIONS[@]}"
```

If configuration reports `Cuba not found`, inspect
`HOAnalysis/normalization/ihixs-ssc40/build-upstream/configure.log` and
check that `--cuba-dir` points to the installed prefix above. A literal
`/path/to/your/cuba-prefix` is a placeholder and will not work. After
installing Cuba, repeat `--stage build` in the same work directory;
the failed configure has not produced a completed build manifest.
The message `Prepared 109 integrations` means input cards were written,
not that the integrations ran. LHAPDF's Python `tarfile` extraction
warnings are separate from the CMake dependency failure.

If a saved production point fails the 0.05% numerical-error check, resume
with precision refinement enabled:

```bash
python3 hgammagamma/run_ihixs_normalization.py \
  --stage calculate "${IHIXS_OPTIONS[@]}" --resume --refine-failed
```

Keep the original settings and work directory. Passing points and the
benchmark are reused. Only a point that exceeds the total-rate precision
limit is retried: its first retry uses `epsrel = min(original/10, 1e-6)`,
then each further retry tightens it by another factor of ten, with at
most three retries. `epsabs` is reduced by the same ratio. Each attempt
has its own `runs/LABEL__precision_N/` directory; the original result and
all retry outputs are preserved. Changed input, executable, PDF or
coupling provenance is refused. A Cuba convergence failure stops the
calculation and leaves its log for inspection.

This targets the actual integration controls in the pinned ihixs source:
the pure EFT terms use Cuhre, which fixes its minimum evaluations to
1,000 and does not use the Vegas `nstart`/`nincrease` controls. Its
per-term relative targets also have relaxation factors, reaching 10,000
for the N3LO `q1q2` term. Increasing the Vegas statistics alone therefore
need not improve the final error. The 0.05% check applies to the complete
raw EFT cross section and remains mandatory after refinement. The final
record stores each selected run's effective integration settings and
input/output hashes, including the selected retry directory.

The benchmark reproduces the unmodified ihixs example's raw
`eftn3lo = 45.1816 pb` at 13 TeV, within 0.5%; this is **not** the 40 TeV
rate. Production uses `sqrt(s) = 40000 GeV`, `mH = 125 GeV`, an on-shell
top mass of 173.2 GeV, the native POWHEG Fermi constant and unit conversion,
and central `muR = muF = mH/2 = 62.5 GeV`. These fixed inclusive scales
replace the native event-dependent MiNNLO scales for the total-rate
calculation. The finite-mass Born rescaling, exact mass corrections,
electroweak corrections and resummation are disabled to match the signal's
pure HEFT hard model. In particular the wrapper reads the raw `eftn3lo`
entry, not ihixs's Born-rescaled `Higgs XS`.

The production source copy takes the hard coupling directly from
`LHAPDF::PDF::alphasQ(muR)`, checked against the linked LHAPDF library for
every run. This deliberately upgrades the native NNLO-PDF/HOPPET
three-loop coupling to the aN3LO PDF's four-loop coupling. The wrapper
also checks Cuba convergence and corrects ihixs's rounded prefactor to
the native POWHEG convention. See
[`ihixs-ssc40.json`](ihixs-ssc40.json) for the full settings and the
[ihixs reference](https://arxiv.org/abs/1802.00827).

The calculation runs seven scale points (central plus the standard
six-point variations), all 100 PDF replicas, and NNLO/N3LO comparisons
with the native NNLO PDF: 109 production integrations in total. Each must
have a relative numerical error of at most 0.05%. Only after the benchmark,
all integrations and provenance checks pass is the record written to:

```text
hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json
```

Builds, inputs, logs, raw results and intermediate manifests are preserved
under `HOAnalysis/normalization/ihixs-ssc40/`. Add `--resume` to repeat a
benchmark or calculation after interruption; verified completed points
are reused and interrupted directories are archived. Changed physics,
PDF files or base integration settings require a fresh `--work-dir` and
`--record` path. The opt-in precision retries above preserve the base
plan and use separate, individually hashed attempt directories.
Preserve the record and raw calculation directory
together when transferring them between hosts. The signal sidecar also
embeds the complete record so reports remain portable.

## Prepare And Run

Run these commands from the repository root. The default `prepare` stage
writes cards and a manifest without starting generators:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --run-tag ho_run_01 --nevents 10000 --nb-core 4 \
  --herwig-module herwig/stable \
  --signal-lhe "$PWD/hgammagamma/HOAnalysis/runs/ho_run_01/powheg/powheg-hjminnlo-merged.lhe"
```

For the signal, first generate LHE files with the existing wrapper described
below, or supply an existing complete merged file:

```bash
python3 hgammagamma/run_powheg_hjminnlo.py \
  --nevents 10000 --jobs 4 --herwig-module herwig/stable \
  --run-dir "$PWD/hgammagamma/HOAnalysis/runs/ho_run_01/powheg"
```

Then run the signal and all three backgrounds:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage all --run-tag ho_run_01 --nevents 10000 --nb-core 4 \
  --herwig-module herwig/stable \
  --signal-normalization hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json \
  --signal-lhe "$PWD/hgammagamma/HOAnalysis/runs/ho_run_01/powheg/powheg-hjminnlo-merged.lhe"
```

Use the same runtime options at preparation and execution. On macOS use
`--herwig-module herwig/730`, or replace the module option with
`--herwig-env /path/to/Herwig/prefix` (an activation script also works).
`--no-herwig-module` uses the current environment. Select a compatible
MG5 Python with `--mg5-python`; Python 3.11 needs the `six` package. If an
activation script points to a removed compiler, explicitly set
`--mg5-fortran /path/to/gfortran --mg5-cxx /path/to/g++`. These overrides
also replace FC/F77/CXX after activation. `--analysis-cxx` selects a compiler
compatible with ROOT; the macOS default is Apple clang.
MG5 exports dependencies internally by default so that CutTools/IREGI can
be rebuilt in the process directory after a compiler upgrade. A consistent
site installation can use `--mg5-dependencies external` to share those libraries.
The runner adds the missing `external virtgranny_red` declaration to the
exported MG5 3.5.15 `genps_fks.f`, needed by gfortran 16. This declaration
does not change its phase-space algorithm or modify the shared MG5 source.
It also resolves `fastjet-config` and `lhapdf-config` inside the selected
runtime and records their absolute paths in `cards/runtime-*.mg5`.
For resumed exports it also updates these paths in the process-local
`Cards/amcatnlo_configuration.txt`, which MG5 reloads at launch.
MadLoop needs loop-reduction libraries for the NLO virtual corrections.
MG5 can try downloading optional Ninja and Collier libraries at the first
loop-process export, even when Ninja is already installed elsewhere.
Configure the local MG5 installation's `input/mg5_configuration.txt` with
the installed Ninja **library directory**, containing `libninja.a` and
the OneLOop library, and the corresponding `../include/mninja.mod`.
The installed library paths and recovery command for Timur are documented
in [`TIMUR.md`](TIMUR.md#recover-from-the-ninja-installer-error).
See [MG5's loop-library selection](https://github.com/mg5amcnlo/mg5amcnlo/blob/3.x/madgraph/interface/loop_interface.py).

For a first test, use a separate tag with `--nevents 30`. The signal may
consume a prefix of a larger merged file; its reference cross section is
calculated from the complete LHE sample. NLO integration still takes time
even when requesting few events. `--run-samples backgrounds` selects only
the three backgrounds; individual names or comma-separated lists also work.

Other stages are `build` (MG5 process export only), `generate` (background
LHE generation and signal LHE validation), `shower`, `analyze`, and
`normalize` (existing HO signal rate only; see below).
`--dry-run` shows the selected stages without writing files. To continue a
campaign, repeat the identical configuration with `--resume`; completed
stages are reused. Changed physics or runtime configuration requires a new
tag or output directory. Existing partial ROOT shower products are preserved
and are never silently overwritten. A finalized shower that consumed its
entire finite LHE input can be recovered with `--stage shower --resume`
after auditing its end-of-file exception, Herwig statistics, source weights
and ROOT files. Other partial or interrupted outputs remain blocked for
inspection.

Outputs live under `HOAnalysis/runs/RUN_TAG/{Signal,Backgrounds}/events/SAMPLE/`.
Each sample has its generator cards, logs, `campaign.json`, Herwig ROOT
files, analysis `.dat`/`.top`/`_var.root` outputs, and
`normalization-RUN_TAG.json`. `--output-dir` can relocate the campaign.
The shared analysis executable is built from `LOAnalysis/Code`; the LO
campaign directories and their results are not used as HO outputs.
New showers also record `shower_completion` in the campaign manifest.
This keeps the requested success target, attempted LHE records, saved
ROOT entries and discarded events separate. `HwSim:OnTheFlyAnalysis No`
saves every delivered event; its historical "pass basic cuts" footer is
the saved-event count in this mode. The analysis must read exactly the
recorded saved count, which can be smaller than `--nevents` after shower
failures. The original event request and generation fingerprint are retained.

## Normalization And Reports

The signal LHE must contain one undecayed status-1 Higgs per event, the
requested beam energy, and the NNLO PDF. Herwig selects only the diphoton
decay with a unit proposal branching fraction. The physical
`BR(H -> gamma gamma) = 0.00227` is applied once by the response analysis
(`--higgs-br` overrides it). No LO signal K-factor is used. Backgrounds
have unit extra weight scale.

The LHE checker records positive/negative weight sums, sum of squared
weights, effective event count, beam/PDF and matching provenance, and rejects
incomplete files or mismatched inputs. HJMiNNLO can leave XSECUP and PDF IDs
as `-1` in `<init>` with `IDWTUP = -4`. In that case the signed mean XWGTUP
sets the native production cross section, with its sampling error, and the PDF IDs
come from the embedded POWHEG card. It never uses the absolute-weight sum as
the physical cross section. For MG5, the integrated `<init>` cross section
is used. The post-analysis additionally records negative-weight diagnostics
and checks response-hypothesis weight closure.

For the HO signal, the native LHE rate is retained as a diagnostic. The
sidecar's production cross section comes from the validated ihixs record.
The runner checks the 40 TeV, 125 GeV, native NNLO QED PDF and pure HEFT
settings before accepting that rate. Negative event weights are retained;
the absolute-weight sum never supplies a physical efficiency or rate.
The signal normalization record has a separate fingerprint, so changing
the inclusive rate does not invalidate completed generation and showers.

To update an already analyzed HO run without rebuilding, regenerating,
showering or rerunning the detector analysis, use:

```bash
RUN_TAG=ho_run_01  # replace with your completed HO run tag at its original location
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage normalize \
  --output-dir "hgammagamma/HOAnalysis/runs/$RUN_TAG" \
  --run-tag "$RUN_TAG" --run-samples signal_gg_h_aa \
  --signal-normalization hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json
```

This reads the event count, BR and LHE path from the existing manifest,
checks the saved analysis and inputs, then updates only the signal
normalization sidecar and campaign rate provenance. It needs neither
`--resume` nor a Herwig runtime option. Old native-rate HO sidecars must
be migrated this way before reporting or resuming their analysis.

Generate the usual report using the HO campaign as its analysis root:

```bash
python3 hgammagamma/make_gammagamma_report.py \
  --analysis-root hgammagamma/HOAnalysis/runs/ho_run_01 \
  --run-tag ho_run_01 \
  --output-dir hgammagamma/HOAnalysis/plots/ho_run_01 \
  --no-density --normalization event_xsec --luminosity-fb 100
```

The report reads the HO normalization sidecars rather than requiring a
MadGraph LO banner. The selected rate is
`sigma_production * weight_scale * sum_selected_signed_weight / normalization_sum_weight`.
For a new shower, `normalization_sum_weight` is the signed weight sum of
all attempted source LHE records in the same Herwig/BR weight convention
as the analysis. The reader scans the full file to determine its weight
maximum; only consumed records enter the denominator. Discarded shower
events therefore have zero simulated response rather than being silently
removed from the normalization population. No unweighted count correction
is applied. Existing campaigns without the new provenance retain their
saved-analysis denominator and exact requested-count checks.
For the HO signal, `sigma_production` is the inclusive ihixs N3LO rate and
`weight_scale` is the physical diphoton BR, included exactly once. With
`--luminosity-fb L`, the HTML and summary CSV give physical selected yields
as `1000 * L * sigma_selected_pb`. Replace 100 in the example with your
chosen luminosity in inverse femtobarns. Plot signal magnification and
unit-area display options do not change these yields. The summary records
the native LHE rate, ihixs fingerprint and separate scale/PDF uncertainties;
the cross-section error column contains numerical integration error only.
New manifests also supply the source denominator and shower event counts,
exception statistics and momentum-consistency diagnostics. These losses
are distinct from detector selection. Failed decays and consistency
warnings remain simulation-quality issues even when the finite input
has been completely processed.
Reports for these campaigns verify the recorded LHE, Herwig card, logs
and ROOT files at their original paths; run them from the canonical
Timur checkout with those artifacts present.
These uncertainties describe the inclusive rate; they do not include
acceptance, BR or detector uncertainties.
Small samples can have negative histogram bins; increasing statistics is
necessary before interpreting efficiencies or distributions. Existing LO
analysis cards with K-factors, or classifiers requiring positive training
weights, should not be reused blindly for these samples.

## Tests To Run Yourself

The new Python tests exercise raw-EFT parsing, completeness/provenance,
HO-only normalization, migration, signed weights, BR handling and yields.
They need no ihixs build, generators or ROOT:

```bash
python3 -m unittest discover -s tests -p 'test_ihixs_normalization.py' -v
python3 -m unittest discover -s tests -p 'test_ho_normalization_*.py' -v
python3 -m unittest discover -s tests -p 'test_ho_shower_completion.py' -v
```

Then run the existing gamma-gamma and HO regression suites with the
runtime/dependencies described in `VALIDATION.md`:

```bash
python3 -m unittest discover -s tests -p 'test_gammagamma*.py' -v
```

Start the new normalization and analysis workflow with a separate small
signal pilot before the full campaign. The upstream benchmark and all
109 integrations are required before a production normalization record
is accepted. No builds, tests, integrations or analyses were run while
implementing this change.

Configuration references are the installed Herwig 7.3 LHE examples and the
MG5 3.5.15 `Template/NLO/Cards/run_card.dat` and
`Template/NLO/MCatNLO/Scripts/MCatNLO_MadFKS_HERWIGPP.Script`.
See also the matching discussion in the
[MG5_aMC@NLO FAQ](https://amcatnlo.web.cern.ch/list_detailed2.htm)
(its old Python/platform instructions do not apply to MG5 3.5.15),
the [Herwig 7 release paper](https://arxiv.org/abs/1512.01178), and the
[smooth-cone isolation definition](https://arxiv.org/abs/hep-ph/9801442).

## What This Compiles

The POWHEG executable to build is:

```text
POWHEG-BOX-V2/HJ/HJMiNNLO/pwhg_main
```

Use this executable to generate Les Houches events. For Herwig 7, shower those
LHE files externally with Herwig. Do not use POWHEG's old
`main-HERWIG-lhef` target; that is for the legacy Herwig interface, not the
Herwig 7 workflow.

## Required Tools

The compile needs:

- `git`;
- `gfortran`;
- a C++ compiler, for example `c++` on macOS or `g++` on Linux;
- LHAPDF 6 with `lhapdf-config` in `PATH`;
- FastJet with `fastjet-config` in `PATH`;
- `zlib`.

Check the toolchain first:

```bash
command -v git
command -v gfortran
command -v lhapdf-config
command -v fastjet-config

lhapdf-config --version
fastjet-config --version
```

On Apple Silicon, make sure `gfortran`, LHAPDF, and FastJet are all built for
the same architecture:

```bash
file "$(lhapdf-config --libdir)/libLHAPDF.dylib"
file "$(fastjet-config --prefix)/lib/libfastjet.dylib"
```

If these show `x86_64` while the compiler is producing `arm64` objects, either
use a consistent x86_64 shell/toolchain or install arm64 builds of LHAPDF and
FastJet.

## Get POWHEG-BOX-V2

From the `HiggsSSC` repository root:

```bash
cd /path/to/HiggsSSC

git clone --filter=blob:none --sparse \
  https://gitlab.com/POWHEG-BOX/V2/POWHEG-BOX-V2.git \
  POWHEG-BOX-V2
```

Materialize only the pieces needed by `HJ/HJMiNNLO`:

```bash
git -C POWHEG-BOX-V2 sparse-checkout add \
  include \
  svnversion \
  MiNNLOStuff \
  HJ

git -C POWHEG-BOX-V2 submodule update --init HJ
```

Check that the expected files are present:

```bash
test -f POWHEG-BOX-V2/include/LesHouches.h
test -f POWHEG-BOX-V2/svnversion/svnversion.sh
test -d POWHEG-BOX-V2/MiNNLOStuff
test -d POWHEG-BOX-V2/HJ/HJMiNNLO
```

If `POWHEG-BOX-V2` already exists as a sparse checkout, run the
`sparse-checkout add` and `submodule update` commands above inside the existing
checkout instead of cloning again.

## Compile HJMiNNLO

Enter the process directory:

```bash
cd /path/to/HiggsSSC/POWHEG-BOX-V2/HJ/HJMiNNLO
mkdir -p obj-gfortran
```

On Linux, the default Makefile settings are usually sufficient:

```bash
make pwhg_main CXX=g++ STDCLIB=-lstdc++
```

On macOS, match the C++ standard library used to build LHAPDF and FastJet. For
the Herwig GCC stack in this workspace, use Homebrew GCC and `libstdc++`:

```bash
source /Users/apapaefs/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate
make pwhg_main \
  CXX=/opt/homebrew/bin/g++-15 \
  CC=/opt/homebrew/bin/gcc-15 \
  STDCLIB=-lstdc++
```

If LHAPDF and FastJet were instead built with Apple clang/`libc++`, use
`CXX=c++ STDCLIB=-lc++`.

For a clean rebuild:

```bash
make clean
rm -f pwhg_main obj-gfortran/*.o obj-gfortran/libfiles.a

# Linux
make pwhg_main CXX=g++ STDCLIB=-lstdc++

# macOS with the Herwig GCC stack in this workspace
source /Users/apapaefs/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate
make pwhg_main \
  CXX=/opt/homebrew/bin/g++-15 \
  CC=/opt/homebrew/bin/gcc-15 \
  STDCLIB=-lstdc++
```

Only run one of the final two `make` commands, depending on the platform.

## Common Compile Failures

Missing `LesHouches.h` means the top-level POWHEG `include` directory is not in
the sparse checkout:

```bash
git -C /path/to/HiggsSSC/POWHEG-BOX-V2 sparse-checkout add include
```

Missing `MiNNLOStuff` or HOPPET-related files means the MiNNLO auxiliary code
is not in the sparse checkout:

```bash
git -C /path/to/HiggsSSC/POWHEG-BOX-V2 sparse-checkout add MiNNLOStuff
```

Undefined symbols involving `LHAPDF::mkPDF`, `fastjet::sorted_by_pt`,
`std::__1`, or `std::__cxx11` usually mean inconsistent C++ link settings or
inconsistent library architectures. On macOS with the Herwig GCC stack in this
workspace, rebuild with:

```bash
make clean
rm -f pwhg_main obj-gfortran/*.o obj-gfortran/libfiles.a
source /Users/apapaefs/Projects/Herwig/Herwig-REAL-stable-gcc-full/bin/activate
make pwhg_main \
  CXX=/opt/homebrew/bin/g++-15 \
  CC=/opt/homebrew/bin/gcc-15 \
  STDCLIB=-lstdc++
```

Then check the linked libraries:

```bash
lhapdf-config --libs
fastjet-config --libs --plugins
```

## Minimal SSC Test Run

The pipeline wrapper for `HJMiNNLO` is:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 10000
```

By default this creates a run directory under
`POWHEG-BOX-V2/HJ/HJMiNNLO`, patches `numevts` from `--nevents`, enforces
`ebeam1 = ebeam2 = 20000d0` for `pp` collisions at 40 TeV, runs the POWHEG
stages, merges the per-seed `pwgevents-*.lhe` files into
`powheg-hjminnlo-merged.lhe`, and writes that merged path to
`powheg-lhe-files.txt`. The individual seed files are still listed in
`powheg-lhe-seed-files.txt`. Use `--jobs N` to split the requested total event
count over `N` POWHEG seed jobs. If `--nevents` is not divisible by `--jobs`,
the wrapper uses two stage-4 groups so the requested total is still exact. Use
`--ebeam` only for deliberate non-SSC studies, and `--no-merge-lhe` only if
you want downstream tools to consume the per-seed files directly.

To check progress from another terminal while the production command is
running, use the same run-defining options with `--status`:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --herwig-module herwig/730 \
  --status
```

For an auto-refreshing view, use:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --herwig-module herwig/730 \
  --watch-status 60
```

The status view reports completed POWHEG stages, latest log activity,
per-seed LHE files, event blocks written so far, and the merged LHE/manifest
once the run has finished.

If a resumed run reaches stage 4 and POWHEG reports that a
`pwgevents-000N.lhe` file already exists, rerun with `--resume`. The wrapper
will now skip complete seed outputs and archive incomplete/stale stage-4 LHE
files under `resume-backups/` before rerunning those seeds, so POWHEG gets a
clean output filename without deleting the partial file:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --herwig-module herwig/730 \
  --resume
```

For example:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --run-dir /path/to/HiggsSSC/POWHEG-BOX-V2/HJ/HJMiNNLO/run-ssc40-hjminnlo-100k
```

The wrapper can configure the Herwig/LHAPDF runtime either through an
activation script/prefix or through environment modules. On the laptop, use:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --herwig-module herwig/730
```

On `timur`, use:

```bash
python3 /path/to/HiggsSSC/hgammagamma/run_powheg_hjminnlo.py \
  --nevents 100000 \
  --jobs 8 \
  --herwig-module herwig/stable
```

Equivalently, set `HERWIG_MODULE=herwig/730` or
`HERWIG_MODULE=herwig/stable` before running the wrapper. If `HERWIG_ENV` is
set instead, it may point either to the Herwig stack prefix or to the activation
script itself; the wrapper normalizes a prefix such as
`/path/to/Herwig-REAL-stable-gcc-full` to
`/path/to/Herwig-REAL-stable-gcc-full/bin/activate`. By default, the wrapper
uses `herwig/730` on macOS and `herwig/stable` on Linux when no explicit
Herwig environment is supplied. After a failed setup-stage attempt, rerun the
same directory with `--resume`.

The recommended starting card for SSC 40 TeV `HJMiNNLO` production is:

```text
hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input
```

It uses `NNPDF40_nnlo_as_01180_qed` through LHAPDF ID `336100`, regenerates
grids by default, keeps negative weights, and sets the central HJMiNNLO options
for an inclusive `gg -> H` NNLO+PS signal. Copy it into a run directory as
`powheg.input`.

After `pwhg_main` is built, make a small 40 TeV test run directory:

```bash
cd /path/to/HiggsSSC/POWHEG-BOX-V2/HJ/HJMiNNLO

mkdir -p run-ssc-hjminnlo-test
cp suggested_run/pwgseeds.dat \
  run-ssc-hjminnlo-test/
cp /path/to/HiggsSSC/hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input \
  run-ssc-hjminnlo-test/powheg.input

cd run-ssc-hjminnlo-test

../pwhg_main
```

The full suggested run script uses multiple integration stages and many
parallel jobs. Use the small test first to catch path, PDF, and linker issues
before launching a large production run.

## Notes For The Gamma-Gamma Analysis

`HJMiNNLO` generates the Higgs production process. The HO campaign runner
handles `H -> gamma gamma` in Herwig and includes the physical branching
ratio once in the response analysis, as specified above.

For rate comparisons, do not apply the simple LO `ggH` K-factor used by the
current `LOAnalysis` signal sample. The HO signal's event shapes come from
HJMiNNLO and its inclusive production rate comes from the ihixs N3LO record.
