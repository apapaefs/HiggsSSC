# Student Guide: ihixs In HiggsSSC

This guide calculates the inclusive **pure-HEFT ggF production cross section
at 40 TeV, at N3LO in QCD**, for the HO diphoton analysis. The rate excludes
the Higgs branching fraction, generation cuts and detector response. It
normalizes only the HO `signal_gg_h_aa` sample. Event shapes remain those of
the HJMiNNLO NNLO+PS sample; the LO analyses and backgrounds keep their
existing rates.

Keep everything in one `HiggsSSC` checkout. The pinned ihixs source is its
`external/ihixs` submodule, and the wrapper applies the maintained fixes in
separate build-source copies. Run the commands below yourself on Timur, in
order. Integrations take time; use `screen` or `tmux` before starting them.

## 1. Clone Into Your Own Account

Use your own Timur username and home directory. No GitHub SSH key or access
to the paper submodule is needed for ihixs:

```bash
ssh YourUsername@timur.kennesaw.edu
mkdir -p "$HOME/Projects"
cd "$HOME/Projects"
git clone https://github.com/apapaefs/HiggsSSC.git
cd HiggsSSC
git submodule update --init external/ihixs
```

For an existing checkout, use its repository root instead of cloning again:

```bash
cd "$HOME/Projects/HiggsSSC"
git pull --ff-only origin main
git submodule sync -- external/ihixs
git submodule update --init external/ihixs
```

Do not run `git submodule add`, change the pinned ihixs commit, or initialize
`paper/` for this calculation. The required settings and POWHEG reference
card are already tracked; a signal LHE file is not needed to calculate the
inclusive normalization.

## 2. Load The Timur Runtime And Install Cuba

The wrapper needs Python 3.9 or newer, CMake, C/C++ compilers, Boost headers,
LHAPDF 6 and Cuba 4.2. Timur's `herwig/stable` stack provides LHAPDF in
`/home/shared/Herwig`. Use the matching GCC 11 system compilers shown below
for both Cuba and ihixs.

If you already installed Cuba with this compiler and the standard
double-precision interface, reuse that prefix. Otherwise install
[Cuba 4.2.2](https://feynarts.de/cuba/) in your own account:

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
cd "$HOME/Projects/HiggsSSC"
```

Build Cuba serially because its archive recipes can update `libcuba.a`
concurrently. `--with-real=8` selects the interface ihixs expects;
`TOOLS_DEFAULT=` skips the optional Qt viewer. Select an actual installed
prefix when running ihixs. The literal `/path/to/your/cuba-prefix` is not
an installation.

On another host, use a compatible LHAPDF/Cuba/compiler stack and replace
`--herwig-module herwig/stable` with `--herwig-env /your/Herwig/prefix` in
the later options. `--lhapdf-dir` and `--cuba-dir` select installation
prefixes; `--boost-dir` can select the directory containing `boost/` headers.

## 3. Install The PDF Sets In Your Checkout

Run from the repository root after loading the runtime:

```bash
cd "$HOME/Projects/HiggsSSC"
source /etc/profile.d/modules.sh
module load herwig/stable
PDF_INSTALL_DIR="$PWD/hgammagamma/HOAnalysis/inputs/lhapdf"
mkdir -p "$PDF_INSTALL_DIR"
export LHAPDF_DATA_PATH="$PDF_INSTALL_DIR:/home/shared/Herwig/share/LHAPDF"
lhapdf --listdir "$PDF_INSTALL_DIR" --pdfdir "$PDF_INSTALL_DIR" update
lhapdf --listdir "$PDF_INSTALL_DIR" --pdfdir "$PDF_INSTALL_DIR" install NNPDF40_an3lo_as_01180_qed_mhou
lhapdf --listdir "$PDF_INSTALL_DIR" --pdfdir "$PDF_INSTALL_DIR" install NNPDF40_nnlo_as_01180_qed
lhapdf --listdir "$PDF_INSTALL_DIR" --pdfdir "$PDF_INSTALL_DIR" install PDF4LHC15_nnlo_100
```

Keep the same `LHAPDF_DATA_PATH` in every calculation shell. `--listdir`
and `--pdfdir` put both the index and downloaded PDFs in your own checkout;
do not install PDFs into the shared Herwig prefix.
The aN3LO NNPDF4.0 QED set provides the production central value and 100
replicas, with `alpha_s(MZ) = 0.118`. The NNLO QED set supplies the two
comparison rates and remains the signal event-generation/shower PDF.
PDF4LHC15 supplies only the published 13 TeV benchmark. The wrapper checks
the installed PDF metadata, member files and hard coupling for each run.

## 4. Build, Benchmark, Then Calculate

The tracked profile is `hgammagamma/HOAnalysis/ihixs-ssc40.json`. It uses
`mH = 125 GeV`, an on-shell top mass of `173.2 GeV`, and central
`muR = muF = 62.5 GeV`. It matches the tracked POWHEG reference's pure
HEFT model, Fermi constant and unit conversion. Finite-mass corrections,
the finite-top Born rescaling, electroweak corrections and resummation are
disabled. The wrapper extracts raw `eftn3lo`, rather than the Born-rescaled
`Higgs XS`, and uses the PDF's `alphasQ(muR)` in the production build.

Use a distinct work directory for your calculation. The default final
record path below is what the HO analysis reads. If that record already
exists, preserve it and its work directory. To make a separate calculation,
choose a new `--record` path and pass that same path as
`--signal-normalization` to your HO analysis; a different existing record
will not be overwritten.

```bash
cd "$HOME/Projects/HiggsSSC"
CUBA_PREFIX="$HOME/.local/cuba-4.2.2-gcc11"
IHIXS_OPTIONS=(
  --settings hgammagamma/HOAnalysis/ihixs-ssc40.json
  --work-dir hgammagamma/HOAnalysis/normalization/ihixs-student01
  --record hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json
  --powheg-input hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf40nnloqed.input
  --herwig-module herwig/stable
  --lhapdf-dir /home/shared/Herwig
  --cuba-dir "$CUBA_PREFIX"
  --cc /usr/bin/gcc --cxx /usr/bin/g++ --jobs 8
)

# Inspect the plan without running or writing anything.
python3 hgammagamma/run_ihixs_normalization.py --dry-run "${IHIXS_OPTIONS[@]}"

for ihixs_stage in build benchmark calculate; do
  IHIXS_STAGE_OPTIONS=()
  if [ "$ihixs_stage" = calculate ]; then
    IHIXS_STAGE_OPTIONS=(--refine-failed)
  fi
  python3 hgammagamma/run_ihixs_normalization.py \
    --stage "$ihixs_stage" "${IHIXS_OPTIONS[@]}" \
    --resume "${IHIXS_STAGE_OPTIONS[@]}" || break
done
```

The loop stops at the first failed stage. `--jobs 8` controls compilation;
integrations run sequentially. `Prepared 109 integrations` means cards were
prepared, not that any cross section has been calculated.

Fresh builds apply both option-parser bounds fixes automatically to the
benchmark and production copies. They also retain the production coupling
adapter and strict Cuba convergence checks. You do not need
`--repair-parser` for a new student workspace. The pinned submodule itself
is kept unchanged. The benchmark retains the published example's physics
and must reproduce raw `eftn3lo = 45.1816 pb` at 13 TeV within 0.5%; this
number is not the SSC cross section.

The production calculation needs all **109 points**: central plus six
scale variations, 100 PDF replicas, and NNLO/N3LO comparisons using the
native NNLO PDF. Every selected result must meet the **0.05% relative
numerical-error limit** and the PDF/coupling checks. `--refine-failed`
reuses passing points and makes up to three tighter precision attempts
for a point that exceeds that limit. Each attempt has a separate directory
and provenance; the tolerance is never relaxed. Cuba failures or invalid
outputs stop the run for inspection.

Only this final message means the normalization was installed successfully:

```text
Installed validated inclusive production normalization: .../ggf-ssc40-n3lo.json
```

## 5. Inspect The Installed Record And Use It In HOAnalysis

After the calculation succeeds, inspect its saved summary without starting
another integration. If you selected a different record path, use it in
`Path(...)` below as well:

```bash
python3 - <<'PY'
import json
from pathlib import Path

path = Path("hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json")
record = json.loads(path.read_text())
print("Record:", path)
print("Cross section [pb]:", record["cross_section_pb"])
print("Numerical error [pb]:", record["cross_section_error_pb"])
print("Production points:", len(record["runs"]))
print("Benchmark passed:", record["validations"]["benchmark"]["passed"])
print("Coupling checked points:", record["validations"]["alpha_s"]["checked_runs"])
print("Theory uncertainties:", record["uncertainties"])
print("Artifacts relative to record:", record["provenance"]["artifact_root"])
PY
```

Expect 109 production points and 109 coupling checks. This command displays
the JSON; the wrapper and HO analysis perform the full validation. Preserve
the JSON **and its entire work directory** together, including build
manifests, input cards, logs, raw outputs and completion hashes. These
generated products, dependency installations and PDF files are ignored by
Git; a fresh GitHub clone does not include someone else's calculated rate.

Use the record in an HO campaign with:

```text
--signal-normalization hgammagamma/HOAnalysis/normalization/ggf-ssc40-n3lo.json
```

The response analysis applies `BR(H -> gamma gamma) = 0.00227` once, and
yields use your chosen luminosity and signed event weights. Follow the
[HO campaign README](README.md#prepare-and-run) to generate or supply
signal events, generate backgrounds, shower, analyze and make the report.
The [Timur guide](TIMUR.md#after-ihixs-finishes) records the instructor's
existing campaign; its `/home/apapaefs` paths are not your student paths.

The instructor's October 5, 2026 refined Timur run installed
`226.82244161097535 +/- 0.06711554812896317 pb` for this pure-HEFT profile,
with all 109 points, both original and repaired benchmarks, coupling checks
and the repaired central comparison passed. The quoted error is numerical;
scale and PDF uncertainties are separate. This documents that completed
run, rather than supplying a substitute record or certifying a fresh
student build. Builds, calculations and regression tests are left for you
to run in your own account.

## Resume And Troubleshoot

- After an interruption, reload the runtime and PDF path, restore the same
  `IHIXS_OPTIONS`, and repeat the loop above. Verified completed points are
  reused; interrupted attempts are archived. Keep the work directory and
  base settings unchanged.
- If CMake cannot find Cuba, inspect
  `hgammagamma/HOAnalysis/normalization/ihixs-student01/build-upstream/configure.log` and verify
  the prefix contains `include/cuba.h` and `lib/libcuba.a`. The configure
  failure can be retried after installing the missing dependency.
- A written `ihixs.out` after a nonzero exit is not a completed result.
  Inspect the named log. Do not manufacture a `complete.json` or bypass
  precision, coupling or provenance checks.
- A legacy workspace created before automatic parser fixes may need
  `--repair-parser --resume`. Preserve its original manifest and results;
  follow the [legacy recovery instructions](TIMUR.md#recover-ihixs-after-a-parser-heap-abort).
  Do not rebuild its completed binaries in place.
- Changed physics, PDF files, compiler/runtime or base integration settings
  require new work and record paths. The optional precision retries keep
  the original plan intact and record each numerical refinement.

Run the provided ihixs regression suites yourself from the repository root:

```bash
python3 -m unittest discover -s tests -p 'test_ihixs*.py' -v
```

These are regression tests for the wrapper and record handling; the build,
published benchmark and 109-point calculation provide the numerical checks.
