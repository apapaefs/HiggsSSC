# Higher-Order Gamma-Gamma Signal Setup

This directory records the higher-order signal workflow for the future
`hgammagamma` analysis. The first target is inclusive gluon-fusion Higgs
production with `POWHEG-BOX-V2/HJ/HJMiNNLO`, showered later with Herwig 7 and
an `H -> gamma gamma` decay setup.

`HJ/HJMiNNLO` is built from a Higgs-plus-jet process, but the MiNNLOPS
construction makes it appropriate as an inclusive `gg -> H` NNLO+PS signal
sample. Do not treat it as requiring an analysis jet.

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

The recommended starting card for SSC 40 TeV `HJMiNNLO` production is:

```text
hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf31.input
```

It uses `NNPDF31_nnlo_as_0118` through LHAPDF ID `303600`, regenerates grids by
default, keeps negative weights, and sets the central HJMiNNLO options for an
inclusive `gg -> H` NNLO+PS signal. Copy it into a run directory as
`powheg.input`.

After `pwhg_main` is built, make a small 40 TeV test run directory:

```bash
cd /path/to/HiggsSSC/POWHEG-BOX-V2/HJ/HJMiNNLO

mkdir -p run-ssc-hjminnlo-test
cp suggested_run/pwgseeds.dat \
  run-ssc-hjminnlo-test/
cp /path/to/HiggsSSC/hgammagamma/HOAnalysis/powheg-hjminnlo-ssc40-nnpdf31.input \
  run-ssc-hjminnlo-test/powheg.input

cd run-ssc-hjminnlo-test

../pwhg_main
```

The full suggested run script uses multiple integration stages and many
parallel jobs. Use the small test first to catch path, PDF, and linker issues
before launching a large production run.

## Notes For The Gamma-Gamma Analysis

`HJMiNNLO` generates the Higgs production process. The `H -> gamma gamma`
decay should be handled in the later shower/decay step, or by a controlled
post-processing step, and the final normalization should include
`BR(H -> gamma gamma)` exactly once.

For rate comparisons, do not apply the simple LO `ggH` K-factor used by the
current `LOAnalysis` signal sample. `HJMiNNLO` is already the higher-order
production prediction.
