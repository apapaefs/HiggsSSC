# The ihixs normalization calculation

The [common HO workflow](../../README.md) installs Cuba/PDFs, builds ihixs,
benchmarks it and calculates the signal normalization. No manually entered
cross section or instructor record is required. A small event pilot still
needs the complete normalization calculation.

## What is calculated

The pinned [ihixs source](https://github.com/dulatf/ihixs) is initialized as
`external/ihixs`; only this public submodule is required. Settings are in
[`ihixs-ssc40.json`](ihixs-ssc40.json): 40 TeV pp, `mH=125 GeV`,
`muR=muF=mH/2`, pure HEFT N3LO QCD. Finite-mass and electroweak corrections
are not added. This profile is matched to the signal's POWHEG conventions.

The 109 production points cover the central result, seven-point scale scan,
100 PDF replicas and NNLO-PDF comparison. The inclusive rate uses
`NNPDF40_an3lo_as_01180_qed_mhou`; HJMiNNLO generation keeps
`NNPDF40_nnlo_as_01180_qed`. `PDF4LHC15_nnlo_100` is used for the separate
13 TeV upstream benchmark. See the [NNPDF QED description](https://nnpdf.mi.infn.it/nnpdf4-0-qed/)
and [aN3LO metadata](https://lhapdfsets.web.cern.ch/current/NNPDF40_an3lo_as_01180_qed_mhou/NNPDF40_an3lo_as_01180_qed_mhou.info).

The build uses separate benchmark and production source copies. Both receive
the option-parser bounds fixes; benchmark physics stays upstream. Production
also uses the LHAPDF hard-coupling adapter and strict Cuba error checks.
The validator requires the benchmarks, all production points, coupling checks,
precision checks and file provenance before accepting a record.

Use the raw `eftn3lo` result, not `Higgs XS` or `R_LO*eftn3lo`, which include
a finite-top Born rescaling. Production converts the ihixs prefactor to the
POWHEG GF convention using the factor specified in the settings. The inclusive
rate upgrades to the aN3LO PDF's four-loop coupling; the native HJMiNNLO
NNLO-PDF/HOPPET evolution remains its recorded three-loop convention.

## Inspect and preserve your result

For the common installation, inspect the record without simulation libraries:

```bash
python3 - <<'PY'
import json
from pathlib import Path
p = Path('inputs/pipeline/normalization/ggf-ssc40-n3lo.json')
r = json.loads(p.read_text())
print('Cross section [pb]:', r['cross_section_pb'])
print('Numerical error [pb]:', r['cross_section_error_pb'])
print('Production points:', len(r['runs']))
print('Benchmark:', r['validations']['benchmark']['passed'])
print('Coupling checks:', r['validations']['alpha_s']['checked_runs'])
print('Theory uncertainties:', r['uncertainties'])
PY
```

Expect 109 completed points and coupling checks. Displaying JSON is not a
substitute for the full validation performed by the pipeline. Preserve the
record **and** `inputs/pipeline/normalization/ihixs/`, including raw outputs,
logs, cards, build manifests and hashes. Generated calculations are ignored
by Git; `git pull` cannot supply a missing result.

## Advanced use

[`run_ihixs_normalization.py`](../run_ihixs_normalization.py) exposes
`prepare`, `build`, `benchmark`, `calculate` stages and `--dry-run`.
Use a new work directory and record for changed physics, compiler/runtime,
PDF files or base integration settings. `--resume --refine-failed` reuses
validated points and records numerical precision retries without changing the
original plan. Legacy parser-repair workspaces use `--repair-parser --resume`;
keep their original binaries and completion evidence.

Cuba must provide the double-precision API (`--with-real=8`); common setup
builds its archive serially and omits the optional Qt tools. See
[recovery](../../docs/TROUBLESHOOTING.md) for failure interpretation and
[validation](../../VALIDATION.md) for measured checks.
