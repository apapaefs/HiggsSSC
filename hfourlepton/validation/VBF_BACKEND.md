# VBF backend and failure recovery

The higher-order four-lepton workflow uses POWHEG `VBF_H` for VBF and separate
MG5 samples for WH, qqZH, ggZH and ttH. The previous MG5 VBF diagram exclusion
is no longer generated. Use a new campaign tag when moving to this backend.

## Evidence from Timur, 26 October 2026

The failed `fourl_lhc13_pilot` VBF export used MG5 3.5.15, `loop_sm-no_b_mass`,
unitary gauge and
`p p > h j j $$ w+ w- z QCD^2=0 QED^2=6 [QCD]` with mixed expansion enabled.
Several high-rank virtual integrals could not be reduced by Ninja/CutTools.
Adding IREGI allowed further progress, but the run eventually reported
`double free or corruption (out)` in `P0_us_hus/GF1.0/log.txt`, at
21:43:38 UTC. Its MG5 parent stayed alive after `Stop all computation`.
The stripped executable did not establish which library routine corrupted
memory. A slow parent after this explicit failure is not a working integration.

An authorized diagnostic used a separate copy of MG5, Feynman gauge and the
standard `6|7|1` reducer list. The process definitions kept the squared orders,
13 TeV proton beams, mH=125 GeV, NNPDF4.0 NLO QED (335900) and unchanged
pole-test tolerances. Four representative channels were tested:

| Configuration | uu→Huu | us→Hus | uū→Hdd̄ | bb→Hbb |
| --- | --- | --- | --- | --- |
| Feynman gauge, original W/Z widths, VH exclusion | 0/20 | 20/20 | 0/20 | 0/20 |
| Feynman gauge, W/Z widths zero, VH exclusion | 20/20 | 20/20 | 0/20 | 20/20 |
| Feynman gauge, W/Z widths zero, full EW Hjj | — | — | 20/20 | — |

Entries count passed pole-check points. The last trial was compile/test only;
it did not integrate or generate events. The sequence recorded 392.54 seconds
within a 900-second diagnostic limit. The original campaign manifest, managed
state and MG5 configuration had identical checksums before and after it.
Detailed cards and logs are on Timur under
`hfourlepton/diagnostics/vbf-feynman-20261026/`. The `passed` field in its
`results.json` refers only to the last, inclusive subprocess pole test.

This establishes that changing reducer or gauge alone is insufficient for
the old sample definition. Removing the VH exclusion would change the physics
sample and overlap with the separate VH accounting. The selected replacement
is the dedicated VBF approximation described by
[Nason and Oleari](https://arxiv.org/abs/0911.5299). It omits s-channel VH and
identical-fermion interference. The omitted contributions are not evaluated
by this workflow. No claim of numerical CMS/ATLAS agreement follows from these
diagnostics.

## Reproducible installation and input checks

The managed installation pins:

| Component | Revision |
| --- | --- |
| POWHEG BOX V2 | `d20df70ba6708da3978d5ab32d91b2a069457c1a` |
| VBF_H | `b9afa1e52f8e90d20d645ea48226588e35642793` |
| MG5 | `1d7559658bdc2b7eba476d4897aa594994fa82af` (3.5.15) |

`setup --workflow fourlepton-ho` initializes the VBF_H submodule in the managed
POWHEG checkout. Its old makefile uses `F77`; setup supplies the configured
Fortran compiler and `-fallow-argument-mismatch` for the legacy interfaces.
The supported default analysis linkage supplies the histogram routines and
FastJet wrappers. A generated `svn.version` include embeds the exact Git
revisions in the LHE header. The calculation sources are unchanged. Setup
records the build command and executable hash; check/run reject a changed
installation.

An isolated build succeeded on Timur with `herwig/stable`, GCC/gfortran 11.5,
Herwig 7.3.0 and the above revisions. Its executable SHA-256 was
`ba8a7efb288d3ac0cbc5e08e9c2066f48bda83667004369b2fbdc2a183ed066c`.
Build evidence is under `hfourlepton/diagnostics/vbf-h-build-20261026/`.
This build was not used to generate a VBF event sample. Binary hashes can
differ across build locations and toolchains; each installation records its
own hash.

The VBF runner uses the diphoton POWHEG integration/seed driver with separate
cards, seeds and grids. It checks the source pins, executable, input card,
event count, beam energy, PDF, native rate, SCALUP, stable-H state, finite-width
mass bounds and seven finite scale weights. Resume requires the same
configuration; completed output is checked by SHA-256. VBF rates come from
the native NLO LHE integral, multiplied once by the Higgs decay branching
fraction. ihixs applies only to ggF. For explicit settings and the approximation
used for the decay branching fractions, see [the campaign guide](../HIGHER_ORDER.md).

## Student recovery steps

1. Run `python3 run_pipeline.py status --workflow fourlepton-ho --energy-tev 13
   --run-tag TAG` and read the named log. A fatal message takes precedence over
   a stale `running` state. Long silence by itself does not prove failure.
2. After an explicit failure, stop the old command with Ctrl-C in its terminal
   or screen session. Do not start setup while that command owns the installation
   lock. Keep the failed export and logs for diagnosis.
3. Update the checkout and run `python3 run_pipeline.py setup --workflow
   fourlepton-ho`, then `check` with the same workflow. Setup prepares software;
   it does not run production.
4. Prepare and run a fresh 30-event pilot tag. The
   [recovery command](../HIGHER_ORDER.md#prepare-and-launch) reuses the previous
   ggF LHE and ihixs record while producing the new VBF sample. Do not reuse an
   MG5 VBF LHE as POWHEG input or change its declared matching scheme.
5. Inspect the pilot report and weight/event checks before running larger
   samples. A 30-event pilot does not provide enough statistics for the numerical
   LHC validation gate.

For MC@NLO, `logs/compile.log` contains the numerical checks and
`SubProcesses/P*/check_poles.log` contains each subprocess result. Missing,
zero-point, partial or aborted checks stop the workflow before integration.
Explicit worker failure messages stop that stage's owned process group and
write `nlo_runtime.json`; unrelated processes are not signalled. Do not disable
the checks, loosen tolerances or edit completion markers to bypass a failure.

Synthetic tests exercise wrong-process/PDF/energy rejection, finite-width
boundaries, signed rates, disjoint decay inputs, optional weights, resume
protection and process cleanup. Full VBF generation, MadSpin/Herwig closure and
the 13/40 TeV campaign results must still be checked when the user runs them.
