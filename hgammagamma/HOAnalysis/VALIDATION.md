# Local validation, 17 September 2026

The setup was checked with MadGraph5_aMC@NLO 3.5.15, Herwig 7.3.0,
ThePEG 2.3.0, LHAPDF 6.5.6, ROOT 6.40.02 and the installed HwSim plugin.
The local GCC installation is now version 16; the Herwig activation script
still exports paths to GCC 15. Explicit MG5 compiler overrides are therefore
needed on this machine. ROOT analysis compilation uses Apple clang.

Completed checks:

- 53 gamma-gamma unit/regression tests, including signed LHE normalization,
  matching/PDF/beam validation, undeclared or already-decayed Higgs rejection,
  MC@NLO card ordering, photon-isolation settings, resume protection, and
  preservation of negative/cancelling histogram bins.
- NLO process export for all three backgrounds, including their virtual and
  real matrix elements, using the five-flavour `loop_sm-no_b_mass` model.
- A 30-event HJMiNNLO → Herwig → HwSim → SSC/GEM signal run in
  `runs/ho_validation`, using the existing complete 100,000-event merged LHE.
  All 30 events reached the response analysis; three had negative weights.
  The response output-tree weight sum agrees with the input sum to floating
  point precision. The LHE-to-HwSim weight ratio is constant across all 30
  events, including the negative ones.
- A 30-event MG5 MC@NLO → Herwig → HwSim → SSC/GEM Drell--Yan run in
  `runs/ho_dy_smoke_v3`. All four subprocess directories passed their
  matrix-element and MC subtraction checks; virtual poles cancelled at
  20/20 test points per directory. All 30 generated events reached the
  response analysis, including four negative-weight events. The LHE-to-HwSim
  weight ratio is constant event by event, and response-hypothesis closure
  differs from the input weight sum by only `1.4e-14`. The dielectron
  response produced eight valid diphoton fake hypotheses.
- Repeating both completed campaigns with `--resume` reused all three
  stages after checking their source and output products.
- HTML/PNG reports were generated for the signal and Drell--Yan pilots,
  and the diphoton-mass plots inspected. Labels and normalization come from
  the higher-order campaign metadata.

These are software and interface checks, not a validation of production
statistics, matching/isolation dependence, or detector-model systematics.
The photon-background process exports have not yet undergone NLO integration
or showering. Full scale/PDF ensembles are not propagated into HwSim outputs.

ThePEG emits its generic `IDWTUP = -4` reader warning for both LHE files.
The `VarNegWeight` mode retains the event weights, as checked event by event
above; rates in the report use the independently checked LHE normalization.
Herwig also reported two EvtGen momentum-conservation warnings in the signal
pilot, with a maximum event-level violation of about 5 MeV; all 30 events
were written and analyzed.

## Prepared Local Campaign

`runs/ho_run_01` contains cards and manifests for 10,000 events per sample.
No production events have been generated for that tag. To run that exact
prepared configuration from the repository root:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage all --run-tag ho_run_01 --nevents 10000 \
  --signal-lhe "$PWD/POWHEG-BOX-V2/HJ/HJMiNNLO/run-ssc40-hjminnlo-nnpdf40nnloqed-100000ev/powheg-hjminnlo-merged.lhe" \
  --herwig-env /Users/apapaefs/Projects/Herwig/Herwig-REAL-stable-gcc-full \
  --mg5-fortran /opt/homebrew/bin/gfortran \
  --mg5-cxx /opt/homebrew/bin/g++-16
```

Append `--resume` to continue completed stages. Use a different tag if
changing event counts, cores, cuts, isolation, PDFs, or runtime options.
The general portable workflow and normalization conventions are in
[`README.md`](README.md).
