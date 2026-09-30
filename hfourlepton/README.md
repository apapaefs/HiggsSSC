# Higgs to four leptons at the 40 TeV SSC

This directory contains the LO campaign for

```text
pp -> H -> ZZ* -> 4e, 4mu, 2e2mu
```

and its neutral-current and heavy-flavour backgrounds.  The same generated
events can be analyzed with the `perfect` or `ssc` response profile.

## Quick start

Inspect the complete irreducible campaign without running external programs:

```bash
python3 hfourlepton/run_four_lepton_campaign.py \
  --sample-set irreducible \
  --detector-response perfect \
  --nevents 100 \
  --run-tag smoke \
  --dry-run
```

Run it after MG5, Herwig, HwSim, ROOT, and LHAPDF are available:

```bash
python3 hfourlepton/run_four_lepton_campaign.py \
  --sample-set irreducible \
  --detector-response perfect \
  --nevents 10000 \
  --run-tag lo_4l_v1
```

Repeat with `--detector-response ssc` to apply the GEM response in the common
post-analyzer.  Products live under `hfourlepton/runs/RUN_TAG/`, and each
response profile has a versioned JSON manifest.

`--sample-nevents NAME=N` may be repeated to override the global event count.
`--run-samples NAME1,NAME2` restricts the selected sample set.  Stable-H signal
production is shared by the three MadSpin decay channels rather than generated
three times.  The aliases `signal`, `bkg_qq`, `bkg_gg`, `bkg_zbb`, and
`bkg_ttz` select the corresponding groups.

Build the combined perfect-versus-SSC report with:

```bash
python3 hfourlepton/make_fourlepton_report.py \
  hfourlepton/runs/lo_4l_v1/manifest_perfect.json \
  hfourlepton/runs/lo_4l_v1/manifest_ssc.json \
  --output-dir hfourlepton/runs/lo_4l_v1/report \
  --luminosity-fb 100 \
  --signal-region 120 130 \
  --mass-range 70 200 \
  --mass-bin-width 5 \
  --strict
```

`--strict` still writes the report but exits with status 2 when an audit error
is present.  Semileptonically biased samples are displayed but excluded from
combined totals until closure is established; `--allow-unvalidated-bias`
provides an explicit diagnostic override.  Other useful options are
`--closure-tolerance EPS` and `--title TEXT`.

## External LHE and future NLO input

Use the external backend with one `KEY=PATH` assignment per sample or shared
production group:

```bash
python3 hfourlepton/run_four_lepton_campaign.py \
  --source-backend external-lhe \
  --external-lhe signal_gg_h_stable=/path/to/stable_higgs.lhe.gz \
  --external-cross-section signal_gg_h_stable=SIGMA_PB \
  --external-order NLO \
  --run-samples signal_gg_h_4e \
  --detector-response perfect \
  --run-tag external_h
```

The signal input must contain exactly one status-one stable Higgs with
`mH=125.09 GeV` in every consumed event; the event records, rather than a
possibly stale banner mass, are authoritative.  Banner-light POWHEG-style
inputs are augmented with the minimal MG5 process/run metadata and a
standalone-complete restricted-SM SLHA card needed by MadSpin. Existing SLHA
blocks are preserved only when all MadSpin-required SM entries are present;
the augmentation strategy and source/prepared checksums are recorded in the
manifest.

The runner then applies the same direct four-body MadSpin decay used for LO
production. Its
normalization is deliberately not a literal branching-ratio multiplier:
MadSpin already scales event weights by its internal branching ratio.  The
runner reads the stable and decayed LHE `<init>` cross sections (or an
explicit `--external-cross-section SAMPLE_OR_GROUP=PB` override) and passes

```text
weight_scale = [sigma(stable H) * target BR(channel)] / sigma(decayed LHE)
```

to the analyzer.  A real run fails closed if either cross section is absent
or a placeholder.  It also validates the proton beams and beam energies and
records the LHE `IDWTUP`, PDF identifiers, and declared perturbative order.
For signed or variable-weight LHE input, the normalization denominator is the
sum of the event weights that ThePEG's LHE handler will use,
`XWGTUP/XMAXUP[IDPRUP]`; for non-unit-weight schemes the full file is scanned
to reproduce the handler's maximum-weight normalization.  This is the
implemented seam for future stable-H NLO input.

The total reference branching fraction is
`BR(H -> 4e,4mu,2e2mu)=1.251e-4`, following the PROPHECY4F value quoted in
[ATLAS-CONF-2017-046](https://cds.cern.ch/record/2273853/files/ATLAS-CONF-2017-046.pdf).
Default channel fractions were integrated
with MG5_aMC 3.5.15 using direct SM `1 -> 4` matrix elements at
`mH=125.09 GeV`; they are stored with their provenance in every manifest and
can be overridden with `--signal-channel-fractions`.

## Reducible backgrounds and decay bias

`--sample-set reducible` selects `Zbb`, `ttZ`, and `ttbar` samples.  The
`ttZ` contribution is represented by eight disjoint strata: the two
`Z -> ee,mumu` modes crossed with the four charge-labelled
`W+ W- -> e/mu` flavour combinations.  Herwig's branching-ratio reweighter
logic provides their explicit normalization through the campaign's
LHE-compatible event-weight plugin.  Prompt `W/Z` factors are evaluated
after the perturbative cascade but before weak hadron decays; this avoids
mistaking virtual current resonances inside heavy-hadron decays for extra
prompt branching factors.  Nominal
samples always retain unbiased heavy-hadron decays.  The optional

```bash
--heavy-flavour-bias semileptonic
```

adds separate `_hfbiased` strata.  The runner discovers the direct
one-hadron electron and muon modes in Herwig's `HerwigBDecays.in`, multiplies
their decay-selector weights by four while leaving every other active mode
enabled, and installs
`HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter`.  For each
decayed parent this realizes the exact full-support proposal

```text
q_m = b_m p_m / sum_n(b_n p_n),
w_decay = sum_n(b_n p_n) / (b_m sum_n(p_n)),
```

where `b_m=4` for the selected direct semileptonic modes and `b_m=1`
otherwise.  Thus no Standard Model decay route is assigned zero proposal
probability.  The event weight is the product of the parent-level factors and
is applied after recursive hadron decays by the heavy-flavour-scoped
`HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter`, to both the nominal and
optional LHE weights.  The unbiased nominal
sample is retained as the closure control (its size can be set with
`--bias-control-events`).

The manifest always records `closure_validated: false` for accelerated strata.
There is intentionally no command-line switch that can claim closure.  These
full-support proposals remain excluded from nominal use until a weighted
comparison with their unbiased controls succeeds.  The report calculates the
biased/control ratios and statistical pulls as diagnostics, but never promotes
that diagnostic to a closure claim; accelerated strata are excluded from
physics totals by default.  Generic light-jet fake backgrounds are not modeled.

## Generator policy

- beams: 20 TeV per beam (`sqrt(s)=40 TeV`);
- hard-process accuracy: LO;
- PDF: `NNPDF40_lo_as_01180`, LHAPDF ID 331900;
- loose generation cuts: `pT(l)>5 GeV`, `|eta(l)|<3`,
  `m(l+l-)>4 GeV`;
- Higgs mass: explicitly patched to 125.09 GeV;
- HwSim saves stable objects, reconstructed objects, hard partons, and
  optional weights.

The SSC profile enables the documented pileup-noise term by default; it is
recorded as configured but inactive in the perfect profile.
`--no-pileup-noise` disables that term for a controlled variation while
retaining the thermal-noise contribution.

The generated continuum samples are full `Z/gamma*` matrix elements with
same-flavour interference.  Both the tree-level and loop-induced continuum
cards explicitly exclude the Higgs (`/ h`) so signal is not counted twice.

## Common reconstruction

Both detector profiles use the same particle-level reconstruction and cuts:

- assign each final-state photon within `DeltaR<0.1` to its nearest electron
  or muon and form dressed leptons;
- require at least four electrons or muons with `pT>10 GeV` and
  `|eta|<2.5`, vetoing electrons in `1.01<|eta|<1.16`;
- require the stable-particle isolation proxy
  `sum ET(DeltaR<0.35)<5 GeV`, using constituents above 0.5 GeV and excluding
  the candidate, its dressing photons, and neutrinos;
- consider every disjoint OSSF pairing, define `Z1` as the pair closest to
  the Z mass, and require `70<mZ1<100 GeV` and `10<mZ2<100 GeV`;
- require every OSSF pair in the selected quadruplet to have
  `mll>4 GeV`;
- resolve extra-lepton ambiguities by `|mZ1-mZ|`, then decreasing quadruplet
  scalar `pT` sum, then stable source indices.

The analyzer does not cut on `m4l`.  The report's default signal region is the
open interval `120<m4l<130 GeV`.

## Response profiles

`perfect` keeps the common dressing, isolation, pairing, and fiducial
selection but uses identity momentum response and unit reconstruction and
trigger efficiencies.

`ssc` uses the GEM response recorded in the manifest:

- electron efficiency 0.90; sampling terms `0.060/sqrt(E)` in the barrel and
  `0.085/sqrt(E)` outside it, a 0.004 constant term, and the GEM thermal and
  pileup noise terms;
- muon efficiency `0.85*0.95=0.8075`; deterministic curvature smearing from
  the digitized GEM TDR table with linear `|eta|` and logarithmic-`pT`
  interpolation, clamped below 10 and above 500 GeV;
- trigger efficiency 0.98 for reconstructed `4e` and 0.99 for `4mu` and
  `2e2mu`;
- `--muon-resolution-scale 0.8` and `1.2` provide the declared down/up
  resolution variations around the default 1.0.

Charge misidentification, explicit pileup-event overlay, conversions, and
impact-parameter response are deliberately unsupported future effects.  The
noise parameterization is not an event-overlay model.

## Manifest and provenance

Every response run writes a schema-versioned manifest containing the
configured process/decay cards, target and native cross sections, signed
pre-decay LHE-handler denominator, detector constants, actual stage-specific
seeds, generation cuts,
software versions, and SHA-256 checksums for the runner, Herwig template,
analyzer source, cards, input LHE, analyzer summary, and output ROOT file.
Unavailable products in a dry run are marked pending rather than represented
as validated artifacts.
