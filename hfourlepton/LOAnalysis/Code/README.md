# Four-lepton HwSim analyzer

Build the analyzer with:

```bash
make
```

Run the committed core regression test with:

```bash
make check
```

The core test directly checks the GEM muon-resolution anchors, clamping and
log-pT interpolation; deterministic electron and muon smearing; nearest-photon
dressing and isolation exclusions; OSSF pairing and channel assignment; and
exclusive response-subset probability closure.

The mandatory input is a `Data` tree containing `numparticles/I`,
`objects[8][10000]/D`, and `evweight/D`. The analyzer also consumes
`numJets` with `theJets[5][100]`, `numbJets` with `thebJets[5][100]`,
`theETmiss[4]`, and `theOptWeights` with `theOptWeightsNames` when those
branches are present. All optional groups may be absent.

The primary command-line contract is:

```bash
./HwSimPostAnalysis_fourlepton INPUT.root \
  --response-profile perfect \
  --seed 14101983 \
  --weight-scale 1 \
  --muon-resolution-scale 1 \
  --tag smoke
```

`--input-list FILE.input` is an alternative to the positional input.
`--detector-response` aliases `--response-profile`, and `--run-tag` aliases
`--tag`. Campaign calls may additionally pass `--output-dir`, `--sample`,
`--category`, `--channel`, and `--optional-weight-names-file`. The last
option restores ordered semantic LHE variation names when HwSim writes the
weight values but leaves its name branch empty; populated HwSim names are
cross-checked rather than overwritten.

With `--output-dir DIR --sample SAMPLE --tag TAG`, outputs are
`DIR/SAMPLE_TAG_PROFILE.root` and
`DIR/SAMPLE_TAG_PROFILE.summary.json`. Without `--sample`, the input stem is
used instead.

The ROOT output contains:

- `AnalysisMetadata`: schema, response, channel, and cut-mask definitions.
- `SourceEvents`: one entry per input event, including raw and scaled weights,
  optional weights, reconstruction-subset closure, selected probability,
  post-baseline accepted-lepton count, and source-level cut probabilities.
- `FourLepton`: every nonzero-probability reconstructed subset containing at
  least four leptons, including subsets that fail the physics cuts.

Lepton reconstruction subsets are formed before the analysis cuts, so the
cutflow stages are meaningful rather than prefiltered constants. The
cumulative cut-mask bits are: 0 at least four reconstructed leptons, 1
kinematic acceptance, 2 electron crack acceptance, 3 isolation, 4 disjoint
OSSF pairs, 5 all OSSF masses above 4 GeV, 6 the Z1 window, 7 the Z2 window, 8
trigger weight applied, and 9 final selection. Channels are 1 for `4e`, 2 for
`4mu`, 3 for `2e2mu`, and 0 when no selected topology can be classified.

The event-weight invariant is
`event_weight = generator_weight * sample_weight * response_weight`.
`generator_weight` is the raw HwSim weight and is not interpreted as a cross
section.

Enumeration is exact for as many as ten baseline-accepted source leptons. If
an event contains additional leptons that deterministically fail the baseline,
their reconstruction states are analytically marginalized; up to ten
candidates (all baseline-passing candidates first) are retained to represent
the staged cutflow. A run containing an event with more than ten
baseline-accepted leptons fails closed rather than silently truncating its
response space.

## Synthetic smoke input

Build the small fixture generator and run both the full optional schema and
the required-branches-only schema:

```bash
make fixture
./make_synthetic_hw_sim_fourlepton /tmp/fourlepton-full.root
./make_synthetic_hw_sim_fourlepton /tmp/fourlepton-minimal.root minimal
./HwSimPostAnalysis_fourlepton /tmp/fourlepton-full.root \
  --response-profile perfect --output-dir /tmp --sample synthetic --tag perfect
./HwSimPostAnalysis_fourlepton /tmp/fourlepton-minimal.root \
  --response-profile ssc --output-dir /tmp --sample synthetic --tag minimal
```

The fixture contains selected `2e2mu`, `4e` with an extra lepton, and
negative-weight `4mu` events; it also exercises photon dressing, optional
weights, the strict pT boundary, the electron crack, isolation, missing OSSF
pairs, and a failing Z1 window.
