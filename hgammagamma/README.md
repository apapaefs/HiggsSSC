# Higgs To Gamma Gamma Study

This directory is for the Higgs boson decay study

```text
h -> gamma gamma
```

using the `loop_sm_haa` MadGraph5_aMC@NLO model in this repository.

## Physics Target

The goal is to generate and study Higgs boson decays to two photons, using an
effective Higgs-photon-photon coupling added to the Standard Model loop model.

In MG5 notation the process is

```text
h > a a
```

where `h` is the Higgs boson and `a` is the photon.

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

From the repository root:

```bash
cd MG5_aMC_v3_5_15
./bin/mg5_aMC
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

For a production-and-decay study, start from a Higgs production mode and attach
the decay. For example:

```text
import model loop_sm_haa
generate g g > h, h > a a
output gg_h_gammagamma
```

## Notes To Record

When adding results to this directory, record:

- the MG5 version;
- the exact process card;
- the Higgs mass and width used in the param card;
- any generation cuts;
- the number of generated events;
- the random seed;
- whether events are parton-level only or showered/hadronized.

## Suggested Outputs

Useful first plots or checks include:

- photon transverse momentum;
- photon pseudorapidity;
- diphoton invariant mass;
- photon separation;
- total generated rate or partial width, depending on the setup.
