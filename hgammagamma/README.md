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

## Download MadGraph

Work in your home directory and download MadGraph5_aMC@NLO version `v3.5.15`
from

```text
https://launchpad.net/mg5amcnlo
```

For example:

```bash
cd ~
wget https://launchpad.net/mg5amcnlo/3.0/3.7.x/+download/LTS_MG5aMC_v3.5.15.tgz
```

This downloads the tarball

```text
LTS_MG5aMC_v3.5.15.tgz
```

in your home directory. Unpack it with:

```bash
tar -xzf LTS_MG5aMC_v3.5.15.tgz
```

You should now have:

```text
~/MG5_aMC_v3_5_15
```

If your copy of this repository is also in your home directory, copy the
`loop_sm_haa` model into the MadGraph model directory:

```bash
cp -r ~/HiggsSSC/MG5_aMC_v3_5_15/models/loop_sm_haa ~/MG5_aMC_v3_5_15/models/
```

If you cloned this repository somewhere else, change the first path to point to
your `HiggsSSC` checkout.

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

For a production-and-decay study, generate the Higgs production process first
and let MadSpin handle the decay. For gluon-fusion Higgs production use:

```text
import model loop_sm_haa
generate g g > h [noborn=QCD]
output gg_h_gammagamma
```

The `[noborn=QCD]` option tells MG5 to use the loop-induced gluon-fusion
process rather than looking for a QCD Born contribution.

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

## Suggested Outputs

Useful first plots or checks include:

- photon transverse momentum;
- photon pseudorapidity;
- diphoton invariant mass;
- photon separation;
- total generated rate or partial width, depending on the setup.
