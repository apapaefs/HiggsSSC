# Four-Lepton Generator Validation

This directory records small, generator-level checks used by the
`hfourlepton` campaign.  It contains command cards, not generated events.

## Signal channel fractions

The default channel fractions were obtained with MG5_aMC 3.5.15 and the
restricted Standard Model by integrating the direct tree-level decays

```text
h > e+ e- e+ e-
h > e+ e- mu+ mu-
```

at `mH = 125.09 GeV`.  The reference survey returned

```text
Gamma(4e)       = 1.353e-7 GeV
Gamma(2e2mu)    = 2.436e-7 GeV
Gamma(4mu)      = Gamma(4e)
```

and hence

```text
f(4e)       = 0.263127188
f(4mu)      = 0.263127188
f(2e2mu)    = 0.473745624
```

The survey uncertainty is at the percent level.  These fractions set the
relative normalization of the three LO decay samples; their sum is normalized
to the external `BR(H -> 4l) = 1.251e-4`.  They are not presented as a
replacement for a PROPHECY4F calculation.

Run the cards from a disposable working directory.  Generated process
directories and events must not be committed.

## Interference validation

The nominal campaign keeps Higgs signal and the loop-induced gluon continuum
separate.  The `gg_4l_coherent.mg5` and `gg_4l_continuum.mg5` cards provide a
small validation pair for estimating the omitted interference near the Higgs
peak.  They are intentionally not part of the default background sum.

## LHE forced-decay weighting

The custom LHE handlers were exercised with ten stable-top events through
Herwig 7.3 and HwSim.  A prompt-only two-event check returned the same
per-event factor in both events,

```text
(BR(W -> e nu or mu nu))^2 = 0.0453183201610.
```

The handler is run after hadronization, before weak hadron decays.  This
timing was checked explicitly against an event containing a virtual `W-`
phase-space line from a `Lambda_b` decay: that line receives no additional
prompt factor.

The final full-support semileptonic proposal was also read successfully by
Herwig with the complete `HerwigBDecays.in` mode table. A one-event `ttbar`
smoke run returned a finite total event factor `0.0286571568`; the same
decay-importance factor was propagated consistently to all 145 optional LHE
weights. The heavy-flavour-scoped handler runs after recursive decays and
counts terminal `B0`, `B+`, and `Bs` copies, including those produced by
excited states while suppressing same-particle and neutral-meson-mixing
ancestors. This validates handler placement and weight propagation, not
physics closure: the accelerated samples remain marked
`closure_validated: false` and excluded from nominal totals.
