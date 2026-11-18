# Additional 40 TeV diphoton backgrounds

The two NPZ tables in `ssc40_v1/` contain reconstructed SSC event
hypotheses at 40 TeV. The manifest lists the feature order, input hashes,
normalization and reference cross sections. Every row has ten variables,
a source-event ID within its sample, and a physical weight in pb.
Sum the weights of the rows passing a selection. For expected counts,
multiply by 1,000 and by the integrated luminosity in fb⁻¹.
Group hypotheses with the same source ID before computing the event-sampling
variance. The source IDs can have gaps.

The student classifier assigns each represented source and all its hypotheses
to one training, validation or test partition. It projects a partition's yield
by its fraction of the represented sources in that sample. Collisions without
a stored contribution already enter the generated exposure used to calculate
each row's physical weight.

| Sample | Generated collisions | Two-candidate rows | 120–130 GeV [pb] | Original classifier at 0.992 [pb] |
| --- | ---: | ---: | ---: | ---: |
| Continuum gg → γγ, loop-induced LO + Herwig | 25,000 | 2,199 | 0.739098 | 0.089928 |
| Two jets faking photons, independent quark-rate upper endpoint | 1,320,000 | 10,629 | 3.010560 | 0.138381 |
| Same dijet collisions, gluon-rate lower endpoint | 1,320,000 | 10,629 available | 0.113334 | 0.005223 |

The loop sample combines broad production outside the generated
100–150 GeV diphoton mass interval with independent focused production
inside it. Each contribution uses its own native cross section divided
by all attempted hard events. It includes the zero-, one- and two-candidate
SSC response hypotheses. An event in the broad sample with a large shower
momentum residual reconstructs at most one photon.

The dijet contribution consists of eligible pairs in 132 equal-size
minimum-bias samples. Each sample has its own generated inelastic cross
section. The default `weight_pb` uses the independent GEM quark fake rate
for each jet. Each event's sum over pairs is below one, so the weights add
for any bin or cut. The alternative `lower_weight_pb` in the same file uses
the gluon fake rate for events containing exactly two eligible jets. These
are alternative assumptions for the same events. Use one column per
prediction. For inclusive single-photon observables, the double-fake sample
contributes both candidates of each pair. It describes no zero- or one-fake
component. The model gives each fake the EM-smeared jet energy.

These reference classifier rates use the frozen classifier of the
original four-sample study. Training a classifier with the added samples
can yield a different threshold and selected rate.

The maintainer exporter checks the original source mapping, verifies the
checksums of the classifier and detector analyzer, and replays both
selections for every minimum-bias seed. It runs on the analysis host with
uproot, NumPy, XGBoost and the existing compiled SSC analyzer. Given the
completed loop-combination and minimum-bias reports, the invocation is

```bash
python3 hgammagamma/package_student_backgrounds.py \
  --loop-combination /path/to/disjoint_loop_combination.json \
  --minbias-report /path/to/mb_132_envelopes.json \
  --minbias-dir /path/to/directory/containing/mb_part_0 \
  --analyzer /path/to/dijet_envelopes \
  --model /path/to/signal_background_xgboost.json \
  --output /path/to/new/ssc40_v1
```
