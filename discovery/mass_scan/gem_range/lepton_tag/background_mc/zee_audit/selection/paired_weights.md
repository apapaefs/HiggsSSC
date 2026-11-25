# Reproducing the paired hard weights and sampling errors

Every collision with a reconstructed candidate in the detector study is matched to its original hard event. The Higgs-window baseline can be cross-checked using the 149 original selected collisions. To include migrations induced by a changed isolation cut, the input to the matching is the full list of candidates **before isolation**. An input row represents a genuine photon, one electron assigned as a photon, and the other electron assigned as a lepton tag. The order of these rows in the compressed detector file defines the order of all probability arrays.

## Hard weights

For each event, the input includes the generator integration, index in the LHE file, parton flavours and momenta, hard photon transverse momentum and scale, ratio of EHLQ1 to NNPDF4.0 for the incident quark densities, and the ratio of the squared production-Zγ and full eeγ matrix elements. The PDF ratio is evaluated at the same factorization scale for both sets. The generator card, archived EHLQ1 source and MG5 current were independently verified in the hard-process audit. [Input conventions and matrix element checks](../hard_process.md)

A collision from integration i in generation region r receives the weight

```text
σ_i / N_r × (|M_prod Zγ|² / |M_full eeγ|²) × (f_EHLQ1 f_EHLQ1 / f_NNPDF f_NNPDF)
    × I(pT_hard_gamma > 10 GeV) × 1.20 × 0.99 .
```

The cross section is in fb and N_r counts all attempted generated events in that phase-space region. The 1.20 factor implements GEM's gluon-fusion normalization. The 0.99 factor converts the simulated 1/33 tracker probability to 3.0%. The other object probabilities, including the trigger, enter each alternative through its detector rows. The core and acceptance-extension regions are disjoint and their rates add. The code keys each collision by its immutable source and LHE index. In shower blocks with lost events the incoming flavours and beam energies identify the generator event. Failed showers contribute zero detector weight.

To regenerate the weights, first extract the core and complementary raw archives and reproduce the shower trees. Build the archived EHLQ1 shared library as in the hard-process audit, and generate `herwig_candidates.json.gz` with the detector reader. Set `work` to the absolute work directory holding these inputs and `MG5` to the path of the MG5 3.5.15 `bin/mg5_aMC` executable. From the repository root in a Python environment with NumPy, LHAPDF and uproot, run the following commands. [Shower replay and input archives](../detector-replay.md) · [EHLQ1 library build](../README.md)

```bash
release=/home/shared/HiggsSSC/student-results/ssc40-tag-background-20261125
audit=discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit
selection="$audit/selection"
features="$selection/herwig_candidates.json.gz"
python3 "$selection/paired_hard.py" collect "$features" "$work/core" "$work/wide" \
    "$release/background.json.gz" "$release/zee-confirmation-spectrum.json.gz" \
    "$work/showers" "$work/generator/libehlq1.so" "$work/paired-momenta.json.gz"
python3 "$selection/paired_hard.py" weight "$work/paired-momenta.json.gz" "$MG5" \
    "$audit/hard_parameters.dat" "$work/matrices" "$work/paired-all.json.gz"
python3 "$selection/paired_hard.py" seal "$work/paired-all.json.gz" "$features" "$work/paired_weights.json.gz"
```

The last command ties the hard input to the SHA-256 of the actual detector file. It also supports computing the hard weights first on a larger kinematic sample and then restricting to all detector candidates. The study used this option to check 677 distinct collisions passing the three-object kinematic preselection. All 149 original selected events reproduce the archived PDF ratio and hard four-momenta. Both the historical-input and original modern rates and their Monte Carlo standard errors reproduce the published baseline.

## Detector variations and paired errors

In Python, add the `selection` directory to the import path. Each variation comprises one unconditional detector probability for each entry in `packet["assignments"]`, in file order. The following example reproduces the nominal rates and their covariance.

```python
import numpy as np
from paired_fit import evaluate, load

packet = load("herwig_candidates.json.gz", "paired_weights.json.gz")
nominal = np.array([r["nominal_probability"] for r in packet["assignments"]])
report = evaluate(packet, {"nominal": nominal})
print(report["results"]["nominal"])
```

The event probability sums mutually exclusive assignments in a collision before forming second moments. Within an integration the covariance is calculated from each collision's vector of alternative weights and all generated zero-weight collisions. The numerical uncertainty on each independent hard cross section is common across alternatives. This yields a paired error on the absolute change, the ratio of two rates and the fraction of the original discrepancy removed, (B_nominal − B_variant)/(B_nominal − 0.0374 fb). The GEM reference and model parameters are held fixed in these sampling errors. The same calculation with unit hard reweighting provides an independent modern-physics result for every alternative.

## How much of the isolation change comes from each object?

For each of the six possible orders, the calculation replaces the scalar 5 GeV isolation decision for one object at a time with the cell-threshold, cluster-subtraction and noise model. It assigns each change in the selected event probability to the object just changed and averages over all six orders. This allocates interactions between three isolation requirements symmetrically. The three signed changes sum to the measured isolation shift for every candidate and for the total selected rate. Errors use the correlations between candidates in a collision and substitutions on the same sample. These attributions apply to the chosen independent-cone noise model.

```bash
python3 "$selection/paired_decomposition.py" "$features" "$work/paired_weights.json.gz" "$work/paired_roles.json"
```
