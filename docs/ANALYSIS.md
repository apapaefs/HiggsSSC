# Cuts and classifiers

For four leptons, use the [HO analysis and run guide](../hfourlepton/HIGHER_ORDER.md).
Both `fourlepton-lo` and `fourlepton-ho` support `cuts` and `xgboost`; HO also
provides the energy-specific `validate` action. The sections below describe
the diphoton interfaces.

Start with a completed diphoton campaign and inspect its ordinary report.
The default cards require all four samples and SSC detector-response metadata;
missing backgrounds must not silently improve a significance estimate.
Single-photon plots include one-candidate outcomes; use diphoton observables
such as `m_gg` when comparing with the selected pair rate.

## Rectangular cuts

```bash
python3 run_pipeline.py cuts --workflow gammagamma-lo --run-tag lo_10k_01
```

The baseline requires two selected photon candidates and
`120 <= m_gg <= 130 GeV`. Outputs are `summary.csv`, `summary.json` and
`index.html` in `hgammagamma/LOAnalysis/analyses/lo_10k_01/cuts/`.
For HO use `--workflow gammagamma-ho` and its completed run tag.

Copy the appropriate [LO](../hgammagamma/analysis_cards/baseline_cuts.yaml) or
[HO](../hgammagamma/analysis_cards/ho_baseline_cuts.yaml) card to make a new
selection, then pass `--config your_card.yaml`. The pipeline supplies campaign
paths, tag and luminosity; these do not need editing in the card. Available
variables are `m_gg`, `pt_gamma1`, `eta_gamma1`, `pt_gamma2`, `eta_gamma2`,
`deltaR_gg`, `deltaPhi_gg`, `pt_gg`, `y_gg`, `n_selected_photons`.

## Optional XGBoost

Install the additional packages once, then use a sufficiently populated run:

```bash
python3 run_pipeline.py setup --workflow gammagamma-lo --with-ml
python3 run_pipeline.py xgboost --workflow gammagamma-lo --run-tag lo_10k_01
```

The LO classifier is exploratory: threshold selection and reported yields do
not constitute an independent discovery test. Compare against the cut result,
inspect overtraining and mass sculpting, and include background systematics
before drawing conclusions. Thirty-event pilots cannot establish classifier
performance.

HO uses the [HO classifier card](../hgammagamma/analysis_cards/ho_xgboost_baseline.yaml).
It groups every response hypothesis from the same source event into one
partition: 60% training, 20% validation, 20% independent testing. Thresholds are
chosen on validation and frozen before test evaluation. Signed weights remain
in physical yields; absolute weights used by training never replace them.
Insufficient effective background statistics produce an unsupported result,
not a discovery claim. Event caps are rejected for HO yields.

## Normalization and older campaigns

LO defaults to `K(ggH)=2` and `BR(H -> gamma gamma)=0.00227`; this is a teaching
rate model. HO requires its validated ihixs sidecar, contains the physical BR
once and accepts no LO K-factor. Detector fake probabilities are already in
the event weights. For luminosity `L` in fb⁻¹, `N = 1000 L sigma_selected` with
`sigma_selected` in pb. `S/sqrt(B)` is only a counting diagnostic.

The direct CLI supports existing campaign layouts without moving them:

```bash
python3 analyze_lo_varfiles.py cuts \
  --config hgammagamma/analysis_cards/baseline_cuts.yaml \
  --analysis-root hgammagamma/LOAnalysis \
  --run-tag YOUR_EXISTING_TAG --output-dir results/your_cuts \
  --luminosity-fb 100
```

Run it in a Python environment with PyROOT and
[the basic requirements](../requirements-analysis.txt); classifiers additionally
need [ML requirements](../requirements-ml.txt). CLI path and luminosity options
override card values. Full trees and original HO provenance files must remain
available for normalization validation.
