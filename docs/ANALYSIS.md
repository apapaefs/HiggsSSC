# Cuts and classifiers

For four leptons, use the [HO analysis and run guide](../hfourlepton/HIGHER_ORDER.md).
Both `fourlepton-lo` and `fourlepton-ho` support `cuts` and `xgboost`; HO also
provides the energy-specific `validate` action. The sections below describe
the diphoton interfaces.

Start with a completed diphoton campaign and inspect its ordinary report.
The default 40 TeV higher-order diphoton cards require six samples with the
SSC response. The signal, prompt diphoton, photon + jet and Drell–Yan samples
come from the campaign. The loop-induced diphoton and dijet double-fake
distributions come from the [weighted event tables](../hgammagamma/HOAnalysis/student_backgrounds/README.md).
Across the six samples, single-photon plots count each selected photon,
and `m_gg` counts each selected pair. The dijet table contains two-fake pairs
and uses the upper estimate from independent quark-to-photon fake rates.

## Additional backgrounds at 40 TeV

At 40 TeV the additional samples cover the loop-induced diphoton continuum,
nonprompt leptons in heavy-flavour events, and QCD dijets with two
photon fakes. The signal uses its existing iHixs normalization. Rare
prompt-photon processes are outside the chosen scope. The standard diphoton
classifier now trains and selects its validation threshold using all five
backgrounds. The separate studies below evaluate the original frozen threshold
and report the effect of each added background on that earlier selection.

### Loop-induced diphoton continuum

[run_loop_gammagamma.py](../hgammagamma/run_loop_gammagamma.py) generates a
continuum-only gluon-fusion sample at loop-induced LO and showers it with
Herwig. The process excludes the Higgs amplitude. The native generator rate
enters once, with no branching-ratio or K-factor rescaling. The
[evaluator](../hgammagamma/evaluate_loop_gammagamma.py) applies the original
mass window and the frozen classifier to the SSC photon response.
See the [HO guide](../hgammagamma/HOAnalysis/README.md) for the pilot inputs
and the [MadGraph paper](https://arxiv.org/abs/1507.00020) for the generation
method.

### Four-lepton reducible controls

The unbiased LO samples are top-pair production with both W decays to a
direct electron or muon, and bottom-pair production with a dilepton
neutral current in each of the electron and muon channels. Each sample uses
the SSC response and its own cross section and decay normalization.
[run_weighted_reducible_pilot.py](../hfourlepton/run_weighted_reducible_pilot.py)
tests the efficiency of enhancing semileptonic heavy-hadron decays.
It compares disjoint parent-event ranges with weighted and ordinary native
Herwig decays. A third, disjoint range showered with EvtGen tests the change
in decay model. Report the rate and the source-event MC uncertainty in the
isolation sidebands as well as in the mass window and after the classifier.

The pilot used 100,000 hard events per process in each of the two native
Herwig decay samples. All 60 showers completed. After applying the decay
weights, inclusive rates agreed within 0.14%. In the broad region with two
of four leptons failing isolation, weighting gave 0.83, 0.90 and 0.97 times
the ordinary-decay effective events per worker-hour for top pairs,
electron-channel bottom pairs and muon-channel bottom pairs, respectively.
Every mass-window and classifier yield was zero. The three-million-event
unbiased sample therefore supplies the current joint 95% upper limit of
0.546 fb under the simulated response, counting failed showers as zero.
The disjoint native-versus-EvtGen comparison shows 14% and 17% lower
native fiducial rates in the top and muon-channel bottom samples. The final
selection has too few events to measure that dependence.

These samples model nonprompt leptons from heavy-flavour decays. Hard
charm-pair production with a Z, generic jet-to-lepton misidentification and
material conversions require other inputs. The control regions used by CMS
illustrate how to check the isolation extrapolation. [CMS, section 7.2](https://arxiv.org/pdf/2103.04956)

### Dijet double fakes

The [SSC response guide](../hgammagamma/SSC_DETECTOR_RESPONSE.md#dijet-double-fake-bounds)
describes the pair weights, source normalization and the limits of transferring
GEM jet-to-photon probabilities to dijets. Evaluate two selections independently.
The cut analysis uses 120–130 GeV. The classifier requires a score of at least
0.992 over the accepted mass range. The hard-dijet samples with and without
multiple parton interactions and the independently generated inclusive
minimum-bias samples describe alternative models of the same contribution.

The inclusive Herwig sample at 40 TeV has 1,320,000 independent collisions
split across 132 seeds. Under independent per-jet fakes, the table shows the
model endpoints in pb. Quoted errors are Monte Carlo standard errors across
seeds. The lower endpoint uses only events with exactly two eligible jets,
both assigned the gluon fake probability. The upper endpoint sums over every
eligible pair using the quark probability. A pair of fake photons is
constructed from the two EM-smeared jet energies and directions.

| Selection | Original background [pb] | Lower dijet endpoint [pb] | Upper dijet endpoint [pb] |
| --- | ---: | ---: | ---: |
| Diphoton mass 120–130 GeV | 5.423 | 0.113 ± 0.008 | 3.011 ± 0.162 |
| Frozen classifier, score ≥ 0.992 | 0.507 | 0.0052 ± 0.0018 | 0.138 ± 0.039 |

The upper endpoints have 319 and 15 contributing collisions, respectively.
If only the single-jet fake probabilities are specified, the joint fake
probability can reach the smaller of the two marginals. In this model the
corresponding classifier envelope rises to 377 ± 105 pb. The
independent-fake interval therefore depends critically on the joint response
of two jets. These endpoints describe the sampled generator and detector
models. They carry no confidence-level interpretation.

For historical comparison, GEM Table 14 gives 60 fb/GeV from two fake
photons at 120 GeV before the angular cut, and 16 fb/GeV after it.
At 140 GeV these rates are 34 and 6.8 fb/GeV. GEM used flavour-dependent jet
rejection, a jet-background normalization factor of 1.5 and a final
requirement of |cos θ*| < 0.7. The present hard-dijet samples use their LO
rates. Compare rates per GeV at a specified selection. [GEM Higgs study, section 4.3](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf)

### Conditional discovery sensitivity

The loop rate combines the 100–150 GeV Born-mass contribution from a
5,000-event focused sample with the complementary region of an independent
20,000-event broad sample. Matching the two hard photons back to the event
record establishes the generated mass for every reconstructed event.
All 16 mass-window and three classifier selections in the broad sample came
from the central Born interval. The combined estimates are 0.739 ± 0.044 pb
for the 120–130 GeV window and 0.0899 ± 0.0158 pb for the classifier. The
one-sided 95% finite-MC upper rates are 1.320 and 0.376 pb, respectively,
allowing for an unsampled contribution outside the central interval. These
bounds hold for each selection separately with the generator rates fixed.
One broad shower had a momentum residual of about 3 TeV. Its hard diphoton
mass was 70.32 GeV and it yielded at most one reconstructed photon. Its
observed contribution to both selections is zero. Giving that event a unit
response would add 0.0559 pb to either selection as a single-event check.
This was the only broad event above the 0.4 GeV threshold of the Herwig
BasicConsistency check. The residual appeared during the shower. Its
underlying cause is unresolved. The original signal and photon+jet
classifier reports also record shower residuals above 4 TeV. The
conditional calculation uses the saved nominal yields for these samples.

The two-channel calculation uses the frozen classifiers at 0.992 in
diphotons and 0.488 in four leptons. Treat their event counts as independent
Poisson variables. At fixed signal and background rates, order the outcomes
by their likelihood ratio for SM signal plus background versus background
alone. The entries below are the first luminosities at which the probability
of crossing a one-sided local 5σ threshold reaches 50%. The diphoton signal
is 0.1084 pb and the four-lepton signal is 1.669 fb. The independence,
background normalizations and jet-fake assumptions are model inputs.

| Dijet endpoint | Added four-lepton background [fb] | Combined luminosity [fb⁻¹] | Months at 10 fb⁻¹/year |
| --- | ---: | ---: | ---: |
| Lower | 0 | 1.02 | 1.22 |
| Lower | 0.546 | 1.19 | 1.43 |
| Upper | 0 | 1.17 | 1.41 |
| Upper | 0.546 | 1.41 | 1.69 |

With both backgrounds at their upper endpoints, raising the loop contribution
from its central classifier rate to its 0.376 pb finite-MC upper rate gives
1.84 fb⁻¹, or 2.20 months. Each row is a conditional scenario.

The original diphoton classifier background is 0.507 ± 0.092 pb from 37
selected test sources. A discovery projection with background uncertainties
requires a constraint on the total background, for example from mass
sidebands, and a joint fake-rate model. The calendar conversion assumes
10 fb⁻¹ per year from the beginning of running. For the historical comparison,
GEM reported 7.8σ at 120 GeV in diphotons alone with 10 fb⁻¹, a narrow
mass bin and its photon angular cut. A square-root luminosity extrapolation
gives about 4.1 fb⁻¹ for 5σ at that mass. [GEM, Table 18](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf)

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

### Diphoton mass after cuts and after the classifier

The HO classifier includes the reconstructed `m_gg` as an input. Its only
preselection is two selected photon candidates; the baseline mass window is
used for the separate cut comparison. A score cut can therefore sculpt both
signal and background mass spectra. Its selected signal width is not a
measurement of detector resolution.

The normal student analysis commands generate these plots automatically:

```bash
python3 run_pipeline.py cuts --workflow gammagamma-ho --run-tag YOUR_TAG
python3 run_pipeline.py xgboost --workflow gammagamma-ho --run-tag YOUR_TAG
```

The `cuts/index.html` report includes `m_gg_cut_based.png` and its PDF/CSV
downloads. The `xgboost/index.html` report includes that cut spectrum,
`m_gg_xgboost.png` and `mass_acceptance.png`, each with PDF/CSV downloads,
alongside the existing classifier diagnostics. No separate plotting command
or extra switch is required. If no validation threshold has sufficient
background statistics, the mass plot shows the spectrum before score
selection and the acceptance plot explains why no selection is available.

To replot an older HO classifier report without retraining or changing its
validation-selected threshold, the standalone command is still available:

```bash
python3 hgammagamma/plot_ho_mass_diagnostics.py \
  --xgboost-dir hgammagamma/HOAnalysis/analyses/YOUR_TAG/xgboost \
  --output-dir results/YOUR_TAG_mass
```

The cut plot uses the full sample, applies every baseline cut except `m_gg`,
and overlays the result including the mass window. The XGBoost plot uses only
the independent test events and their saved projection factors. Both include
a 100–150 GeV overview and a 120–130 GeV view in **1 GeV bins**, comparable
to the approximately 0.85–0.9 GeV signal core resolution. The mass acceptance
also uses 1 GeV bins. Physical spectra use signed yields and MC uncertainty
grouped by source event. The standalone command additionally writes an HTML
report and unbinned width measurements. Central quantiles use explicitly labelled absolute weights because
a signed-weight cumulative distribution need not be monotonic; signed
moments and the negative weight fraction are reported alongside them.

New HO `scores.csv` exports contain all ten reconstructed variables. Older
reports require their original `_var.root` files: the tool joins by sample
and tree entry, checks source-event identities, weight closure and the
complete eligible population. For local copies, add
`--input-dir PATH` with files named `signal_gg_h_aa.root`,
`bkg_prompt_aa.root`, `bkg_gamma_j.root`, and `bkg_dy_ee.root`.
The existing 10 GeV campaign mass histogram cannot be rebinned to resolve the
peak; the continuous event-level masses are needed.

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
