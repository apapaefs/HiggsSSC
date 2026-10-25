# Four-lepton campaigns: 13 TeV validation and 40 TeV SSC

The workflows are prepared for you to run. No new event generation, shower,
ihixs calculation, classifier training or numerical LHC comparison was run
while implementing this extension. Local synthetic checks test the code;
they do not establish agreement with an experiment. Start with the pilots
below before committing resources to larger samples.

## Prepare and launch

Run these commands from the repository root. `prepare` writes review cards
under `hfourlepton/plans/TAG/` and invokes no simulation programs. The managed
`run` action generates its executable cards from the same configuration; editing
a preview card does not change that action.

```bash
python3 run_pipeline.py prepare --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_pilot --nevents 30
python3 run_pipeline.py prepare --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_pilot --nevents 30
```

Install/check the same higher-order environment used by the diphoton channel,
including the optional analysis dependencies:

```bash
python3 run_pipeline.py setup --workflow fourlepton-ho --with-ml
python3 run_pipeline.py check --workflow fourlepton-ho --with-ml
```

Setup defaults to Timur's `herwig/stable` module. For another installation,
use the environment and compiler options in the [main guide](../README.md).
Launch each pilot yourself, preferably inside a persistent screen session:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_pilot --nevents 30 --jobs 4
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_pilot --nevents 30 --jobs 4
```

Each energy needs its own pilot evidence. Even these small pilots include
costly integration and a full 109-point ihixs calculation. They check execution
and bookkeeping, not physics sensitivity. The 13 TeV run produces one
`perfect` response manifest with CMS-like particle selection; the 40 TeV run
produces `perfect` and `ssc` manifests with identical SSC selection cuts.

Inspect `hfourlepton/runs/TAG/report/index.html`. The report checks normalization,
response probabilities, optional weights, event losses and the saved products.
Use `status` with the same workflow, energy and tag to inspect progress. Repeat
the original `run` command to resume completed stages. Existing mismatched or
partial outputs cause an error rather than automatic replacement.

An older runner can stop at the first ggF event with a `status-1 Higgs mass`
error because it applies an on-shell check to POWHEG's finite-width sample.
After updating the code, preserve the failed four-lepton manifest as a backup
only after confirming that no four-lepton production, decay, shower or analysis
outputs exist. The next invocation can then create that stage's manifest with
the corrected code hashes. Keep the completed POWHEG directory, ihixs record,
managed campaign state and failure log, and repeat the original run command.

After successful pilots, choose fresh tags and sufficient statistics:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_10k --nevents 10000 --jobs 4
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_10k --nevents 10000 --jobs 4
```

Ten thousand events per decay sample is a starting point, not a guarantee of
adequate selected-background statistics. The continuum generation is broad;
inspect effective counts before choosing the next sample sizes.

## Samples and rates

There are 27 samples from 13 hard-process integrations. Every row below is
generated for `4e`, `4mu` and `2e2mu` Higgs decays or continuum final states.

| Contribution | Hard process and shower | Inclusive normalization |
| --- | --- | --- |
| ggF | POWHEG HJMiNNLO, Herwig; pure HEFT | Energy-specific ihixs N3LO |
| VBF | MG5 MC@NLO, Herwig; s-channel W/Z excluded | Native NLO |
| W+H, W−H | MG5 MC@NLO, Herwig | Native NLO, both charges |
| qqZH | MG5 MC@NLO, Herwig | Native NLO |
| ggZH | MG5 loop-induced LO, Herwig | Native LO, separate from qqZH |
| ttH | MG5 MC@NLO, Herwig | Native NLO |
| qq→4ℓ | Full neutral-current matrix elements, Higgs excluded; MC@NLO | Native NLO |
| gg→4ℓ | Continuum loops, Higgs excluded; Herwig | Native loop-induced LO |

The ggF calculation follows the diphoton convention, including the full ihixs
scale/PDF scan and its numerical checks. Its hard PDF is NNPDF4.0 NNLO with QED;
the target rate uses the diphoton N3LO PDF convention. Other samples use
`NNPDF40_nlo_as_01180_qed` (LHAPDF 335900), with native NLO scales set to
MG5's dynamical choice 3, HT/2. Loop samples are explicitly LO despite their
NLO PDFs. No common correction factor is assigned to gg→4ℓ or non-ggF signals.
Electroweak and finite-quark-mass corrections to the ihixs rate are disabled,
as in the existing pure-HEFT diphoton calculation.

All HO Higgs samples use the pole mass mH=125 GeV. HJMiNNLO samples the
event mass within its recorded Breit-Wigner window: the default hwidth is
0.004152 GeV and bwcutoff is 1000, giving 120.848–129.152 GeV. The
`higgsfixedwidth` switch fixes the width prescription, not the event mass.
Validation checks the pole inputs against ihixs and each generated Higgs mass
against those bounds; MG5 stable-H inputs still require an on-shell Higgs.
MadSpin uses a direct scalar 1→4 decay
matrix element. Associated W/Z/top particles decay inclusively in Herwig.
The target branching fraction is applied once through

```text
analysis weight scale = sigma(production) × BR(H→channel) / sigma(MadSpin LHE)
selected sigma = sigma(MadSpin LHE) × sum(selected analysis weights)
                 / sum(consumed LHE handler weights)
```

For ggF, `sigma(production)` is the ihixs value. The common total branching
fraction is `1.251e-4`; the existing channel fractions were calculated at
125.09 GeV and are used as an explicit approximation at 125 GeV. Their original
survey uncertainty and provenance remain in the manifest; see the
[decay validation notes](validation/README.md).

With N requested events per decay channel, each signal production group needs
3N stable Higgs events. Three nonoverlapping source ranges feed the independent
decay samples; integration grids are shared. Source ranges, hashes and the
resulting decay inputs are recorded. This avoids cross-channel train/test
leakage even when shower events are discarded. The response hypotheses of one
event always belong to one statistical group. `perfect` and `ssc` outputs are
alternative responses and must never be added together.

The shower follows the diphoton matching settings, respects SCALUP, reads a
finite LHE stream without replay, and keeps signed weights. Failed shower
events contribute zero response but remain in the consumed-source denominator.
This prevents a loss of events from increasing the quoted efficiency. Large
losses still require investigation; accounting for them does not validate the
shower. All finite optional weights must survive the decay boundary.

Compact LHE `<weights>` values are matched to IDs in `<initrwgt>` declaration
order. MG5 3.5.15 MadSpin rescales explicit `<rwgt>` weights but can copy a
compact vector without the decay factor. The runner prepares a separate,
profile-independent file under `weight_inputs/SAMPLE/` with explicit IDs.
It checks the nominal decay factor, source event order and production momenta,
and corrects only compact values verified to equal their original source values.
Already scaled values are not scaled twice. Original LHE files remain unchanged;
the preparation record includes input/output hashes and conversion counts.

The continuum generation cuts, pT(ℓ)>2 GeV, |η(ℓ)|<5 and dilepton mass>2 GeV,
are below the analysis cuts. Their stability must be checked with looser
generation in a future systematic study. Zbb, ttZ and ttbar remain separate
[reducible-background controls](README.md#reducible-backgrounds-and-decay-bias).
The nominal cards do not include data-driven light-jet fakes, bbH/tH,
electroweak ZZjj, tribosons or signal–continuum interference. The resulting
counting sensitivity therefore covers the stated irreducible model.

For separate reducible-control pilots, the existing LO runner supplies the
unbiased Zbb, ttZ and ttbar definitions. First extend the environment with
`python3 run_pipeline.py setup --workflow fourlepton-lo` to install its LO PDF.
Then, in that activated runtime, the following commands launch the controls:

```bash
python3 hfourlepton/run_four_lepton_campaign.py --order lo --sample-set reducible \
  --ebeam 6500 --selection-profile cms --detector-response perfect \
  --run-tag fourl_lhc13_controls --nevents 30
python3 hfourlepton/run_four_lepton_campaign.py --order lo --sample-set reducible \
  --ebeam 20000 --selection-profile ssc --detector-response ssc \
  --run-tag fourl_ssc40_controls --nevents 30
```

Supply `--mg5-dir`, `--mg5-python` and the Herwig environment options for your
installation, as in the LO guide. These tiny control pilots only check the
chain: rare fourth-lepton tails need much larger unbiased samples or validated
importance sampling. They are not automatically combined with the HO totals.

## Cuts and detector response

The 13 TeV selection follows the fiducial definitions in Tables 7–8 of
[CMS HIG-19-001](https://arxiv.org/abs/2103.04956). The SSC version uses the same
simple pairing/mass topology with the GEM acceptance and isolation.

| Quantity | 13 TeV particle comparison | 40 TeV SSC selection |
| --- | --- | --- |
| Individual lepton pT [GeV] | electrons >7; muons >5 | all >10 |
| Absolute pseudorapidity | electrons <2.5; muons <2.4 | all <2.5 |
| Two hardest lepton pT [GeV] | >20, >10 | >20, >10 |
| Electron crack | none | exclude 1.01<|η|<1.16 |
| Isolation | scalar pT in R=0.3 below 0.35 pT(ℓ) | scalar ET in R=0.35 below 5 GeV |
| Candidate | two disjoint OSSF pairs; nearest-Z first | same |
| Pair masses [GeV] | 40<mZ1<120; 12<mZ2<120 | same |
| Other pairs | every OS mass >4 GeV; every ΔR>0.02 | same |
| Analysis range [GeV] | 105<m4ℓ<140 | same |
| Baseline counting window [GeV] | 120<m4ℓ<130 | same |

Photons within ΔR=0.1 dress the nearest lepton. The original LO `gem` selection
remains available for comparison. The new `ssc` selection name denotes the
cuts, while `--detector-response ssc` independently enables smearing and
efficiencies.

Electron energy resolution uses the GEM sampling terms 0.060/√E in the barrel
and 0.085/√E in the endcap, a 0.004 constant term and the existing thermal and
pileup noise terms in quadrature. Muons are smeared in signed curvature q/pT
using the approximate digitization of GEM TDR Fig. 4-20, interpolated linearly
in |η| and logarithmically in pT. Values outside 10–500 GeV use the endpoint
resolution and are counted in diagnostics. At 500 GeV the central and
|η|=2.5 fractional resolutions are 5% and 12%. Electron efficiency is 0.90;
muon efficiency is 0.85×0.95. Directions are unchanged. These are a simplified
GEM response, not CMS or ATLAS detector calibrations.

Sources: [GEM TDR](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-262.pdf),
[SSC Higgs study](https://lss.fnal.gov/archive/other/calt-68-1856.pdf), and the
existing [detector implementation notes](../hgammagamma/SSC_DETECTOR_RESPONSE.md).

## Cut analysis, XGBoost and LHC comparison

After completing the larger campaigns, launch the analyses yourself:

```bash
python3 run_pipeline.py cuts --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_10k --luminosity-fb 100
python3 run_pipeline.py xgboost --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_10k --luminosity-fb 100
python3 run_pipeline.py validate --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_10k
```

Use `--energy-tev 13` and the LHC tag for its cut/XGBoost analyses. Both modes
also support legacy `fourlepton-lo` samples. The
[analysis cards](analysis_cards/) define the sample inventory and cuts.
Copy a card and pass `--config PATH` to change rectangular cuts or training
settings. To repeat an analysis with changed settings, call the
[`analyze_fourlepton.py`](analyze_fourlepton.py) CLI with a fresh `--output-dir`;
nonempty analysis outputs are protected.

The nine classifier inputs are mZ1, mZ2, pT(4ℓ), |y(4ℓ)|, the four ordered
lepton transverse momenta and minimum pair ΔR. There are no jet categories,
matrix-element discriminants or m4ℓ training input. Training uses the broad
mass range; threshold selection and the cut/classifier comparison use the
same baseline mass window. Splits are 60% training, 20% validation and 20%
test, grouped before event selection. Class-balanced absolute weights are
used only for fitting. Physical yields, variances and threshold selection
use signed weights. The threshold is fixed on validation before the test
sample is evaluated. A minimum of 25 effective background events is required
at the chosen threshold. The default 10% background uncertainty is an
illustrative nuisance assumption, not an experimental estimate.

Outputs under `hfourlepton/runs/TAG/cuts/` and `xgboost/` include HTML,
JSON/CSV yields and mass plots. XGBoost also writes the model, per-event
partition/score CSV, training/test score plots and feature importances.
Inspect mass sculpting even though m4ℓ is absent from training.

`validate` writes `lhc_validation/comparison.json` and `index.html`, comparing
mode-by-mode fiducial acceptances with CMS Table 8. It requires at least 1000
selected effective events per production mode. Its 10% relative acceptance
tolerance is a broad implementation check, not the CMS uncertainty. Too few
events or a disagreement produces `not_validated` and exit status 2. No
normalization is tuned to the published measurement.

This is an approximate particle-level benchmark. HwSim lacks prompt/hadron
ancestry and stores particles above 0.4 GeV, so dressing, isolation and
associated-production combinatorics differ from the experiment. CMS detector
ID, impact-parameter and vertex selection are not emulated. An
[ATLAS four-lepton measurement](https://arxiv.org/abs/2004.03969) also uses
different kinematic thresholds, a matrix-element choice with extra leptons,
and lepton–jet overlap rules. No separate numerical ATLAS comparison or
experimental likelihood is implemented. Agreement in this CMS benchmark
would support the basic selection and acceptance implementation; it would
not validate the historical SSC resolution or establish discovery reach.

## Local implementation checks

These commands use synthetic fixtures and do not generate physics samples
or fit a classifier:

```bash
make -C hfourlepton/LOAnalysis/Code check
python3 -m unittest discover -s tests -p test_fourlepton_ho.py
python3 -m unittest discover -s tests -p test_fourlepton_campaign.py
python3 -m unittest discover -s tests -p test_fourlepton_report.py
python3 -m unittest discover -s tests -p test_ihixs_normalization.py
python3 -m unittest discover -s tests -p test_student_pipeline.py
```
