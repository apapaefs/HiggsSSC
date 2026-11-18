# Higher-order diphoton physics

The [common workflow](../../README.md) installs and runs the complete
pipeline in your account. This reference explains its physics; use
[recovery](../../docs/TROUBLESHOOTING.md) for interrupted jobs.

## Hard processes and shower

| Sample | Process and matching | Hard/shower PDF |
| --- | --- | --- |
| `signal_gg_h_aa` | POWHEG HJMiNNLO; Herwig `h -> gamma gamma` | NNPDF40_nnlo_as_01180_qed, 336100 |
| `bkg_prompt_aa` | `p p > a a [QCD]`, MC@NLO | NNPDF40_nlo_as_01180_qed, 335900 |
| `bkg_gamma_j` | `p p > a j [QCD]`, MC@NLO | same NLO PDF |
| `bkg_dy_ee` | `p p > e+ e- [QCD]`, MC@NLO | same NLO PDF |

Defaults are 20 TeV per proton beam, `mH=125 GeV`, Herwig's angular-ordered
shower, hadronization and installed underlying-event tune. The MPI PDF remains
that of the tune. HJMiNNLO supplies inclusive NNLO+PS Higgs shapes; its
Higgs-plus-jet construction does **not** require an analysis jet. An independent
[ihixs](IHIXS_STUDENT_GUIDE.md) pure-HEFT N3LO calculation fixes the inclusive
signal rate. This changes normalization, not the accuracy of every shape.

Backgrounds use `loop_sm-no_b_mass` and five-flavour parton lists including
`b` and `b~`. Drell–Yan contains gamma*/Z interference. MG5 uses
`parton_shower=HERWIGPP`, `ickkw=0`, `event_norm=average`, and central scales
`muR=muF=HT/2` (`dynamical_scale_choice=3`). Its `--parton` step postpones the
shower while retaining MC@NLO subtraction terms.

The templates follow Herwig 7.3's `LHE-MCatNLO.in` recoil/spin settings and
`LHE-POWHEG.in` SCALUP veto. Both preserve signed weights and prohibit replay
of finite LHE input. No additional truncated-shower implementation is supplied.
Scale-reweight information stays in the LHE; the HwSim response currently uses
the nominal weight rather than the full variation ensemble.

## Generation cuts and omissions

Prompt photons use smooth-cone isolation: `gamma_is_j=False`, `ptgmin=10 GeV`,
`R0gamma=0.4`, `epsgamma=1`, `xn=1`, `isoEM=True`, and `|eta|<6`.
Anti-kt jets have `R=0.4`; a 10 GeV Born-jet threshold applies only to gamma+jet.
Drell–Yan uses `pT(l)>10 GeV` and `m(ll)>30 GeV`.

Check generation-cut and isolation dependence before a physics result. The
smooth-cone definition differs from the detector's jet-based isolation proxy;
photon fragmentation is absent. For the default 40 TeV SSC plots, the four
generator samples above are combined with the published loop-induced
gg → γγ continuum and dijet double-fake event tables. The dijet distribution
uses the upper independent-fake scenario. See the
[student inputs and normalizations](student_backgrounds/README.md).
Signal–background interference is outside the calculation. Gamma+jet
uses one hard fake candidate with the [documented flavour model](../SSC_DETECTOR_RESPONSE.md),
whose ambiguity grows with extra NLO partons.

## Separate loop-induced continuum

[`run_loop_gammagamma.py`](../run_loop_gammagamma.py) adds the continuum
`g g > a a / h [noborn=QCD]` as a separate `bkg_loop_gg_aa` sample. The Higgs
is excluded from the amplitude and checked in the LHE provenance. This is the
one-loop amplitude squared, of order `alpha_s^2 alpha_EM^2`, followed by an
ordinary LO shower. That contribution is absent from the nominal NLO diphoton
calculation. It is not an NLO prediction for the gluon channel, and no signal
or background K-factor is applied.

The model is `loop_sm-no_b_mass`, with massless bottom quarks and a massive
top quark. The NLO PDF and coupling convention is shared with the nominal
continuum. The Born cuts are `pT(gamma)>10 GeV`, `|eta(gamma)|<6`; there is no
outgoing QCD parton at this order. The same SSC reconstruction, isolation and
photon efficiencies are applied after showering. Rare prompt-photon processes
and signal–continuum interference are outside this addition.

Start with a technical smoke test, then use a new tag and seed for a larger
pilot; a successful smoke test alone does not establish an acceptance result:

```bash
python3 hgammagamma/run_loop_gammagamma.py --stage all \
  --run-tag gg_aa_loop_smoke --nevents 100 --nb-core 20 \
  --seed-base 181126031 --ebeam 20000
```

The usual runtime overrides (`--mg5-dir`, `--mg5-python`, compiler paths,
`--no-herwig-module`) also apply. The native MG5 cross section supplies the
rate. Source EOF must be proven without replay, and all attempted source
weights form the normalization denominator, including failed shower events.
Cards, logs, ROOT weight moments and the completion record are checked before
the SSC analysis is accepted.

Use a saved nominal cut summary to add the new contribution without changing
the cuts or nominal rates:

```bash
python3 hgammagamma/evaluate_loop_gammagamma.py \
  --analysis-root hgammagamma/HOAnalysis/runs/gg_aa_loop_pilot \
  --run-tag gg_aa_loop_pilot \
  --nominal-summary results/ho_bg_1m_01/cuts/summary.json \
  --output results/gg_aa_loop_pilot/frozen_cuts.json
```

The output includes the additive rate, finite-source MC error and combined
totals, with hashes of the frozen summary and loop manifest. It also gives a
conservative one-sided 95% upper bound: count hard events with a positive
selected response, use the [exact binomial limit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html),
and multiply by the native rate. This bounds the weighted yield because each
hard event's selected response weight is between zero and one; the complete
response weights are checked to sum to one per saved source. With zero
selected sources the bound is `sigma * (1 - 0.05**(1/N_attempted))`.
It covers finite sampling with a fixed native rate, not theory or detector
uncertainties. For the saved
classifier, first package its model and a reference of the original scores:

```bash
python3 hgammagamma/frozen_gammagamma_classifier.py \
  --report-dir results/ho_bg_1m_01/xgboost \
  --input-dir results/ho_bg_1m_01/mass_diagnostics/inputs \
  --output-dir results/frozen_diphoton_classifier
```

Pass `--frozen-classifier results/frozen_diphoton_classifier` to the evaluator
to obtain a separate classifier contribution. It replays the original scores
and test yields before using the saved feature order and threshold. New loop
events use full-exposure weights and are neither split nor used for training.
The classifier requires two photon candidates and applies its saved score
threshold; it does not impose the rectangular mass window as an extra cut.
Assess generation-cut, scale, PDF and shower uncertainties before treating a
pilot rate as a physics prediction.

If the broad sample is sparse, a concrete follow-up is three independent
Born-mass samples using native MG5 `mmaa`/`mmaamax` values `(0,100)`,
`(100,150)` and `(150,-1)` GeV, where `-1` disables the upper cut. Allocate
most events to the middle interval and use a distinct seed for each. Sum
`sigma_i * selected_weight_i / all_attempted_weight_i` across the three
disjoint intervals, keeping the same photon cuts, PDFs and SSC response.
The focused sample must not simply be added to the inclusive sample.
The dedicated runner accepts `--gen-mgg-min 100 --gen-mgg-max 150` and
checks the Born mass of every LHE event. A new isolated campaign can copy a
completed export using `--reuse-loop-export /path/to/completed/mg5_process`;
the copy has independent writable files, and its origin is recorded. This
does not reuse events from the old run. Use a new tag and seed.
Alternatively, use the inclusive sample outside the central interval and the
independent focused sample inside it. For each completed campaign, run
[`born_mass_mapping.py`](../born_mass_mapping.py) with `--manifest`,
`--response-root` and a fresh `--output-dir`. It matches the saved hard-photon
four-vectors to the original LHE sequence, checks any shower losses and binds
each detector-response row to its Born mass. Ambiguous matches or mass-window
boundaries stop the audit.

After both mappings and frozen evaluations pass, combine them with:

```bash
python3 hgammagamma/combine_loop_gammagamma.py \
  --broad-report results/loop_broad/frozen_selections.json \
  --focused-report results/loop_focused/frozen_selections.json \
  --broad-mapping results/loop_broad_born_map/summary.json \
  --focused-mapping results/loop_focused_born_map/summary.json \
  --frozen-classifier results/frozen_diphoton_classifier \
  --output results/loop_combined.json
```

The combination checks independent seeds, matching physics cards and response
sources, and reproduces both unpartitioned reports before applying the Born
partition. Each component keeps its full attempted-source denominator. Source
MC variances add across the independent samples. For each selection separately,
the sum of two one-sided 97.5% component limits gives at least 95% coverage;
this does not claim joint coverage across the cut and classifier results.

## Signed normalization

The signal input must contain one undecayed status-one Higgs per event with
compatible beams, mass, PDF and matching provenance. Herwig forces the
diphoton decay with unit proposal branching fraction; the response applies
physical `BR(H -> gamma gamma)=0.00227` once. There is no LO K-factor.
Backgrounds have unit additional weight scale.

HJMiNNLO may write `XSECUP=-1`, missing PDF IDs and `IDWTUP=-4`. The checker
then obtains the native rate from the signed mean `XWGTUP` and PDF IDs from
the embedded card. MG5 supplies the integrated `<init>` rate. The native
signal rate is diagnostic: its physical production rate comes from the
validated ihixs record.

For each sample,

```text
sigma_selected = sigma_production * weight_scale
                 * sum_selected_signed_weight / normalization_sum_weight
N_selected = 1000 * luminosity_fb * sigma_selected_pb
```

For new HO showers, the denominator includes **all attempted source records**
in the same Herwig/BR weight convention as the analysis. Discarded shower
events contribute zero response; they are not removed from the normalization
population. The full LHE is scanned for its weight maximum, while only consumed
records enter the denominator. Legacy campaigns retain their recorded
saved-analysis denominator and exact-count checks.

Reports validate the original LHE, cards, logs, ROOT files and sidecars. Keep
those artifacts at their recorded paths. Negative histogram bins are possible,
especially in small samples; absolute-weight sums cannot replace physical
signed rates. Numerical integration error, scale and PDF uncertainty are
separate. None includes detector, acceptance or branching-fraction uncertainty.

## Existing inputs and individual stages

The common `run` accepts `--signal-lhe PATH` and
`--signal-normalization PATH` for compatible inputs you already calculated.
Without these options it calculates both. Invalid supplied inputs are rejected
before expensive upstream work.

The direct [`run_gammagamma_ho_campaign.py`](../run_gammagamma_ho_campaign.py)
handles `prepare`, `build`, `generate`, `shower`, `analyze`, `normalize`, `all`;
it requires the LHE and normalization record. Its default stage only prepares
cards. Use its `--help` for controlled variations and preserve the original
arguments for `--resume`.

A completed analysis with an older signal-rate sidecar can be updated without
regenerating events:

```bash
python3 hgammagamma/run_gammagamma_ho_campaign.py \
  --stage normalize --run-samples signal_gg_h_aa \
  --run-tag YOUR_EXISTING_TAG \
  --output-dir hgammagamma/HOAnalysis/runs/YOUR_EXISTING_TAG \
  --signal-normalization inputs/pipeline/normalization/ggf-ssc40-n3lo.json
```

This validates saved inputs and updates only signal normalization provenance.
Then remake the report. See [cuts/ML](../../docs/ANALYSIS.md) and
[validation](../../VALIDATION.md) before interpreting results.
