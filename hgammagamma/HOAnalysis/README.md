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
photon fragmentation is absent. Separate loop-induced `gg -> gamma gamma`,
dijet double fakes and signal–background interference are omitted. Gamma+jet
uses one hard fake candidate with the [documented flavour model](../SSC_DETECTOR_RESPONSE.md),
whose ambiguity grows with extra NLO partons.

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
