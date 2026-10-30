# Four-lepton campaigns: 13 TeV validation and 40 TeV SSC

The workflows are prepared for you to run. Synthetic checks and an isolated
VBF build test check the implementation; full production and numerical LHC
validation remain separate steps. Start with the pilots below before choosing
larger samples. The [VBF investigation](validation/VBF_BACKEND.md) explains the
switch from MG5 VBF to POWHEG `VBF_H` and records what was tested.

## Minimal student preset

Use `--preset minimal` for a smaller complete campaign, designed around a
40-core allocation and a 1–2 day **target per collider energy**. The runtime
target needs a full timing run; it is not a measured completion time or a
deadline enforced by the runner. No production run is launched by `prepare`.

All 27 samples are included. The default event counts below apply separately
to each of `4e`, `4mu` and `2e2mu`:

| Process | Events per decay channel | Production treatment |
| --- | ---: | --- |
| ggF Higgs | 20,000 | HJMiNNLO + ihixs N3LO normalization |
| VBF Higgs | 10,000 | POWHEG VBF_H NLO+PS |
| W+H, W−H, qqZH, ttH (each) | 5,000 | MC@NLO |
| ggZH | 5,000 | Loop-induced LO |
| qq → 4ℓ | 20,000 | MC@NLO |
| gg → 4ℓ | 2,000 | Loop-induced LO |

This totals **231,000 hard events**, with disjoint Higgs-event slices for the
three decay channels. Zbb, ttZ and ttbar remain separate reducible-background
control studies. The minimal preset does not change PDFs, matching, signal
branching fractions, detector resolutions, analysis observables or the full
109-point ihixs calculation. An existing valid normalization record is reused.

Install/check the loop libraries and optional XGBoost dependencies first:

```bash
python3 run_pipeline.py setup --workflow fourlepton-ho --preset minimal --with-ml
python3 run_pipeline.py check --workflow fourlepton-ho --preset minimal --with-ml
```

Setup preserves a configured COLLIER installation, otherwise looks for one
locally and in the shared runtime. If none is available, minimal setup uses
MG5's COLLIER installer inside the private installation. To choose an existing
build, add `--collier-dir /path/to/collier` to `setup`; that directory must contain
`libcollier.a` and `modules/collier.mod` (or `../include/collier.mod`). Ninja/OneLOop
is discovered from the shared Herwig installation or selected with `--ninja-dir`.
The libraries must use compatible Fortran compilers. Library hashes, generated
cards and the resolved event budget are recorded with the campaign.

For a fresh installation, first run a small execution check:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal_check --jobs 40 --event-multiplier 0.0015
```

That check generates at most 30 events per sample and cannot establish physics
precision. An existing completed pilot at the same energy also satisfies the
managed pilot prerequisite. Then prepare and launch the nominal minimal set:

```bash
python3 run_pipeline.py prepare --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal --jobs 40
python3 run_pipeline.py run --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal --jobs 40
```

For 40 TeV, use `--energy-tev 40` and a different tag, such as
`fourl40_minimal`; that energy requires its own pilot. The `perfect` and `ssc`
analyses share one shower output per sample, with separate detector analyses.
Hard production uses up to 40 cores, one production group at a time. Afterwards,
up to eight sample workers each receive five MadSpin cores; their showers and
analyses run independently. With other `--jobs` values, the worker allocation
is reduced automatically. Logs are in `runs/TAG/workers/PROFILE/SAMPLE.log`;
only the coordinator updates the main manifest. On a worker failure, other
active workers finish, queued work stops, and completed outputs are preserved.

Increase **all** counts uniformly with a multiplier and a new tag:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal_x2 --jobs 40 --event-multiplier 2
```

Multiplier 2 means 462,000 events; multiplier 10 means 2,310,000. Fractional
counts round upwards per sample. `--nevents` and per-sample overrides cannot
be combined with this preset. Repeat the original command to resume; changing
the multiplier requires a new tag. New generated campaigns use tag-dependent
seeds. They are complete replacement datasets, not automatic additions to old
campaigns. Do not concatenate campaigns sharing an external LHE input. Integration
grids are shared by the decay channels within a production group, but are not
automatically transferred between campaign tags.

The runtime settings use external loop libraries, COLLIER first for gg → 4ℓ
(`7|6|1`) and Ninja first elsewhere (`6|7|1`). Loop-induced generation uses
sampled helicities (`nhel=1`) and explicit survey splitting across the assigned
cores. Loop survey accuracy
is 5%; MC@NLO uses `req_acc=0.03`, independently of the event multiplier.
Numerical-stability and pole-cancellation checks are unchanged. These integration
targets are not uncertainties on final selected yields.

Continuum generation uses **90 ≤ m4ℓ ≤ 180 GeV**, summing the four hard charged
leptons and excluding QCD radiation. The usual off-shell amplitudes and
same-flavour interference are included. This window needs a migration check
against an independent campaign with `--wide-continuum`, particularly for the
SSC response and QED radiation. That option removes only the m4ℓ generation
window; the loose lepton generation cuts still apply. No numerical closure or
LHC validation is claimed by selecting the preset.

Run the usual `cuts`, `xgboost` and 13 TeV `validate` actions on the completed tag.
The saved preset selects default background effective-count requirements of
400 after cuts (about 5% MC precision) and 100 for the XGBoost test selection
(about 10%). These use signed weights. The analysis suppresses a sensitivity
claim below those thresholds; increase statistics if needed. An explicit
`--config` uses that analysis card's own thresholds. Omitting `--preset minimal`
for a new campaign keeps the original uniform-count workflow.

## Prepare and launch

Run these commands from the repository root. `prepare` writes review cards
under `hfourlepton/plans/TAG/` and invokes no simulation programs. The managed
`run` action generates its executable cards from the same configuration; editing
a preview card does not change that action.

```bash
python3 run_pipeline.py prepare --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_vbf_pilot --nevents 30
python3 run_pipeline.py prepare --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_vbf_pilot --nevents 30
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
  --energy-tev 13 --run-tag fourl_lhc13_vbf_pilot --nevents 30 --jobs 4
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl_ssc40_vbf_pilot --nevents 30 --jobs 4
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

If an older checkout failed at the first WH decay with
`InvalidParamCard: parameter yukawa: 5 is not at zero`, revision `31b4a30`
corrected that failure.
The old banner validator mistook the `loop_sm-no_b_mass` card for incomplete
metadata and inserted massive-bottom parameters into a massless-bottom model.
The corrected validator accepts the model's zero bottom mass and omitted
fixed-zero Yukawa entry. When replacement metadata is needed, it uses the
declared model's restriction card. This also covers W−H, qqZH, ggZH and ttH.

For the affected 13 TeV pilot, with only the six ggF/VBF samples complete and
the first WH decay failed, run the bounded recovery command using the original
options (including any input overrides) **from that audited revision**.
The recovery deliberately rejects later source versions, including the minimal
preset implementation. A completed pilot does not need this repair. Do not
change its hash allowlist to bypass the version check.

```bash
python3 run_pipeline.py repair-madspin --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_vbf_pilot --nevents 30 --jobs 4 \
  --signal-lhe hfourlepton/runs/fourl_lhc13_pilot/powheg/powheg-hjminnlo-merged.lhe \
  --signal-normalization inputs/pipeline/normalization/ggf-lhc13-n3lo.json
```

This command runs no generators. It checks the specific software transition,
original configuration, completed output hashes, WH production and decay slice.
It backs up the manifest and failed private MadSpin directory under
`recovery/madspin-no-b-mass-v1`, records the migration and makes the original
`run` command resumable. It preserves completed samples and WH hard events.
Other failure states or changed inputs are rejected. Do not use `--force` or
manually remove provenance checks to work around a rejected repair.

For a campaign that failed with the old MG5 VBF definition, first stop its
parent with Ctrl-C in the original terminal/screen session. Keep that directory
and its logs. Run `setup` again to install the pinned `VBF_H` executable, then
start a new tag. Completed ggF production and ihixs normalization can be reused:

```bash
python3 run_pipeline.py setup --workflow fourlepton-ho
python3 run_pipeline.py prepare --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_vbf_pilot --nevents 30 --jobs 4
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl_lhc13_vbf_pilot --nevents 30 --jobs 4 \
  --signal-lhe hfourlepton/runs/fourl_lhc13_pilot/powheg/powheg-hjminnlo-merged.lhe \
  --signal-normalization inputs/pipeline/normalization/ggf-lhc13-n3lo.json
```

These paths refer to the previous 13 TeV pilot. The new tag rebuilds the
four-lepton decay/shower/analysis products without repeating ggF integration or
ihixs. Do not erase manifests or edit their fingerprints to force a resume.
Repeat the complete command, including input overrides, to resume the new tag.
An existing VBF LHE made by `run_powheg_vbf.py` can also be supplied using
`--vbf-lhe PATH`; its process, PDF, pole inputs, scale weights and beam energy
are checked. Use energy-matched inputs at 40 TeV.

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
| VBF | POWHEG VBF_H, Herwig; conventional VBF approximation | Native NLO |
| W+H, W−H | MG5 MC@NLO, Herwig | Native NLO, both charges |
| qqZH | MG5 MC@NLO, Herwig | Native NLO |
| ggZH | MG5 loop-induced LO, Herwig | Native LO, separate from qqZH |
| ttH | MG5 MC@NLO, Herwig | Native NLO |
| qq→4ℓ | Full neutral-current matrix elements, Higgs excluded; MC@NLO | Native NLO |
| gg→4ℓ | Continuum loops, Higgs excluded; Herwig | Native loop-induced LO |

VBF uses the dedicated [POWHEG implementation](https://arxiv.org/abs/0911.5299),
with t-channel weak-boson exchange and without s-channel VH or identical-fermion
interference. VH remains a separate sample. This is a VBF approximation, not
the full electroweak H+2-jet calculation; its neglected terms are not assigned
a numerical uncertainty here, especially without VBF jet cuts.

The pinned `VBF_H` source uses five light flavours and fixed
μR=μF=mW=80.398 GeV. Its electroweak inputs are the source defaults:
mZ=91.1876 GeV, ΓZ=2.4952 GeV, ΓW=2.141 GeV,
αEM=1/128.930 and sin²θW=0.23102. The hard PDF and αs come from LHAPDF
335900. Seven scale weights vary μR and μF by 0.5, 1 and 2, excluding the
opposite extreme pairs. The cards impose no Born jet cuts or suppression.
No ihixs factor is applied to VBF.

For the remaining MC@NLO processes, the runner exports, compiles and runs the
MG5 numerical tests before allowing integration. Each virtual subprocess must
complete its requested pole-check points; zero, partial, missing or aborted
checks fail. MG5's numerical tolerances remain unchanged. An explicit fatal
message, nonzero exit or user interrupt stops the process group started by
that stage and records the failure. Quiet but active work is allowed.
Logs are under `production/GROUP/logs/{export,compile,launch}.log`, with
`nlo_runtime.json` recording the command, host, process ID and status.
`status` flags explicit failure text even if an older runner left its saved
state as `running`. See the [recovery notes](validation/VBF_BACKEND.md).

The ggF calculation follows the diphoton convention, including the full ihixs
scale/PDF scan and its numerical checks. Its hard PDF is NNPDF4.0 NNLO with QED;
the target rate uses the diphoton N3LO PDF convention. Other samples use
`NNPDF40_nlo_as_01180_qed` (LHAPDF 335900), with MC@NLO scales set to
MG5's dynamical choice 3, HT/2. VBF uses the fixed scale above. Loop samples are
explicitly LO despite their NLO PDFs. No common correction factor is assigned
to gg→4ℓ or non-ggF signals.
Electroweak and finite-quark-mass corrections to the ihixs rate are disabled,
as in the existing pure-HEFT diphoton calculation.

All HO Higgs samples use the pole mass mH=125 GeV. HJMiNNLO samples the
event mass within its recorded Breit-Wigner window: the default hwidth is
0.004152 GeV and bwcutoff is 1000, giving 120.848–129.152 GeV. The
`higgsfixedwidth` switch fixes the width prescription, not the event mass.
Validation checks the pole inputs against ihixs and each generated Higgs mass
against those bounds. VBF_H uses the same pole mass and width with a fixed-width
Breit-Wigner over 124.87544–125.12456 GeV (30 widths). Its source revision and
window are checked separately. MG5 stable-H inputs require an on-shell Higgs.
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
The explicit entries use spaces around each numeric value because ThePEG 2.3
parses whitespace-separated tokens. XML-valid entries without these spaces
can silently lose their variations. With `thepeg-config` on `PATH`,
`python3 -m unittest discover -s tests -p test_fourlepton_weights.py` also checks
the converted file with the installed reader, without generating or showering
events.

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
python3 -m unittest discover -s tests -p test_fourlepton_powheg_vbf.py
python3 -m unittest discover -s tests -p test_fourlepton_mg5_guard.py
python3 -m unittest discover -s tests -p test_fourlepton_campaign.py
python3 -m unittest discover -s tests -p test_fourlepton_report.py
python3 -m unittest discover -s tests -p test_ihixs_normalization.py
python3 -m unittest discover -s tests -p test_student_pipeline.py
```
