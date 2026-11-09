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

To increase just the three qq continuum backgrounds, use
`--qq-event-multiplier`. A factor of 10 gives 200,000 events per qq channel
and keeps the other minimal counts in the table above: **771,000 hard events**
in total. After the 40 TeV pilot has passed, launch this campaign on Timur with:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl40_minimal_qq200k \
  --preset minimal --qq-event-multiplier 10 --jobs 80
```

This produces both perfect and SSC detector analyses. The existing 40 TeV ihixs
normalization is reused; the new campaign generates its own hard events.
The global `--event-multiplier` also scales the boosted qq counts: adding
`--event-multiplier 2` doubles the whole 771,000-event set. The two factors
are multiplied before event counts are rounded upwards. Their values are
recorded in the plan and campaign identity; repeat the same options to resume
and use a new tag to change either factor. The default qq factor is 1.

After production and the strict report finish, run the two analyses:

```bash
python3 run_pipeline.py cuts --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl40_minimal_qq200k --luminosity-fb 100
python3 run_pipeline.py xgboost --workflow fourlepton-ho \
  --energy-tev 40 --run-tag fourl40_minimal_qq200k --luminosity-fb 100
```

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

### Recover a run stopped after the W+H export

The first minimal-preset release checked NLO exports using the loop-induced
build layout. An error mentioning a missing
`SubProcesses/MadLoop_makefile_definitions` immediately after W+H export came
from this wrapper check. MG5 NLO exports instead put the link settings in
`Source/make_opts` and `SubProcesses/makefile_loop`. The corrected check audits
both NLO files, including their library and include settings.

Keep the failed directory. Since its software fingerprint differs after the
fix, use a new tag and supply the completed ggF and VBF files explicitly:

```bash
python3 run_pipeline.py run --workflow fourlepton-ho --preset minimal \
  --energy-tev 13 --run-tag fourl13_minimal_nlofix --jobs 40 \
  --signal-lhe hfourlepton/runs/fourl13_minimal/powheg/powheg-hjminnlo-merged.lhe \
  --vbf-lhe hfourlepton/runs/fourl13_minimal/powheg-vbf/powheg-vbf-merged.lhe
```

This reuses both completed Higgs productions and the cached ihixs normalization;
it exports W+H afresh and then continues the remaining steps. Use the new tag
for analysis, and do not combine it with outputs from the failed campaign.
For another original tag or event multiplier, adjust these options to match
that campaign. Do not use `--force` or edit provenance records to bypass the
saved configuration checks.

### Recover the MG5 loop-combination failure

An `incomplete production export` message on restart can follow a failed
integration. Inspect the original log before changing tags or event counts.
The observed `bkg_gg_4mu` failure in MG5 3.5.15 was an extreme-weight warning
followed by `gen_ximprove_share.combine_iteration` failing to read
`G*/results.dat`. Its unstable-job recovery deletes those files before a
recursive call tries to read them again. The corresponding
[upstream code](https://github.com/mg5amcnlo/mg5amcnlo/blob/3.x/madgraph/madevent/gen_ximprove.py)
also needs consistent job indices and event-file selection after an exclusion.

Managed minimal runs now install a checksum-checked fix in the private MG5
installation. The original Python source and a patch receipt are saved there;
`runs/TAG/mg5_loop_recovery.json` records the fix and any exports that already
existed. Completed exports and LHE files are unchanged. The fix postpones
cleanup until combination succeeds, tracks repeated exclusions by original
job number, and excludes the same jobs from the integral and the event file.
It preserves MG5's existing numerical thresholds, library order and outlier
criterion. Software tests exercise the actual MG5 methods with synthetic
results; a successful production rerun is a separate validation step.

For this specific failure, audit the existing campaign first:

```bash
python3 hfourlepton/repair_loop_production.py \
  --run-tag fourl13_minimal_nlofix --sample bkg_gg_4mu
```

The audit checks the original configuration, the failure evidence, and all
preceding completed groups, including their hashes and event counts. To archive
only the incomplete group and then continue using the original run command:

```bash
python3 hfourlepton/repair_loop_production.py \
  --run-tag fourl13_minimal_nlofix --sample bkg_gg_4mu --apply
bash hfourlepton/plans/fourl13_minimal_nlofix/run.sh
```

`--apply` runs no generators. It moves the failed group, including partial
grids and logs, into `runs/TAG/recovery/higgsssc-loop-recovery-v1/SAMPLE/` and
saves its audit and original metadata. The run then reuses the eleven completed
groups, regenerates `gg → 4μ`, and continues with `gg → 2e2μ` and postprocessing.
Some intermediate results were deleted before the crash, so the interrupted
refinement cannot be resumed reliably. This repair refuses final LHE files,
unrecognized failures, changed provenance, or campaigns that have begun
postprocessing. Do not use `--force` or remove the completed production groups.

### Recover a MadSpin logging-startup failure

The first bounded MadSpin launcher omitted MG5's formatter-registration import.
If a worker log ends with `No module named 'logging.ColorFormatter'`, it failed
before importing its event card or generating any decays. The corrected launcher
uses the initialization from the stock MadSpin frontend. Managed minimal runs
now check that startup and the private core allocation before hard production:

```bash
inputs/pipeline/venv/bin/python hfourlepton/madspin_worker.py \
  --mg5-dir inputs/pipeline/MG5_aMC_v3_5_15 --cores 5 --check
```

This check generates no events and does not save MG5 configuration changes.
The migration below requires every hard-production group to have its full
event count. The `fourl13_minimal_nlofix` campaign also contains an undersized
`bkg_gg_4mu` result (19 of 2,000 events), which must be regenerated first using
the next subsection. Then audit and apply the specific metadata migration:

```bash
python3 hfourlepton/repair_madspin_startup.py --run-tag fourl13_minimal_nlofix
python3 hfourlepton/repair_madspin_startup.py --run-tag fourl13_minimal_nlofix --apply
bash hfourlepton/plans/fourl13_minimal_nlofix/run.sh
```

The migration accepts only the known old/new launcher hashes and this logging
failure. It verifies every completed hard-production group, requires that no
decays or showers have started, and checks the prepared decay inputs. It archives
the original manifest and worker metadata/logs under
`runs/TAG/recovery/madspin-logging-startup-v1/`, then updates the campaign
fingerprint with an explicit migration record. All hard events, production
markers, decay cards and staged LHE inputs are unchanged. `--apply` runs no
simulation; the original run command resumes the same tag at postprocessing.

### Large SSC reports and analysis memory

SSC efficiency enumeration can produce tens of millions of detector hypotheses
from 200,000 source events. Reporting reads ROOT files in 32 MB chunks using
Uproot, combines rejected hypotheses with the same source, channel and cut mask,
and keeps each selected hypothesis separately. This preserves the cutflows,
source-grouped statistical uncertainties, probability closure and optional-weight
audits. Checks against the analyzer summary use the original tree entry count.
The log shows each sample and progress through its hypothesis tree.

Cuts and XGBoost use the same reader. A zero-weight row for every source event
keeps the full source population in train/validation/test fractions; only selected
hypotheses need additional analysis rows. Generated ROOT files are unchanged.

`setup --workflow fourlepton-ho` (or `fourlepton-lo`) installs the reader in the
private environment without upgrading the existing core packages. After an
interrupted report, repeat the original `run` command: the pipeline verifies the
completed analysis inventory and resumes at reporting, then the usual `cuts` and
`xgboost` commands can run. No event generation or detector reconstruction is
needed. Preserve the interrupted log before restarting.

### Recover an undersized loop-induced sample

MG5 can exit successfully after reporting `fail to reach target`, leaving an
unweighted LHE with far fewer events than requested. Managed minimal runs now
use `checked_campaign.py` to verify native hard-production counts and provenance
before starting any decay/shower workers. Each passing group gets a separate
receipt in `runs/TAG/production_validation/`. Signal counts are checked on the
full production group, before its disjoint decay-channel slices are compared.

For the affected gg → 4μ sample, these commands first audit the saved campaign,
then explicitly regenerate **only that group**:

```bash
python3 hfourlepton/retry_short_production.py \
  --run-tag fourl13_minimal_nlofix --sample bkg_gg_4mu
python3 hfourlepton/retry_short_production.py \
  --run-tag fourl13_minimal_nlofix --sample bkg_gg_4mu --run
```

Use the affected campaign's tag, for example `fourl40_minimal_qq200k` for
the 40 TeV campaign. The helper restores its saved energy, job allocation,
global event multiplier and qq-background multiplier, and verifies the
complete configuration before changing any files.

Without `--run`, no campaign files change and no generator starts. With `--run`,
the original LHE, grids, export and metadata are archived under
`runs/TAG/recovery/loop-short-sample-v1/bkg_gg_4mu/`. The replacement uses the
original managed environment and job allocation, with `sde_strategy=2`
and `hard_survey=1`. These are documented MG5 controls for changing phase-space
sampling and strengthening the initial survey; see Table 8 of
[Speeding up MadGraph5_aMC@NLO](https://pmc.ncbi.nlm.nih.gov/articles/PMC8136271/).
The process, PDF, cuts, seed, requested count and reduction-library order stay
the same. The exact numerical changes and replacement artifact are recorded in
`production/bkg_gg_4mu/numerical_retry.json` and linked from its production marker.

This retry is a numerical remedy whose speed and stability must be checked
for the affected campaign. It must produce the requested event count
(2,000 at the default minimal multiplier) and record the requested settings
in its LHE banner. Inspect numerical warnings as well as the final count.
If it fails again, the helper preserves
the second attempt and refuses another automatic replacement. Do not use
`--force` to bypass that failure. After the replacement passes, repeat the
original run command to reuse completed groups and finish any queued groups.
The historical 13 TeV campaign also needs the MadSpin startup migration above;
that migration does not apply to the 40 TeV campaign.

The first retry helper required `hard_survey` to be present in a freshly
exported run card, although MG5 3.5.15 defines it as a hidden parameter with
default zero. The corrected helper inserts an explicit `hard_survey=1` when
that entry is absent. It still rejects duplicate or malformed settings and
requires both retry settings to be explicit in the output LHE banner.

If the retry stopped with `expected exactly one explicit hard_survey setting`,
repeat the same audit and `--run` commands above after updating the checkout.
The helper recognizes this specific failure only before compilation or
integration, checks the original archive and the other completed groups, and
archives the unfinished export under
`recovery/loop-short-sample-v1/bkg_gg_4mu/hidden-survey-card-startup-v1/` before
making a fresh export. The original 19-event output remains in its first
archive. This does not bypass a failed numerical attempt. On Timur, the existing
`hfourlepton/plans/fourl13_minimal_nlofix/resume_madspin.sh` launcher performs
these steps and continues with the MadSpin migration after successful generation.

## Prepare and launch

### Finite-LHE shower replay recovery

ThePEG 2.3 can silently reopen an exhausted LHE file on the final requested
event even when `AllowedToReOpen No` is set. In the November 2 minimal run,
one discarded event led to one extra input attempt and a repeated first hard
event at the end of each affected shower. The completion validator rejected
these outputs with `inconsistent source/generated/saved shower counts`.

New four-lepton HO showers request `N_LHE + 2` execution slots and require the
first end-of-file exception with exactly `N_LHE` source attempts. The extra
slots are an end-of-file guard, not additional generated events. One extra slot
would still hit the bug when no event is discarded. The execution plan is saved
as `herwig/shower_execution.json` and hashed into the completion record. All
events in the finite source are processed, including a native producer's small
surplus over the requested budget. Discarded shower events have zero response
and their signed source weights still enter the normalization denominator.
Existing valid completion records remain usable.

For the affected managed minimal campaign, audit and then apply the narrow
recovery before repeating its original launch command:

```bash
python3 hfourlepton/repair_shower_replay.py --run-tag fourl13_minimal_nlofix
python3 hfourlepton/repair_shower_replay.py --run-tag fourl13_minimal_nlofix --apply
bash hfourlepton/plans/fourl13_minimal_nlofix/run.sh
```

The first two commands start no simulations. The repair validates hard
production, completed analyses, decay provenance and the specific failed
showers. It archives the invalid showers, worker logs and original metadata
under `runs/TAG/recovery/finite-lhe-shower-v1/`, and migrates only the known
shower-launcher checksum change. Hard events, completed decays and successful
analyses are unchanged. The final command repeats the affected showers and
continues the unfinished samples. Do not use `--force` or delete the run.

### Completed-manifest checksum after a report retry

Older managed launchers re-entered the campaign when only reporting remained.
The minimal workers reused their outputs but overwrote timing fields in the
completed manifest, so the subsequent frozen-artifact check rejected that
manifest. The corrected launcher validates the frozen analysis inventory and
resumes directly at reporting. It leaves completed manifests and worker
metadata unchanged on repeated runs.

If an earlier retry has already stopped with `completed analysis product
changed: .../manifest_perfect.json`, use the narrow recovery tool:

```bash
python3 hfourlepton/repair_completed_manifest.py --run-tag fourl13_minimal_nlofix
python3 hfourlepton/repair_completed_manifest.py --run-tag fourl13_minimal_nlofix --apply
bash hfourlepton/plans/fourl13_minimal_nlofix/run.sh
```

The audit requires all 54 ROOT/summary hashes to match the original frozen
inventory, all workers to have resumed already-complete samples, unchanged
physics fields in their saved job manifests, and valid configuration, hard
production, source and shower provenance. The original manifest bytes are not
reconstructed from an old checksum. Instead, an explicit re-audit authorizes
updating only its managed checksum. The old state, current manifest and worker
evidence are copied to `runs/TAG/recovery/completed-manifest-resume-v1/`.
No LHE, ROOT, analyzer summary, campaign manifest or report is rewritten by the
repair. The last command rebuilds the report and marks the pipeline complete.

### Campaign commands

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

This is an approximate particle-level benchmark. Legacy HwSim outputs lack
prompt/hadron/Higgs ancestry; CMS v3 requires the optional ancestry plugin
described below. Associated-production combinatorics can still differ from
the experiment. The
`objects` tree contains soft stable particles as well: `PTCutParticles=0.4`
controls reconstructed jet/MET inputs, not the saved stable-particle list.
CMS detector
ID, impact-parameter and vertex selection are not emulated. An
[ATLAS four-lepton measurement](https://arxiv.org/abs/2004.03969) also uses
different kinematic thresholds, a matrix-element choice with extra leptons,
and lepton–jet overlap rules. No separate numerical ATLAS comparison or
experimental likelihood is implemented. Agreement in this CMS benchmark
would support the basic selection and acceptance implementation; it would
not validate the historical SSC resolution or establish discovery reach.

### Reprocess the CMS selection using saved showers

The `cms_fiducial_v2` selection corrects two differences found in the first
13 TeV minimal campaign. Photon dressing uses ΔR < 0.3, following the
fiducial definition in [CMS HIG-16-041, section 10.2](https://arxiv.org/pdf/1706.09936),
which [CMS HIG-19-001](https://arxiv.org/pdf/2103.04956) also uses. For four
same-flavour leptons, a candidate is rejected if an alternative pairing has
a Z closer to its nominal mass and a companion pair below 12 GeV
(HIG-19-001, section 4). Other valid quadruplets in the event are still
considered. GEM/SSC selection and detector response are unchanged.

From the repository root on Timur, run:

```bash
python3 run_pipeline.py reprocess --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl13_minimal_nlofix \
  --revision cms_v2 --jobs 4
```

This verifies the completed shower records and original analysis checksums,
builds a private analyzer, and writes a separate reconstruction, strict report,
acceptance comparison, cut analysis and XGBoost analysis. It does not generate
events, run MadSpin, shower events, or replace the original campaign products.
Seeds, signed source-weight normalization, optional weights and hard-event
group identities are preserved. The new manifest and both analyzer metadata
records identify the selection version and its dressing/veto settings.

Results are under
`hfourlepton/runs/fourl13_minimal_nlofix/revisions/cms_v2/index.html`.
Use the same command to resume an interrupted revision; it verifies and skips
completed outputs. Failed attempts are kept in separate directories. A changed
analyzer, input, card or analysis setting requires a new revision name.

```bash
python3 run_pipeline.py status --workflow fourlepton-ho \
  --energy-tev 13 --run-tag fourl13_minimal_nlofix --revision cms_v2
```

An acceptance result of `not_validated` is recorded separately from execution
completion, and the cut/ML analyses still run. Statistical requirements stay
at 400 effective background events for cuts and 100 for the minimal ML card.
Missing sensitivity estimates or an unset classifier threshold still require
more statistics; they are not converted into successful physics checks.

For later custom cut/ML cards, invoke `analyze_fourlepton.py` with the revision's
`manifest_perfect.json`, the chosen card, and a fresh output directory. Reuse
the existing event groups for a consistent train/validation/test split.

The November 2026 read-only ttH diagnostic gave 34.75% acceptance with the
original selection and 35.92% with both corrections, compared with 31.50% in
CMS Table 8. This is still an unresolved acceptance comparison. Correcting
these selection rules does not establish closure. Revisit ancestry and
associated-production combinatorics before interpreting a remaining offset.
The original background variance was dominated by qq→4ℓ. Use the corrected
outputs to revise the event budget for those channels first.

The [3 November CMS validation audit](../docs/FOURLEPTON_VALIDATION_AUDIT.md)
identifies the remaining photon-ancestry and isolation mismatch, quantifies
its effect with read-only variations, and projects the qq-background event
budget. The v2 acceptance comparison remains provisional; counting all
photons in isolation is a diagnostic variation, not a validated correction.

### CMS v3 with HwSim ancestry

The [private HwSim ancestry plugin](hwsim/README.md) records particle history
and enables prompt-lepton selection, dressing by the emitting lepton's lineage,
and isolation with only that lepton's FSR removed. Use
`--ancestry-plugin-dir PATH` when reprocessing a completed campaign or starting
a new 13 TeV campaign. New campaigns using this option build their analyzer
inside the run directory. The existing shared HwSim installation and analyzer
binary are not replaced.

The ancestry smoke test used 200 events each from ggF, qqZH, ttH and qq→4ℓ.
All 800 events reproduced the original stable-particle momenta and weights
exactly. This checks that recording ancestry does not alter the shower; the
full CMS acceptance comparison is a separate test.

For an independent background statistics campaign, the direct runner supports
selecting only the qq channels and scaling their minimal budgets:

```bash
PYTHONNOUSERSITE=1 inputs/pipeline/venv/bin/python hfourlepton/checked_campaign.py \
  --order ho --selection-profile cms --detector-response perfect \
  --run-tag fourl13_qq200k_truth --run-samples bkg_qq \
  --preset minimal --event-multiplier 10 --seed-base 20261103 \
  --nb-core 32 --ebeam 6500 --higgs-mass 125 \
  --mg5-dir inputs/pipeline/MG5_aMC_v3_5_15 \
  --mg5-python "$PWD/inputs/pipeline/venv/bin/python" \
  --no-herwig-env --no-herwig-module \
  --ancestry-plugin-dir inputs/pipeline/hwsim-ancestry-v1
```

Run this inside the Herwig environment. The three qq channels receive 200,000
events each, with the existing 90–180 GeV generation window and NLO matching.
The independent seed avoids repeating the earlier event population. Replace
the old qq samples in a combined analysis manifest; do not add two copies of
their physical cross sections. Thirty-two generation cores leave eight cores
for the ancestry validation showers on a 40-core allocation.

Once both campaigns are complete, assemble them and rerun the strict report,
acceptance check, cuts and XGBoost without changing either input campaign:

```bash
PYTHONNOUSERSITE=1 inputs/pipeline/venv/bin/python hfourlepton/replace_qq_backgrounds.py \
  --base-manifest hfourlepton/runs/YOUR_RUN/revisions/cms_v3/manifest_perfect.json \
  --qq-manifest hfourlepton/runs/fourl13_qq200k_truth/manifest_perfect.json \
  --output-dir hfourlepton/runs/YOUR_RUN/revisions/cms_v3_qq200k
```

This checks the physics settings and source checksums, replaces exactly the
three qq channels, and keeps the cut/ML effective-count requirements at 400/100.

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
