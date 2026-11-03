# Four-lepton CMS validation audit — 3 November 2026

The completed `fourl13_minimal_nlofix/revisions/cms_v2` campaign has 27 samples,
zero strict-report errors or warnings, and an acceptance result of
`not_validated`. This investigation reads its existing showers. It does not
change the selection, normalization, event budgets, or completed products.

## Acceptance diagnosis

The largest identified definition mismatch is the treatment of photon origin
in dressing and isolation. The analyzer dresses every stable electron or muon
with nearby photons, regardless of origin. It then removes photons assigned
to **any** lepton from every lepton's isolation sum. Consequently, photons
from hadron decays can increase a lepton's momentum and disappear from its
isolation activity. Nonprompt leptons also participate in candidate selection.

The published fiducial definition uses dressed leptons and relative isolation
in a cone of 0.3: [CMS HIG-16-041, section 10.2](https://arxiv.org/pdf/1706.09936),
carried forward in [CMS HIG-19-001, section 10.3 and Tables 7–8](https://arxiv.org/pdf/2103.04956).
The public [Run-2 UF analysis implementation at commit 84a1d522](https://github.com/lucien1011/UFHZZAnalysisRun2/blob/84a1d5220218627437a86bf9a2c55c6b7f49796c/UFHZZ4LAna/src/UFHZZ4LAna.cc#L4103)
provides a useful implementation cross-check: it selects leptons using their
boson ancestry, requires a lepton parent for dressing photons, and removes
that lepton's identified FSR photons from isolation. The public code is
supporting evidence for the convention; it is not proof of the exact software
revision used to produce Table 8.

Controlled variations on all 21 saved signal samples give the following
acceptances in percent. The original event weights and denominators are used
in every column.

| Mode | CMS Table 8 | Current v2 | Remove only the candidate's assigned photons from isolation | Count all photons in isolation |
|---|---:|---:|---:|---:|
| ggF | 40.20 | 42.53 | 42.33 | 38.67 |
| VBF | 44.50 | 46.60 | 46.18 | 42.78 |
| WH | 32.90 | 35.64 | 35.29 | 32.29 |
| ZH | 34.00 | 37.91 | 37.56 | 34.35 |
| ttH | 31.50 | 35.92 | 35.09 | 31.32 |

The last column is a sensitivity test, **not a proposed CMS correction**:
it also counts genuine FSR photons in isolation while leaving the dressed
momenta unchanged. Its improved agreement does not establish validation.
These variations show that photon treatment is large enough to account for
the observed offsets; they do not measure the exact correction from proper
ancestry information.

Two other checks narrow the diagnosis:

- Choosing Z1 first by proximity to the Z mass, then Z2 by the largest scalar
  lepton-pT sum, changes ZH from 37.9078% to 37.8973% and ttH from 35.9185%
  to 35.8477%. Candidate ordering alone does not explain the offsets.
- The ggZH component supplies 5.62% of the inclusive ZH rate and has 47.92%
  acceptance. The qqZH component alone has 37.31% acceptance, so changing
  the ggZH mixture cannot explain the full ZH discrepancy. An overall
  normalization cancels within each component; ihixs applies only to ggF.

The saved `partons` branch contains the four Higgs-decay hard leptons, but
direct angular matching to stable leptons within delta-R 0.1 fails for about
60.4% of the signed ZH population and 54.4% of ttH. The hard and showered
momenta need not have the same recoil. Those matches are not a valid ancestry
substitute, and their diagnostic rates must not be quoted as truth-level
acceptances. No claim about the precise combinatorial contribution follows
from this matching attempt.

## What a proper correction needs

Save lepton and photon ancestry in HwSim: the prompt boson/Higgs ancestry of
leptons, the emitting lepton lineage of FSR photons, and hadron/tau ancestry
needed to distinguish other sources. Apply the documented photon and lepton
definitions to those records and version the selection and metadata.

This requires another shower pass to produce the additional information;
it cannot be reconstructed exactly from the current stable-particle arrays.
The existing hard-production and MadSpin LHE files can be reused. First use a
small separate ZH/ttH validation sample to check the new records and the
acceptance, then process the full set. The user launches these jobs. The
current 10% comparison tolerance and statistical requirements should stay
unchanged while the definition is corrected.

## Background event budget

With the current v2 definition, the 120–130 GeV cut region at 100 fb^-1 has
253.04 +/- 2.00 signal events and 89.45 +/- 8.78 irreducible-background events
(Monte Carlo errors only). The background effective count is 103.82, below
the cut card's requirement of 400. The qq continuum accounts for 99.943% of
the background variance. The XGBoost validation partition has only 29.58
effective background events before a score cut, against a requirement of 100.

A useful next planning target is **200,000 events per qq channel**, compared
with the present 20,000; this is 600,000 qq events in total. Keeping the gg
statistics fixed and assuming unchanged distributions, this projects to
about 1,033 effective background events for the full cut sample. The saved
validation scores give these approximate projections:

| Illustrative score cut | Current validation effective background | With 10x qq statistics |
|---|---:|---:|
| none | 29.6 | 293.9 |
| 0.4 | 14.1 | 139.8 |
| 0.5 | 13.0 | 128.9 |
| 0.6 | 11.9 | 117.9 |
| 0.7 | 8.9 | 88.2 |

These are variance-scaling projections at fixed scores, not selected working
points or guarantees after retraining. The tighter rows depend on very few
current events. Reassess the budget after correcting the object definition;
200,000 per channel is a starting point for moderate score cuts. Increasing
all signal and loop-induced samples is not needed to address this variance.

About 71% of the generated qq events have hard four-lepton mass below 100 GeV.
A future stratified mass campaign could improve efficiency. Simply cutting
away this region is premature: with v2, one such mixed-flavour event enters
the 120–130 GeV window, and additional migrations occur in the 105–140 GeV
preselection. Preserve coverage and proper cross-section normalization in
each mass slice, and measure migrations with the corrected object definition.
No generation cuts or event counts were changed in this investigation.

## Reproducing the definition probe

The diagnostic compiles against the analyzer source saved in the revision,
checks its hash and the shower-file hashes, and writes into a fresh directory.
It reads all 27 existing samples, uses four workers by default, and does not
run MG5, MadSpin, Herwig, training, or the campaign pipeline.

```bash
PYTHONNOUSERSITE=1 inputs/pipeline/venv/bin/python \
  hfourlepton/validation_tools/probe_saved_showers.py \
  --manifest hfourlepton/runs/fourl13_minimal_nlofix/revisions/cms_v2/manifest_perfect.json \
  --output-dir /tmp/fourl-cms-definition-audit --jobs 4
```

The all-sample diagnostic took approximately two minutes on Timur. Results
include weighted counts and hard-mass migration bins; they are diagnostic
variations, not new nominal analysis outputs. The baseline reproduces the
five saved acceptance numbers. The original revision manifest SHA-256 is
`13e7d1b263fbfeb7cdaeb5ee3150c7d159b463642e5127038ea4b168ad1f36e4`.
