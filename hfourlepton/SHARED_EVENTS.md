# Shared four-lepton events

The student entry point is documented in the
[root README](../README.md#analyze-the-prepared-four-lepton-events-on-timur).
Use `python3 student_events.py attach` from a student-owned checkout on Timur.
All local commands write to the student's workspace; the release is read-only.

## Publish a completed revision

Publishing copies completed data and metadata. It never calls an event
generator, MadSpin, a shower, or the reconstruction executable. On Linux it
uses copy-on-write copies when the filesystem supports them, falling back to
ordinary copies. It never uses hardlinks into the instructor's working files.

The 2026-11-12 release is assembled with:

```bash
python3 hfourlepton/shared_events.py \
  --perfect hfourlepton/runs/fourl40_minimal_qq200k/revisions/ssc_qq1m/manifest_perfect.json \
  --ssc hfourlepton/runs/fourl40_minimal_qq200k/revisions/ssc_qq1m/manifest_ssc.json \
  --reference-dir hfourlepton/runs/fourl40_minimal_qq200k/revisions/ssc_qq1m \
  --extra-ssc hfourlepton/runs/fourl40_minimal_qq200k/revisions/ssc_qq1m_extended_test/manifest_additional_ssc.json \
  --extended-reference hfourlepton/runs/fourl40_minimal_qq200k/revisions/ssc_qq1m_extended_test \
  --destination /home/shared/HiggsSSC/fourlepton/ssc40-20261112
```

Use a **new version directory** for another release. An interrupted publication
can resume with exactly the same inputs and destination. An existing completed
release cannot be overwritten. Ensure the parent directory is readable and
traversable by student accounts, and allow enough space for ordinary copies if
the filesystem cannot make copy-on-write copies.

The publisher verifies the event-file checksums against the completed source
manifests. It relocates operational paths in manifests, summary sidecars,
input lists and shower evidence, then runs the existing completed-manifest
audits against the shared copies. Event files, event counts, normalization,
signed-weight moments, seeds and physics settings remain unchanged.
Original manifests and the original versions of relocated evidence are saved
under `provenance/`. Historical paths inside those originals, ihixs provenance
or ROOT metadata describe how production was performed; they are not runtime
dependencies for the student commands.

`dataset.json` is written only after the audits pass and lists the release's
files, sizes and SHA-256 hashes. Publication finishes by setting files to
mode `0444` and directories to `0555`. Students cannot alter the shared data.
Access does not require changing permissions on the instructor's home.

## Data layers

- `manifest_perfect.json` and `manifest_ssc.json`: the same 27 nominal samples
  with different detector responses. Use one profile for cuts or training;
  the report compares both.
- `nominal/analysis/`: unchanged reconstruction ROOT trees with relocated
  sidecars and weight-name mappings. This is the normal student starting point.
- `nominal/events/`: LHE inputs, HwSim shower ROOT files and shower-completion
  evidence. These allow a separate reconstruction study without generation.
- `nominal/production/`: stable-Higgs inputs documenting the disjoint decay
  slices used to prevent train/test overlap across channels.
- `cards/`: instructor starting cards; `attach` gives each student editable
  copies. The release versions remain unchanged.
- `supplemental/`: independent qq samples for a separate classifier test.
  They are excluded from the nominal manifests.
- `reference/`: the completed instructor reports, model and test-extension
  results. These are historical reference products, not new student results;
  their original command/manifest paths can refer to the instructor's campaign.

The reference products and source records are copied for inspection, not
rewritten to claim that they were produced in the shared directory. A fresh
student analysis records its own manifest, cards, command and outputs.

## Check access

From a clean checkout in another Timur account:

```bash
python3 student_events.py attach
python3 student_events.py setup
python3 student_events.py check --full
python3 student_events.py cuts
```

The first three commands check access, Python dependencies and all checksums.
The last checks the complete analysis path with new, private output. Running
`report` additionally verifies the paired detector responses and creates the
observable plots; `xgboost` trains a new classifier from the nominal dataset.
