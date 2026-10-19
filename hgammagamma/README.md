# Higgs to two photons

Use the [top-down workflow](../README.md) to install and run the code. This page
records the LO physics choices; [HO](HOAnalysis/README.md) uses the same
[SSC/GEM response](SSC_DETECTOR_RESPONSE.md) with different hard processes and
normalization.

## Samples and conventions

| Sample | MG5 process | Response |
| --- | --- | --- |
| `signal_gg_h_aa` | `g g > h [noborn=QCD]`, MadSpin `h > a a` | genuine photons |
| `bkg_prompt_aa` | `p p > a a` | genuine photons |
| `bkg_gamma_j` | `p p > a j` | one jet promoted to a photon |
| `bkg_dy_ee` | `p p > e+ e-` | two electrons promoted to photons |

The beams are protons at 20 TeV each. The private MG5 installation contains
the repository's `loop_sm_haa` model, which supplies the effective
Higgs–photon coupling for the decay. MadSpin uses `spinmode none` for the
spin-zero Higgs. The LO PDF is `NNPDF40_lo_as_01180` (LHAPDF ID 331900) in
both the hard process and shower. The source process/run cards and campaign
manifest record the mass, scales, cuts and seed used in each run.

The teaching LO signal rate uses `K(ggH)=2` and
`BR(H -> gamma gamma)=0.00227`; background rate factors default to one.
The branching fraction follows the [CERN Yellow Report reference](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR).
These conventions do not make the LO sample a precision prediction. NLO/NNLO
PDFs must not be substituted while keeping the LO setup unchanged.

The default backgrounds omit the separate loop-induced `gg -> gamma gamma`
box contribution, dijet double fakes and signal–background interference.
Gamma+jet flavour assignment and rare fake probabilities carry detector-model
assumptions; they are not measured detector performance.

## Outputs and interpretation

New campaigns live in `hgammagamma/LOAnalysis/runs/TAG/`. Each sample has a
manifest, MG5 LHE input, Herwig/HwSim ROOT output and analysis products:

| File | Contents |
| --- | --- |
| `*.top` | Weighted histograms |
| `*.dat` | Settings, event counters, weights and closure diagnostics |
| `*_var.root` | Per-hypothesis variables and weights for cuts/classifiers |
| `*.evp` | Source indices with a valid diphoton hypothesis |

A source event can produce several exclusive detector hypotheses. Use summed
weights for yields, not ROOT row counts. The physical yield is
`1000 * luminosity_fb * selected_cross_section_pb`. Unit-area plots show shapes
only; display magnification does not increase the signal yield.

LO and HO reports show both views by default: unit-area densities on the left
and stacked cross sections per bin in pb on the right. On narrow screens the
views appear one below the other. Each view has its own PNG, SVG and CSV
downloads; the density files have a `_density` suffix. Cross-section files
keep their existing names. A non-positive signed histogram sum can prevent
unit-area normalization; in that case the report explains why the shape is
unavailable and still shows the physical cross-section plot.

The direct `make_gammagamma_report.py` command accepts
`--plot-mode density` for shapes only and `--plot-mode stacked` (or the
existing `--no-density`) for stacked plots only. Single-view reports use the
original filenames without a suffix. `--normalization` controls the stacked
view; its default is `event_xsec`. `--luminosity-fb` adds physical event yields
to the summary, without changing the plotted cross sections.

The default `ssc` analysis includes smearing, acceptance, isolation,
efficiencies and fake probabilities. The direct runner's legacy `none` mode
uses different cuts and cannot model gamma+jet or dielectron fakes; it is not
an SSC run with smearing simply disabled.

See [cuts and classifiers](../docs/ANALYSIS.md),
[recovery](../docs/TROUBLESHOOTING.md) and
[validation evidence](../VALIDATION.md). The direct
[`run_gammagamma_campaign.py`](run_gammagamma_campaign.py) remains available
for controlled sample, card and detector variations; use a new tag for each
configuration.
