# Replication and production records

## Refit from a clone

From the repository root, `python student_events.py discovery` reproduces the discovery calculation in a few seconds to a minute on an analysis workstation. Install the three packages in [requirements.txt](requirements.txt) first. The same command runs on Timur. The commands below verify the distributed inputs and write a new results directory.

```bash
(cd discovery/data && shasum -a 256 -c SHA256SUMS)
python student_events.py discovery --output results/discovery-check
```

On Linux, `sha256sum -c SHA256SUMS` also works. The published inputs and fit are dated 23 November 2026. The archived fit SHA-256 is `0e118f60f9ef74abde08db9a07792ab93cf4a1c532160d53b3c73a6e484016e4`. [Published numbers](published/results.json) · [Input checksums](data/SHA256SUMS)

| File | Content and use |
| --- | --- |
| [spectra.tar.gz](data/spectra.tar.gz) | Seven reference spectra (six processes with the loop in two Born regions), 51 earlier packets and 696 final packets. Each packet includes its independent source identifier, bin rates in fb, full 55 × 55 covariance in fb² and event exposure. An internal SHA-256 manifest covers every member. |
| `records/` inside the spectrum archive | Original final fit and reach, extended four-lepton test, zero-count heavy-flavour input record and the iHixs normalization. Some historical records identify generator paths in the production account. The portable refit reads the archive members directly. |
| [production-sources.tar.gz](data/production-sources.tar.gz) | Frozen site-specific production and extraction scripts, detector source, example cards, seed settings, the fit used during the campaign, and the earlier baseline and four-lepton analysis recipes. Use the student command for the portable calculation. |
| [Four-lepton mass inputs](data/fourlepton-spectrum.json) | Per-process bins and source covariances projected from the frozen original and extended scores. The output CSV includes the three heavy-flavour processes with zero central rate. |
| [GEM table excerpt](data/gem-table18.csv) | The two neighboring mass points from Table 18. The plotting code converts each quoted significance to a 5σ luminosity at its own mass. |
| [Published outputs](published/) | PNG and PDF figures, CSV tables and the JSON result. The CSV mass-spectrum entries are fb per one-GeV bin. `background_mc_sigma_fb` describes the total simulated background per bin. `residual_mc_sigma_fb` includes its correlation with the fitted sidebands. Full covariance matrices are in the input packets. |

The replay verifies each packet against the original fit manifest and checks both the bin edges and source independence. The physics checks cover signed rates, fitted spurious strengths, their uncertainties, all three diphoton shapes, both combined background scenarios and the blinded-window χ². Floating-point comparisons allow 3 × 10⁻⁶ absolute and 2 × 10⁻⁶ relative differences across libraries. The four-lepton calculation separately reconstructs the rates from the per-process test yields at 100 fb⁻¹ and the joint zero-count upper bound from the three unbiased samples. The single-channel Poisson test searches the successive critical counts for the first luminosity with at least 50% power. [Portable refit](reproduce.py) · [Likelihood functions](likelihood.py) · [Frozen classifier test procedure](../hfourlepton/extend_xgboost_test.py)

## Event-level data on Timur

The release lives at `/home/shared/HiggsSSC/student-results/ssc40-discovery-20261123`. Students can read the snapshots and write results into their own checkouts. The older channel event sets remain the starting point for per-object plots, cuts and classifier training. The discovery mass spectrum additionally pools the large background campaign. [Channel setup and student commands](../README.md#analyze-the-prepared-40-tev-events-on-timur)

| Path | Contents |
| --- | --- |
| `/home/shared/HiggsSSC/diphoton/ssc40-20261118` | Baseline signal, prompt γγ, γ + jet and Drell–Yan event data, cards and normalization. Git also supplies the baseline loop and dijet event tables. |
| `/home/shared/HiggsSSC/fourlepton/ssc40-20261118` | Student four-lepton event release with the three heavy-flavour samples and both detector responses. |
| `/home/shared/HiggsSSC/fourlepton/ssc40-20261112/reference/{xgboost,extended_test}` and `supplemental/` under the same release | Frozen model, original and additional scores, extended-test summary, and the additional qq-continuum events. The discovery fit uses this independent extension. |
| `sources/sideband_20pct_20261120/` under the discovery release | Final campaign with `control/final_fit.json`, `control/final_reach.json`, the fit baseline, individual packets, hard-event parents, weighted response trees, dijet pair tables where produced, per-shard spectra, seeds, run cards and repair records. |
| `sources/sideband_campaign_20261119/`, `sources/sideband-20261119/` and `sources/earlier_ho/` | Earlier fit inputs and surviving generator and response records from the previous production waves. |
| `generator-details.tar.gz` and `source-files.sha256` | The first file contains native dijet normalization logs, loop Born-mass maps and numerical NLO integration-grid files. The second records the SHA-256 of each file in the shared `sources/` tree. |
| `analysis/` and `records/` | The portable analysis, figures and input bundle, plus the three standalone auxiliary records needed to rebuild the bundle. |

To check the connection from the independent fit inputs to the final generator production, run the following from your Git checkout on Timur. The first command pools the NLO and dijet shard spectra, checks their covariances and source coverage, and matches the hard and loop-response file identities against the source manifest. The second command reads and hashes the underlying files, which takes longer.

```bash
release=/home/shared/HiggsSSC/student-results/ssc40-discovery-20261123
python discovery/verify_sources.py "$release"
(cd "$release" && sha256sum --quiet -c source-files.sha256)
```

The following commands reconstruct the portable input bundle from the shared originals.

```bash
mkdir -p results/discovery-rebuild
python discovery/build_inputs.py "$release/sources/sideband_20pct_20261120" \
  "$release/records" results/discovery-rebuild/spectra.tar.gz
cmp discovery/data/spectra.tar.gz results/discovery-rebuild/spectra.tar.gz
```

The selected four-lepton mass histogram can also be rebuilt from the original and additional student score files. The projection checks the integrated rate and signed source variance for each of the 27 signal and continuum samples against the frozen extended-test summary.

```bash
python discovery/prepare_fourlepton.py /home/shared/HiggsSSC/fourlepton/ssc40-20261112/reference \
  > results/discovery-rebuild/fourlepton-spectrum.json
```

The clean HTTPS checkout was tested with Python 3.14.7, NumPy 2.5.5, SciPy 1.18.1 and Matplotlib 3.11.2. The independent Timur replay used Python 3.9.25, NumPy 2.0.2, SciPy 1.13.1 and Matplotlib 3.9.4.

### From generated events to spectra

* **NLO prompt γγ and γ + jet.** Each final `nlo/parents/<tag>/` contains the native compressed Les Houches event file, cross section and cards. Its shower shards carry the analysis ROOT trees and signed-weight records. The source denominator uses the signed sum over the complete hard parent. Shower shards are disjoint ranges of that parent. Add their normalized contributions before pooling independent parents. The two processes have different hard cross sections. [Higher-order generator and shower guide](../hgammagamma/HOAnalysis/README.md)
* **Gluon box.** Each run has a hard-event checksum, the reconstructed analysis tree and a Born-region assignment. The focused sample covers 100–150 GeV in Born diphoton mass. The broad sample supplies the complement, including migration into the reconstructed fit range. For broad samples the maps accompany the native event files. [Born-mass mapping](../hgammagamma/born_mass_mapping.py)
* **Dijet double fakes.** Each minimum-bias seed has a native inelastic cross section. Within each event, sum all eligible pair weights by reconstructed mass bin, then form the source outer product. The final production validated and removed the large native minimum-bias ROOT files for ordinary shards. The compact pair rows, moments, cross sections, seeds and per-shard spectra support the mass-spectrum replay. Repaired shards also have recovery metadata. Re-decays and control runs do not add physical exposure. Generating new exclusive observables from unsaved particles requires rerunning the relevant seed with the original generator and response. [Fake-probability and response conventions](../hgammagamma/HOAnalysis/student_backgrounds/README.md)
* **Four leptons.** The shared score files and model reproduce the 0.488 classifier threshold selection without retraining. The baseline estimator combines the original test set with 200,000 independent qq events per flavour, weighted by the relative exposures. The additional 3 × 1 million unbiased heavy-flavour events supply the 95% upper scenario. [Four-lepton event guide](../hfourlepton/SHARED_EVENTS.md) · [Upper-bound study](../docs/ANALYSIS.md)
* **New generation.** Use the repository's supported pipelines in a new run directory and select the same hard process, PDF, detector response and normalization as in the archived cards. The production snapshot records the original host paths, CPU assignments and seed schedules for inspection. Absolute paths in those archived orchestration scripts refer to the original installation. Use new independent hard seeds when extending the statistics. Preserve the source-event weights, inclusive event counts and phase-space split before combining a new sample with the published one. [Diphoton production](../hgammagamma/HOAnalysis/PIPELINE_STUDENT_GUIDE.md) · [Four-lepton production](../hfourlepton/HIGHER_ORDER.md)
