# The associated-production diphoton signal

At 40 TeV our showered W⁺H, W⁻H and ttH samples give **0.742 ± 0.033 fb** for a 125 GeV Higgs passing the lepton tag before the mass window. In GEM's ±1.70 GeV window the rate is **0.709 ± 0.032 fb**. Errors describe the finite signal simulation. Applying the historical aggregate acceptance to modern production rates instead gives **0.854 fb** before the window. [Direct measurement](direct.json) · [Historical-acceptance calculation](rates.json)

| Process | Production rate [pb] | Direct selected tag [fb] | MC error [fb] |
| --- | ---: | ---: | ---: |
| W⁺H | 3.27050 | 0.1465 | 0.0080 |
| W⁻H | 2.38910 | 0.1170 | 0.0063 |
| ttH | 5.46877 | 0.4790 | 0.0314 |

We showered 15,000 independent MC@NLO hard events per mode using Herwig 7.3, forcing the stable Higgs to decay into two photons. The source LHE files are the original stable-Higgs production events archived with the SSC four-lepton study. Each mode has three independent generator samples. We processed every hard collision once and applied the H → γγ branching fraction of 0.00227 to the production rate. Herwig decayed the associated tops and W bosons with their own branching fractions. Of the 45,000 attempted collisions, 44,996 produced a saved shower. The signed normalization includes the four failed attempts at zero efficiency. [Exact input hashes, seeds, selected events and response probabilities](events.npz)

For the tag, we identify the two Higgs photons by shower ancestry. We use the SSC photon energy smearing, fiducial acceptance including the calorimeter cracks, and photon identification efficiencies. We require both photons and at least one electron or muon to pass 20 GeV and |η| < 2.5, and reconstructed diphoton transverse momentum above 40 GeV. Particle-level isolation sums the transverse momenta of the other calorimeter-visible stable particles within ΔR < 0.45 around a photon or ΔR < 0.30 around a lepton. The threshold is 5 GeV. This approximates GEM's tower isolation with particles from the same proton collision. The photon smearing includes the project's pileup-noise term. The tag includes the 98.8% GEM trigger efficiency. Muons use the SSC acceptance and identification probability and unsmeared momentum. For more than one lepton we calculate the probability of identifying at least one. [GEM cuts and trigger](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=60) · [SSC detector response](../../../../../hgammagamma/SSC_DETECTOR_RESPONSE.md)

The selected cross section splits into **0.6714 fb** from direct W → e/μ and **0.0710 fb** from W → τ → e/μ. This counts the probability of a tag from a τ daughter after any eligible direct lepton fails identification. The simulated selection gave zero tags from heavy-hadron decays or other lepton origins. The rates describe true WH and ttH signal events. The backgrounds to the tag come from GEM's study.

We apply the existing inclusive signal selection to the same showers and check the photon-pair mass and event weight against the reference SSC executable for every selected collision. The tagged and inclusive photon pairs can differ. Shared photons have a common identification outcome. For a source event accepted by both selections we remove the common probability at its *inclusive* photon-pair mass and assign the tag at its *tagged* mass. The total shared cross section is **0.7086 fb**. The untagged associated contribution is **5.67 fb** between 120 and 130 GeV and includes the photon pairs selected by the inclusive analysis. [Selections and overlap](measure.cc) · [Reference comparison](extract.py)

The mass template bins span 70–170 GeV in 0.05 GeV steps. We smooth the signed histograms with a Gaussian kernel of width 0.15 GeV and clip any remaining negative bins, restoring each category's signed integral in that interval. The full-peak selected rate and window count in the table above are measured directly from the signed event weights. The template file also records the negative signed bin content before smoothing. [Template calculation](summarize.py)

## Replication

From the student checkout the following commands regenerate the portable result using [the discovery dependencies](../../../../requirements.txt).

```bash
sig=discovery/mass_scan/gem_range/lepton_tag/signal
python3 "$sig/transfer.py" > /tmp/ssc-tag-transfer.json
python3 "$sig/summarize.py" "$sig/events.npz" > /tmp/ssc-tag-direct.json
cmp /tmp/ssc-tag-transfer.json "$sig/rates.json"
cmp /tmp/ssc-tag-direct.json "$sig/direct.json"
```

On Timur, the immutable student release contains the raw shower ROOT files, native run cards and logs, per-shard completion records, original LHE file identities and the measurement source. The 317 files have a combined size of 1.33 GB. Its SHA256SUMS file has SHA-256 `7e61561e4ca05ee6595578c5368bfcf324db7aa7c88976dd4dfef2bc71a8d545`. The measure-and-check step requires Uproot, ROOT and a C++17 compiler. It compares the original inclusive executable with the tagged analysis before exporting the per-event table.

```bash
sig=discovery/mass_scan/gem_range/lepton_tag/signal
release=/home/shared/HiggsSSC/student-results/ssc40-tag-diphoton-20261125
work=$(mktemp -d)
(cd "$release" && sha256sum --quiet -c SHA256SUMS)
c++ -O2 -std=c++17 -I"$PWD" $(root-config --cflags) "$sig/measure.cc" -o "$work/measure" $(root-config --libs)
c++ -O2 -std=c++17 $(root-config --cflags) hgammagamma/LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC.cc -o "$work/inclusive" $(root-config --libs)
python3 "$sig/extract.py" "$release" "$work/measure" "$work/inclusive" "$work/events.npz"
python3 "$sig/summarize.py" "$work/events.npz" > "$work/direct.json"
```

To regenerate the raw showers, use the same stable Higgs LHE parents and an ancestry-enabled HwSim plugin. The release records the plugin hash `94f7b191523cba59d774bf8711afdc3caec32afa2d023f5a52a4a2819134417a`. Each independent 5,000-event file from the student archive is split into five disjoint groups of 1,000 collisions. The shower script sets unique seeds for all 45 jobs and verifies the finite-source termination and signed weight denominator. The script copies the archived matching, PDF, MPI and shower settings, changes the Higgs decay and enables ancestry. The photon smearing seed for each job is its shower seed plus 99901. The integer seeds are in the compact event file. [Ancestry plugin source and build instructions](../../../../../hfourlepton/hwsim/README.md)

```bash
module load herwig/stable
sig=discovery/mass_scan/gem_range/lepton_tag/signal
release4l=/home/shared/HiggsSSC/fourlepton/ssc40-20261118
work=$(mktemp -d)
python3 "$sig/transfer.py" --manifest "$release4l/manifest_ssc.json"
python3 hfourlepton/hwsim/build_plugin.py --source-repo /home/shared/Herwig/src/Herwig-7.3.0/Contrib/hwsim --output-dir "$work/plugin"
python3 "$sig/shower.py" "$release4l/nominal/events" "$work/plugin/HwSimAncestry.so" "$PWD" "$work/showers" --workers 24
```

The historical-acceptance scenario uses the same production cross sections, the published W → ℓν branching fraction for each charged lepton species, the probability of at least one direct electron or muon from ttH, and the interpolated historical tagged efficiency of 0.700/6.225. That efficiency is shared by WH and ttH. GEM published their combined cut flow. [GEM Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61) · [W branching fraction](https://pdglive.lbl.gov/BranchingRatio.action?pdgid=S043.7)
