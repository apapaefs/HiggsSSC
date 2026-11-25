# Extra leptons in inclusive diphoton backgrounds

GEM's lepton-tagged background model includes heavy-flavour production with two photons, radiative W and Z production, and top events with misidentified photons. The study found the isolated bottom-pair plus diphoton contribution negligible. It did not tabulate separate rates for inclusive diphoton, photon + jet or ordinary dijet events with an extra nonprompt or misidentified lepton. The fraction of the inclusive continuum represented by GEM's tagged sample only measures the scale of possible event overlap. The rate at which other inclusive events acquire a lepton tag requires a separate control. [GEM's tagged background list and treatment](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=59)

We matched the original one-million-event prompt diphoton and photon + jet HwSim trees to the signed, weighted photon response by source-event index. The control requires two selected photon candidates with pair transverse momentum above 40 GeV. An additional electron or muon must have transverse momentum above 20 GeV and absolute pseudorapidity below 2.5. We also checked the angular distance to the nearest stored jet, requiring ΔR > 0.3. The native trees contain reconstructed objects, with jets stored above 10 GeV. The stable-particle information needed for GEM's 5 GeV isolation-energy cut in a cone of radius 0.3 is unavailable in these trees. The check assigns unit lepton identification efficiency.

| Source sample | Selected photon-pair events, 123.3–126.7 GeV | With a fiducial lepton | With a lepton separated from jets | Cross section before lepton selection [fb] |
| --- | ---: | ---: | ---: | ---: |
| Prompt γγ | 130 | 1 | 0 | 399.80 ± 35.06 |
| Photon + jet | 73 | 0 | 0 | 158.97 ± 25.55 |

In the broader 100–150 GeV interval eight prompt diphoton and five photon + jet events have a fiducial lepton. All have a jet within ΔR = 0.3. A single candidate in these samples can carry about 3 fb. GEM's background in the 125 GeV window is 0.408 fb, so the available statistics and approximate isolation cannot establish that an additional contribution is small. These two processes have no electron already assigned to a photon through an electron-fake response. A corresponding Drell–Yan control would have to identify a third lepton after removing both candidate electrons.

As a scale estimate, GEM found an isolated jet-to-electron misidentification probability of 1/123,750 in its 40 GeV quark-jet study. Applying that value to every stored jet above 20 GeV within |η| < 2.5 in the two selected samples gives approximately 0.0064 fb in the 125 GeV window. This deliberately also counts the jet used as a fake photon in photon + jet events. For the dijet double-fake baseline, 31 pairs from 30 collisions give 285.69 ± 52.93 fb after the same kinematic cuts and before demanding a lepton. The complete unordered-pair table determines the number of additional jets passing the photon kinematic acceptance in each collision. Removing the two jets already assigned to photons gives 340.64 ± 96.94 fb after weighting by that multiplicity, or another 0.0028 fb with the same jet-to-electron transfer. The sum of these three illustrative contributions is about 0.0092 fb. This uses a single probability for all jet flavours and momenta and provides no upper bound on jets faking electrons or muons or on semileptonic heavy-flavour decays. [GEM electron identification, Table 3](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=19)

The [source-level control output](inclusive-lepton-control.json) records event counts, signed rates, cross-section errors from independent source events, and the native and detector file hashes. On Timur, run the [control script](check_inclusive_leptons.py) with NumPy and Uproot installed. The second command reproduces the [dijet result](dijet-extra-jets.json) from the student checkout with NumPy alone.

```bash
python3 discovery/mass_scan/gem_range/lepton_tag/check_inclusive_leptons.py
python3 discovery/mass_scan/gem_range/lepton_tag/check_inclusive_leptons.py --dijet
```
