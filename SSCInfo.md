# SSC Detector And Higgs Channel Notes

These notes treat the GEM detector studies as the main numerical SSC baseline.
The GEM Higgs note is the most complete SSC-era source used here for both
channels: it updated the GEM TDR Higgs analysis, used detailed GEANT
simulations for detector response and lepton/photon identification, and covered
the channels $H\to\gamma\gamma$, $H\to 4\ell$, $\ell\ell\nu\nu$, and
$\ell\ell jj$.

The SDC program is still relevant as a comparator. The SDC TDR and SDC
four-lepton/muon-resolution notes are listed in the Fermilab SSC report index,
but the numerical fake-rate and Higgs-channel tables below are GEM-based
because those quantities are explicitly tabulated in the accessible GEM Higgs
study.

## SSC/GEM Detector-Performance Assumptions

| Object | Geometric acceptance / efficiency | Projected resolution | Mis-tag / fake assumptions |
| --- | --- | --- | --- |
| Photons | EMC coverage $0.1 < \lvert\eta\rvert < 3.0$. For high-resolution $H\to\gamma\gamma$, photons were excluded from the barrel/endcap transition $1.01 < \lvert\eta\rvert < 1.16$. Real photon/electron ID was kept at about 90% at the shower-ID stage; including geometry and cuts, photon/electron efficiency was 80-85%. | EMC energy resolution $6\text{-}8\%/\sqrt{E}\oplus 0.4\%$; position resolution $4.4~\mathrm{mm}/\sqrt{E}$; photon pointing about $40\text{-}50~\mathrm{mrad}/\sqrt{E}$ plus a constant term. The resulting $H\to\gamma\gamma$ mass resolution was 0.66-1.0 GeV for $M_H=80\text{-}160~\mathrm{GeV}$ in the baseline design. | Jet$\to\gamma$ after photon ID: quark jets $4.1\to 3.4\times 10^{-4}$ and gluon jets $1.2\to 0.83\times 10^{-4}$ over $M_{\gamma\gamma}=80\text{-}160~\mathrm{GeV}$. Electron$\to\gamma$ was studied separately: $R(\gamma/e)$ ranged from 5.0% to 0.14% depending on the tracker-hit cut; the analysis assumed about 0.15% near $M_Z$ and 2% elsewhere. |
| Electrons | Same EMC fiducial logic as photons. High-resolution analyses excluded $1.01 < \lvert\eta\rvert < 1.16$; overall electron/photon ID efficiency after geometry and cuts was 80-85%. | Electron energy resolution follows the EMC resolution. In the $H\to ZZ^\ast\to 4e$ study, the reconstructed Higgs mass resolution was 1.05-1.33 GeV for $M_H=140\text{-}180~\mathrm{GeV}$. | Jet$\to e$ rejection was tabulated as a total rejection factor of 123,750, corresponding to a fake probability of roughly $8\times 10^{-6}$, consistent with the text's "typically $10^{-5}$" electron fake level. |
| Muons | Muon coverage $0.1 < \lvert\eta\rvert < 2.5$. Typical muon ID efficiency was about 80% within the fiducial volume, including 85% geometrical acceptance from holes/cracks and 95% reconstruction efficiency. | Muon momentum resolution at $p_T=500~\mathrm{GeV}$ was projected as 5% at $\eta=0$ and 12% at $\lvert\eta\rvert=2.5$. In $H\to 4\ell$, $4\mu$ mass resolution was 1.59-2.22 GeV and $2e2\mu$ was 1.36-1.77 GeV for $M_H=140\text{-}180~\mathrm{GeV}$. | GEM did not tabulate a universal jet$\to\mu$ fake probability. It stated that the calorimeter reduced hadron punch-through to well below genuine-muon rates; reducible lepton backgrounds were evaluated at event level through $ZQ\bar{Q}$ and $t\bar{t}$ samples with isolation and track matching. |
| Jets | HCAL coverage $0.1 < \lvert\eta\rvert < 5.5$. Central-jet ID efficiency was not the key object in $H\to\gamma\gamma$ or $4\ell$; jets mainly enter as fake photons/leptons or as forward tags. | Central jet energy resolution: $60\%/\sqrt{E}\oplus 4\%$. A forward-tagging study for very heavy $H\to ZZ\to 4\ell$ assumed forward calorimeter segmentation $0.2\times 0.2$ and forward-jet energy resolution $200\%/\sqrt{E}\oplus 6\%$. | Jet$\to\gamma$ and jet$\to e$ rates are listed above. For heavy-Higgs forward tagging, requiring one jet with $E>0.5~\mathrm{TeV}$ and $2 < \lvert\eta\rvert < 5.2$ tagged 60% of signal and 27% of background, increasing purity but reducing statistics. |

## Scenario A: $H\to\gamma\gamma$

The SSC/GEM study treated $H\to\gamma\gamma$ as the essential low-mass
discovery channel. It explicitly states that below about $M_H<130~\mathrm{GeV}$
the $4\ell$ mode is not sufficient and $H\to\gamma\gamma$ is the rare channel
with adequate signal-to-background and significance. The key detector
requirement is a precision EMC with strong photon ID.

### Selection And Efficiency

The analysis used isolated photon candidates, EMC shower-shape discrimination,
HCAL veto information, tracker-based electron/photon separation, and a
$\lvert\cos\theta^\ast\rvert$ cut.

After the final cuts, the listed accepted signal cross sections correspond to
roughly 18-27% of the tabulated $H\to\gamma\gamma$ production times branching
rate values across the mass points. GEM also reported a 98.8% trigger
efficiency for events passing the photon selection.

### Mass Resolution

The baseline $H\to\gamma\gamma$ mass resolutions in GEM were:

| Quantity | 80 GeV | 100 GeV | 120 GeV | 140 GeV | 160 GeV |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baseline $\Delta M_H$ [GeV] | 0.66 | 0.77 | 0.84 | 0.91 | 1.0 |
| Degraded EMC case [GeV] | 1.3 | 1.5 | 1.7 | 1.8 | 2.0 |

The degraded case corresponds to $a=14/17$, $b=1.0$, compared with the baseline
$a=6/8$, $b=0.4$ energy-resolution assumptions.

### Backgrounds And Fake Rates

The irreducible background is prompt $\gamma\gamma$ production. Reducible
backgrounds are $\gamma j$, $jj$, and Drell-Yan $e^+e^-$ with electrons
misidentified as photons.

GEM's jet$\to\gamma$ study generated $10^8$ $\gamma q$ and $\gamma g$ events,
applied isolation, and then detailed GEANT shower-shape ID. After photon ID,
quark-jet fake rates were approximately $3.4\text{-}4.1\times 10^{-4}$;
gluon-jet fake rates were approximately $0.8\text{-}1.2\times 10^{-4}$.

Electron$\to\gamma$ was controlled with central-tracker hit/no-hit information.
GEM needed about $R(\gamma/e)\sim 0.5\%$ near the $Z$ peak to keep Drell-Yan
below 20% of the Higgs signal, and it assumed 0.15% near $M_Z$ and 2%
elsewhere.

### Discovery Reach

In the summary sensitivity table, $H\to\gamma\gamma$ alone was projected at:

| $M_H$ [GeV] | Luminosity | Projected significance |
| ---: | ---: | ---: |
| 100 | $10~\mathrm{fb}^{-1}$ | $4.6\sigma$ |
| 120 | $10~\mathrm{fb}^{-1}$ | $7.8\sigma$ |
| 140 | $10~\mathrm{fb}^{-1}$ | $9.0\sigma$ |
| 150 | $10~\mathrm{fb}^{-1}$ | $7.3\sigma$ |
| 80-90 | $30~\mathrm{fb}^{-1}$ | $3.9\text{-}4.9\sigma$ |

### Late-1990s/Early-2000s Relevant References

By the 1999 ATLAS Physics TDR, the same experimental logic was standard:
$H\to\gamma\gamma$ was described as a narrow mass peak over continuum
$\gamma\gamma$, requiring excellent EM energy/angular resolution, photon/jet
rejection, and photon/electron separation.

ATLAS used an 80% photon-ID efficiency, photon $p_T$ cuts of 40 and 25 GeV,
$\lvert\eta\rvert<2.4$, a transition-region veto, and mass resolutions of about
1.1-1.7 GeV over 80-150 GeV. The TDR also quantified the required large
$\gamma/j$ rejection and an electron-veto inefficiency as low as 0.19% for
$Z\to ee$ rejection.

Phenomenologically, weak-boson-fusion $H\to\gamma\gamma jj$ with energetic
forward jets was already proposed in 1997 as a low-background alternative using
parton-level matrix-element studies.

## Scenario B: $H\to ZZ^{(\ast)}\to 4\ell$

GEM called $H\to ZZ/ZZ^\ast\to 4\ell$, with $\ell=e,\mu$, the cleanest SSC
Higgs signal because of the four isolated leptons. For $M_H>2M_Z$, both lepton
pairs can be constrained to $M_Z$; below threshold, the narrow reconstructed
four-lepton peak is the main discriminator.

### Selection And Efficiency

For intermediate masses $140<M_H<2M_Z$, the cuts were:

- four leptons with $\lvert\eta\rvert<2.5$;
- $p_T>10~\mathrm{GeV}$;
- electron crack veto $1.01 < \lvert\eta\rvert < 1.16$;
- lepton isolation $R=0.35$, $E_T^{\rm cut}=5~\mathrm{GeV}$;
- lepton ID and track matching;
- dilepton mass windows $10<M_{\ell\ell}^{\rm low}<100~\mathrm{GeV}$ and
  $70<M_{\ell\ell}^{\rm high}<100~\mathrm{GeV}$.

Offline trigger efficiency after cuts was above 98% for $4e$ and above 99% for
the other modes.

### Reducible Backgrounds

The important backgrounds were irreducible $ZZ/ZZ^\ast$, $ZQ\bar{Q}$ with
$Q=b,t$, and $t\bar{t}$, where $b$ jets can produce isolated or fake leptons.

GEM found that lowering the lepton threshold from 10 to 5 GeV increased jet
background by a factor of 3-4 for only a 33% signal gain, while increasing the
isolation cone from $R=0.30$ to $R=0.35$ reduced jet background by more than 50%
at a 7% signal cost.

### Mass Resolution And Reach For $140\text{-}180~\mathrm{GeV}$

| Quantity | 140 GeV | 150 GeV | 160 GeV | 170 GeV | 180 GeV |
| --- | ---: | ---: | ---: | ---: | ---: |
| $4e$ mass resolution [GeV] | 1.05 | 1.06 | 1.13 | 1.23 | 1.33 |
| $4\mu$ mass resolution [GeV] | 1.59 | 1.62 | 1.73 | 1.84 | 2.22 |
| $2e2\mu$ mass resolution [GeV] | 1.36 | 1.46 | 1.56 | 1.71 | 1.77 |
| Combined significance, $10~\mathrm{fb}^{-1}$ | $11\sigma$ | $13\sigma$ | $8.1\sigma$ | $5.7\sigma$ | $10\sigma$ |

These numbers are from GEM Table 22.

### Heavy Higgs Reach

For $M_H=200,400,600,800~\mathrm{GeV}$, GEM projected accepted signal cross
sections of 21, 14, 4.3, and 1.3 fb against $ZZ$ backgrounds of 3.0, 2.3, 1.0,
and 0.6 fb.

The corresponding $10~\mathrm{fb}^{-1}$ significances were $38\sigma$,
$28\sigma$, $9.7\sigma$, and $4.3\sigma$. At 800 GeV, $30~\mathrm{fb}^{-1}$
raised the significance to $7.5\sigma$, though systematics on the broad $ZZ$
background became important. Forward jet tagging was tested for the 800 GeV
case and gave a higher-purity but lower-statistics sample.

### Late-1990s/Early-2000s Relevant References

The 1999 ATLAS TDR used the same core $4\ell$ strategy: optimized $Z$-mass and
off-shell dilepton mass windows, four-lepton mass reconstruction, lepton
isolation, vertex/impact-parameter rejection of heavy-flavor backgrounds, and
mass-window counting.

Its $H\to ZZ^\ast\to 4\ell$ kinematic acceptance was about 27-54% for
$M_H=120\text{-}180~\mathrm{GeV}$, and it used a nominal 90% per-lepton
ID/reconstruction efficiency in significance estimates. Full-simulation
resolutions at $M_H=130~\mathrm{GeV}$ ranged from about 1.42 GeV in $4\mu$ to
1.81 GeV in $4e$ at high luminosity.

For heavy $H\to ZZ\to 4\ell$, ATLAS emphasized a harder-$Z$ $p_T$ cut and
possible forward-jet tagging, with reach up to roughly 800 GeV.

## Signal-Vs-Background Discrimination Techniques Available By The Late 1990s / Early 2000s

For an SSC-style recast, these are the techniques that would be historically
defensible:

| Technique class | $H\to\gamma\gamma$ use | $H\to 4\ell$ use | References |
| --- | --- | --- | --- |
| Cut-based optimized selections | Photon $p_T$, $\eta$, isolation, transition-region veto, $\lvert\cos\theta^\ast\rvert$, mass window around $M_{\gamma\gamma}$. | Lepton $p_T$, $\eta$, isolation, track matching, $Z$-mass constraints, off-shell dilepton thresholds, $4\ell$ mass window. | GEM and ATLAS TDR both used this as the baseline. |
| Shower-shape and object-ID likelihoods | EMC strip/tower shower-shape likelihoods, HCAL veto, photon pointing, conversion/electron-veto logic. | Electron/muon ID, track matching, isolation; later ATLAS also used impact-parameter/vertexing to reject $t\bar{t}$ and $Zb\bar{b}$. | GEM used shower-shape likelihoods and GEANT-based ID; ATLAS TDR quantified photon/jet and photon/electron rejection. |
| Mass-peak / sideband logic | Narrow $M_{\gamma\gamma}$ peak over continuum; resolution directly drives significance. | Narrow $M_{4\ell}$ peak below threshold; $Z$-mass constraints above threshold. | GEM explicitly relates $H\to\gamma\gamma$ significance to mass resolution and uses $M_H\pm 2\Delta M_H$ in $4\ell$. |
| Likelihood-ratio and $CL_s$ statistics | Applicable to binned $M_{\gamma\gamma}$ or multichannel combinations. | Applicable to $4\ell$ mass windows and multichannel $4e/4\mu/2e2\mu$ combinations. | LEP Higgs searches used likelihood-ratio ordering and $CL_s$; Read emphasized applicability to counting, multichannel, and multidimensional discriminants. |
| Neural networks / multivariate discriminants | Not part of the original GEM baseline, but historically available by the late 1990s and early 2000s. | Especially useful where kinematics, isolation, impact parameters, and mass variables can be combined. | LEP/L3 analyses combined event likelihood or neural-network output with reconstructed Higgs mass to build final discriminants; ATLAS studies subsequently used neural-network methods for VBF Higgs channels. |
| Forward-jet tagging / central-jet veto ideas | Useful for VBF $H\to\gamma\gamma jj$, not the original inclusive GEM $\gamma\gamma$ baseline. | Tested in GEM for very heavy $H\to ZZ\to 4\ell$; important for VBF/heavy-Higgs phenomenology. | GEM forward-tag study; Rainwater-Zeppenfeld $H\to\gamma\gamma$ in WBF; Rainwater-Zeppenfeld-Hagiwara $H\to\tau\tau$ WBF showing forward jets as background suppression. |
| Angular and polarization observables | $\lvert\cos\theta^\ast\rvert$ was already used in GEM $H\to\gamma\gamma$. | Harder-$Z$ $p_T$, $Z$-polarization/decay-angle ideas, and mass constraints. GEM tried $Z$ polarization for 800 GeV but found limited discrimination after cuts/statistics. | GEM $H\to\gamma\gamma$ and heavy $4\ell$ sections; ATLAS heavy $H\to ZZ$ used harder-$Z$ $p_T$. |

## Practical Bottom Line

For an SSC-era discovery study, a historically faithful baseline would use
GEM-like detector performance, cut-based selections, GEANT-derived object
ID/fake rates, and Poisson or likelihood-ratio significances.

A late-1990s/early-2000s upgraded analysis could defensibly add LEP-style
likelihood-ratio/$CL_s$ combinations, event-likelihood or neural-network
discriminants, and VBF forward-tag categories. Boosted decision trees and the
modern matrix-element discriminants used in later LHC Higgs analyses would be
less faithful to the immediate late-1990s baseline.

## References

### SSC-Era Detector And Higgs Studies

- GEM Collaboration, *GEM Technical Design Report*, SSC-GEM-TN-93-262,
  SSCL-SR-1219, April 1993. Used here as the detector-performance baseline.
  Fermilab Library Server:
  <https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-262.pdf>.
- GEM Collaboration, *Higgs Boson Studies with the GEM Detector at the SSC*,
  SSC-era GEM Higgs analysis note/update to the GEM TDR. This is the numerical
  source for the GEM $H\to\gamma\gamma$ and $H\to ZZ^{(\ast)}\to 4\ell$ tables
  summarized above; the exact report number should be filled in from the
  archive copy used.
- SDC Collaboration, *Solenoidal Detector Collaboration Technical Design
  Report*, SSC-era detector TDR. Used here as a comparator rather than the
  numerical source; the exact report number should be checked against the
  Fermilab SSC report index.

### LHC TDR And Higgs Search Methodology

- ATLAS Collaboration, *ATLAS detector and physics performance: Technical
  Design Report, 1*, CERN-LHCC-99-014, ATLAS-TDR-14, CERN Document Server:
  <https://cds.cern.ch/record/391176>.
- ATLAS Collaboration, *ATLAS detector and physics performance: Technical
  Design Report, 2*, CERN-LHCC-99-015, ATLAS-TDR-15, CERN Document Server:
  <https://cds.cern.ch/record/391177>.
- LEP Higgs Working Group for Higgs Boson Searches and the ALEPH, DELPHI, L3,
  and OPAL Collaborations, "Search for the Standard Model Higgs Boson at LEP,"
  arXiv:hep-ex/0306033: <https://arxiv.org/abs/hep-ex/0306033>.
- L3 Collaboration, "Search for the Standard Model Higgs boson in $e^+e^-$
  collisions at $\sqrt{s}$ up to 202 GeV," arXiv:hep-ex/0012019:
  <https://arxiv.org/abs/hep-ex/0012019>.
- A. L. Read, "Presentation of search results: the $CL_s$ technique,"
  *J. Phys. G* **28** (2002) 2693-2704,
  DOI: `10.1088/0954-3899/28/10/313`.

### Weak-Boson-Fusion Phenomenology

- D. Rainwater and D. Zeppenfeld, "Searching for $H\to\gamma\gamma$ in weak
  boson fusion at the LHC," arXiv:hep-ph/9712271:
  <https://arxiv.org/abs/hep-ph/9712271>.
- D. Rainwater, D. Zeppenfeld, and K. Hagiwara, "Searching for
  $H\to\tau\tau$ in weak boson fusion at the LHC," arXiv:hep-ph/9808468:
  <https://arxiv.org/abs/hep-ph/9808468>.
