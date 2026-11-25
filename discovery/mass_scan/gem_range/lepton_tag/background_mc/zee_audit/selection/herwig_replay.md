# GEM isolation variables on the tagged Herwig events

The input is the complete 640,000-collision eeγ sample in two independent generator regions. Each of the 498 original Herwig 7.3 shower blocks is read before isolation. All three truth-matched electromagnetic objects must pass the published fiducial acceptance and smeared 20 GeV thresholds, the reconstructed candidate diphoton mass must lie in 123.3–126.7 GeV, and its transverse momentum must exceed 40 GeV. This also admits collisions that failed the previous isolation selection. The alternate selection starts from 267 candidate assignments in 266 independent collisions. Of these, 150 assignments in 149 collisions pass the original scalar isolation. The two charge assignments give exclusive detector outcomes. The program checks that all contributing events have a single eligible genuine photon and two eligible charged electrons. [Nominal selection and packet validation](../detector-replay.md)

Each record in the compressed JSON is a possible assignment, identified by the original shower directory, saved tree entry and reconstructed mass. Records for the same collision must be summed before calculating a variance. The collision numbers have gaps relative to the source LHE in shower blocks that lost an event. The paired hard-weight calculation handles that mapping. The nominal weight includes the genuine photon identification, electron tag, electron fake and trigger. The three-object weight before isolation reduces exactly to 0.988 × (0.90 × 0.96) × 0.90 × (0.90/33). Each originally selected mass, event identity and probability must agree with the archived packet before this file can be written. [Hard weights](paired_weights.md)

## Calorimeter and charged tracks

Stable electrons and photons deposit their entire energy in one electromagnetic tower. Other visible stable particles apart from muons deposit their energy in a hadronic tower. Each compartment sums energies within a tower before applying the cell threshold. Each cell's transverse energy is its energy divided by the hyperbolic cosine of its center pseudorapidity. The truth direction of the particle seeds the cluster. The cone uses a radius of 0.45 for each photon candidate and 0.30 for the tag, measured to tower centers around that particle. Only the EM deposits in the 5 × 5 EM core are subtracted. The primary candidate is always identified using the same smeared energy and direction as the published simulation.

The EM pitches used to approximate the grid of GEMFAST Fig. 3 are 0.026, 0.037, 0.052 and 0.075 in successive absolute pseudorapidity intervals 0–1.56, 1.56–1.90, 1.90–2.20 and 2.20–2.66. A pitch of 0.11 covers the edge of isolation cones outside 2.66. Hadronic cells use three times the corresponding nominal pitch. The numbers of rings and azimuthal bins are rounded to span each interval and the full azimuth. At an interface, a 5 × 5 core includes the central EM ring and two rings on each side and the same azimuthal span as five cells of the seed's ring. The HCAL veto uses hadronic tower centers falling within this core. A uniform 0.026 model with equal EM and hadronic pitch supplies a segmentation variation. The digitized cell geometry, pointlike deposits and common reconstruction are approximations to the GEM detector response. [GEM barrel tower dimensions, Fig. 5](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=23) · [GEMFAST segmentation, Fig. 3](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-376.pdf#page=9)

Each of the genuine photon, fake electron and tag has its own object record. The compressed feature table has SHA-256 `392abbc7426ea20df3e1ba767ec464b71075c7901174b345bc580ce45dae854a`. [Per-assignment feature table](herwig_candidates.json.gz)

| Field | Meaning |
|---|---|
| `true_pt, true_eta, true_phi` and `pt, eta, phi` | Stable-particle and project-smeared kinematics |
| `ieta, iphi, primary_et` | Seed cell indices and transverse energy of the primary in the seed cell |
| `cone` | Original scalar stable-particle transverse momentum in the cone, excluding the candidate, neutrinos and muons |
| `core_iso` | EM plus hadronic transverse energy in the calorimeter cone after 5 × 5 EM-core subtraction |
| `cell_iso_et` and `cell_iso_e` | Same quantity with a 0.5 GeV cell transverse-energy or energy threshold in each compartment |
| `core3_em_et, core_em_et` | Extra EM transverse energy in the 3 × 3 and 5 × 5 regions. Add the chosen containment fraction of the primary energy for a shape estimate |
| `core_had_et, core_had_e` | Hadronic transverse energy and energy in cells mapped into the EM core |
| `extra_track3_pt, extra_track_pt` | Highest transverse momentum of another charged particle in the 3 × 3 or 5 × 5 EM core |
| `extra_em_pt, extra_em_dr, core_merged` | Highest transverse momentum of another EM particle in the 5 × 5 core, its angular distance, and the number of such particles |
| `mean_em_angle_rad, mean92_em_angle_rad` | Energy-weighted mean three-dimensional opening angle of all electrons and photons, including the primary, about their energy-vector sum. The first uses the 5 × 5 EM core, and the second uses the 1992 three-by-three footprint with a uniform pitch of 0.04 |
| `uniform` | Corresponding quantities with the uniform grid |

The charged-track check projects stable particles from their production directions. It does not simulate magnetic bending, conversions, track efficiency or showers in detector material. For the angular observable, a clean isolated electron or photon gives zero and the 1992 threshold is 0.005 radians. It measures the opening of different electromagnetic particles before interactions in detector material. [GEM 1992 angular definition](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf#page=21) · [Application to the tagged channel](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf#page=33)

Cell isolation and sums of neighboring particle energies can be varied independently of intrinsic single-shower containment. GEANT strip likelihoods, longitudinal shower fluctuations and cell-by-cell noise correlations require additional detector information. [GEM object response and sources](gem-response.md)

## Rebuilding the feature table

First follow the [archived-shower instructions](../detector-replay.md) to produce the trees and manifests. From the repository root with ROOT and a C++17 compiler available, run

```bash
audit=discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit
release=/home/shared/HiggsSSC/student-results/ssc40-tag-background-20261125
work="$PWD/zee-detector-replay"
g++ -O2 -std=c++17 -I"$PWD" $(root-config --cflags) "$audit/selection/herwig_cells.cc" \
  $(root-config --libs) -o "$work/herwig_cells"
python3 "$audit/selection/herwig_replay.py" \
  "$release/background.json.gz" "$release/zee-confirmation-spectrum.json.gz" "$work/showers" \
  "$work/herwig_cells" "$work/herwig_candidates.json.gz" --workers 24
```
