# Electroweak continuum at the SSC

The three independent MadGraph5_aMC@NLO 3.5.15 tree-level samples use two 20 TeV proton beams and the NNPDF40_nlo_as_01180_qed central member (LHAPDF ID 335900). This matches the PDF in the tagged-signal shower. The rates use the matrix-element integration for each sample. The Higgs production K-factor plays no role.

| Sample | Hard final state and response | Generator acceptance |
|---|---|---|
| `w_lnu_aa` | Both charges of ℓνγγ, with ℓ = e, μ, τ. Herwig decays the τ. Two prompt matrix-element photons must survive the response. | Each photon has pT > 14 GeV and absolute η < 3.0. The lepton has pT > 13 GeV and absolute η < 3.5. The true photon-pair mass lies between 60 and 170 GeV. |
| `z_ll_aa` | e⁺e⁻, μ⁺μ⁻ and τ⁺τ⁻ accompanied by two photons. The amplitude includes γ* and Z exchange and their interference. | Same photon cuts as for ℓνγγ. At least one lepton must have pT > 13 GeV. There is no generator requirement on the other lepton, its direction, or the dilepton mass. |
| `z_ee_a` | e⁺e⁻γ, with either electron allowed to mimic the second photon. The other electron must independently supply the tag. | The prompt photon and both electrons have pT > 14 GeV and absolute η < 3.0. The photon must be separated from both electrons by ΔR > 0.20. |

The `sm-lepton_masses` restriction keeps the physical electron, muon and tau masses. Both two-photon processes include emission from the charged leptons and interference with the other tree-level diagrams. The neutral current includes off-shell dileptons. Adding separately generated Wγ/Zγ final-state-radiation contributions to these two samples would count collisions twice. The one-photon electron-fake sample is a separate detector hypothesis. The hard processes include incoming u, d, s, c and their antiparticles with diagonal CKM. The nominal calculation omits bottom- and photon-initiated reactions, gluon fusion at one loop and QCD corrections to the hard processes. GEM’s extra 20% on the Zγ components is a separate sensitivity variation. [GEM simulations, section 5.1 and Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

The detector analysis applies its cuts after the QCD and QED shower, MPI, hadronization and unstable-particle decays. Incoming-parton radiation can move events across the generation acceptance. The test in `ew_diagnostics.py` measures the fraction of the reconstructed 125 GeV window within 2 GeV of a generation pT threshold, 0.2 of an η boundary, 0.1 of the eγ ΔR boundary, or 5 GeV of the generated diphoton-mass boundaries. These fractions describe migration *inside* the available generated acceptance. Evaluating the excluded phase space requires an additional integration with looser generation cuts.

The same diagnostic compares the luminosity for massless b quarks to that for d and s with each event's incoming directions, momentum fractions and MG5 factorization scale. The neutral-current matrix elements for the three down-type flavours are identical in the massless approximation. For charged currents the calculation uses the Wolfenstein parameters and CKM elements of MG5’s `sm-ckm` restriction for the u and c rows. The variation keeps the light-quark shower and detector response. It reports the fraction in which the *true* dilepton mass and the mass of the full final state lie within 10 GeV of the Z mass. These truth-level fractions use the blocks in which every LHE event was written to the detector tree, so the original LHE event indices map directly to the saved indices. The central rate uses all complete blocks, including blocks with a documented shower rejection.

## Repeat the generation and response

Load the Herwig 7.3, ROOT and LHAPDF environments. Use a Python 3 environment with Uproot and NumPy. Use the ancestry-enabled HwSim plugin built by the four-lepton pipeline, the response executable built from `response.cc`, and an LHAPDF data directory containing set 335900. The versioned shower template is based on Herwig’s generic LO Les Houches input. It selects angular-ordered QCD and QED showering, uses the hard PDF for the shower and leaves MPI at the installed tune.

Run from the repository root with the input and output paths appropriate to your installation. For the first W integration and all its disjoint shower blocks use

```bash
code=discovery/mass_scan/gem_range/lepton_tag/background_mc
python3 "$code/ew_generate.py" --mg5 "$MG5/bin/mg5_aMC" --process w_lnu_aa \
    --events 16000 --seed 7311001 --cores 14 --work "$work/w_lnu_aa-001"
python3 "$code/ew_campaign.py" --hard "$work/w_lnu_aa-001" --root "$work/showers" \
    --plugin "$plugin" --response "$response" --repo "$PWD" --cores 14 --seed 9300000
```

The two other original integrations used hard seeds 7312001 and 7313001 and allocated 18 and 12 CPUs to Zℓℓγγ and Zeeγ, respectively. The initial Zeeγ smoke block used seed 4730000 and the first 250 collisions. The remaining 500-event block seeds are the start index plus 9500000. The first Zℓℓγγ sample uses start index plus 9400000. The remaining independent integrations run through `ew_series.py` in groups of 48,000 with blocks of 1,250 and new hard and shower seeds. The saved receipts carry each exact interval and seed.

`hard.json` contains the original integrated rate, its MC integration error in pb, beam energies, hard-process identity, seed, cross-section inputs and LHE checksum. The native LHE banner and MG5 process, run and parameter cards are copied to the same directory. For a finished MG5 run with an interrupted receipt step, repeat the original `ew_generate.py` invocation with `--resume-completed`. The recovery validates the immutable event banner. A shower receipt is written atomically after the finite input is exhausted and the response outputs have been written. The first block in each later integration keeps its ROOT file, which permits direct inspection of the detector response. Intermediate ROOT files in the other successfully reduced blocks can be reproduced from the hard source, card, seeds and software identities.

Run the truth diagnostics in the same LHAPDF Python environment with

```bash
python3 "$code/ew_diagnostics.py" "$work" "$work/ew-diagnostics.json"
```

## Complementary generation control

For a direct migration test, specify `--outside-nominal` with the standard process and write to a *separate* campaign root. This runs a different integration. Its events fill only the phase space between the nominal generation cuts and a wider envelope. All photons have pT > 8 GeV and absolute rapidity < 4 in the wider envelope. In ℓνγγ and ℓℓγγ the diphoton mass range expands to 50–180 GeV. In ℓνγγ the charged lepton has pT > 8 GeV and absolute rapidity < 4.5. In ℓℓγγ at least one charged lepton has pT > 8 GeV. In eeγ both electrons have pT > 8 GeV and absolute rapidity < 4, and each electron has ΔR > 0.08 from the photon. MadGraph evaluates the complementary cuts with its native kinematic functions. The producer independently tests each output LHE collision against the original acceptance before creating a completion receipt.

Each physics rate in the wider envelope equals the *sum* of the nominal and complementary rate for that background. The two integrations have different physical cross sections and must be reduced separately. The existing student nominal three-sample analysis can be rerun independently of this diagnostic.
