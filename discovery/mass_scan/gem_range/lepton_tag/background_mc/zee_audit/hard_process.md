# Hard-process comparison for the electron-fake background

The 1992 GEM study quotes 14.4 pb for dilepton plus photon production with transverse momenta of the hard photon and Z above 10 GeV. It includes a 20% increase to the quark-initiated prediction to approximate gluon fusion. The 1993 report provides the final mass-dependent electron-fake rates. The independent PYTHIA 5.6 calculation reproduces the 1992 hard rate when the dilepton cross section includes both electrons and muons. Only the electron decay contributes to the electron-fake background. [1992 study, §7.1](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-92-126.pdf#page=31) · [1993 study, Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61) · [PYTHIA 5.6 inputs and results](vintage.json)

With 20 TeV proton beams, four light quark flavours, the NNPDF4.0 NLO QED central PDF, and the same 10 GeV production cut, the MG5 tree-level cross section for a stable Z and a photon is **369.51 ± 0.36 pb**. For the separately integrated production and decay amplitude with Z → e⁺e⁻ it is **12.426 ± 0.017 pb**. This integration includes the decaying pole up to 15 Z widths and integrates over the full lepton directions with spin correlations. For a massless muon the prediction is equal up to negligible lepton mass effects. Doubling for two flavours and applying GEM’s 1.2 multiplier gives **29.82 ± 0.04 pb**. The matched historical two-flavour value is 14.4 pb. Omitting the multiplier on the modern calculation gives **24.85 pb**. Both modern rates use the frozen electroweak inputs. [Stable Z integration](hard_stable_z.json) · [Electron final state integration](hard_z_decay.json)

The production sample used for the discovery forecast has a different generator acceptance and the full eeγ tree amplitude. Its two disjoint generator regions have cross sections **24.987 ± 0.018 pb** and **87.875 ± 0.069 pb**. Of these cross sections, respectively **8.51 × 10⁻⁶** and **5.09 × 10⁻⁷** pass the complete reconstructed Higgs-window selection, including the electron-fake probability. The ratio of their sum to GEM’s 14.4 pb would mix different generation cuts, amplitudes and lepton flavours.

## Reweighting the selected events

Each of the 149 collisions contributing to the 125 GeV window can be identified in the original LHE event files. For shower blocks with a missing event, the calculation matches the IDs, longitudinal momenta and energies of both incoming partons to the source LHE. The selected event weights from this lookup agree with the independent PDF analysis. The two generator regions contain 640,000 attempted collisions. The default normalization gives **0.25746 ± 0.02616 fb**.

We export the full tree-level amplitude and the two production-Zγ diagrams as separate MadGraph standalone C++ matrix elements. Both evaluations use the same incoming and outgoing momenta, flavour, spin sums, width and couplings. The underlying PDFs and phase-space factors cancel in their ratio. The two amplitudes for photon emission from the initial quark and antiquark are added coherently. Up and charm have identical matrix elements in the input model, as do down and strange. Reweighting the same collisions keeps the Herwig shower and detector response fixed.

| Change to the hard weight | Selected eeγ background [fb] |
| --- | ---: |
| Original full eeγ model | 0.25746 ± 0.02616 |
| Production Zγ diagrams with Z → ee | 0.20064 ± 0.02388 |
| Also require the hard photon pT > 10 GeV | 0.19448 ± 0.02307 |
| Also restrict the dilepton mass to 15 widths around the Z | 0.18987 ± 0.02304 |
| Also change both incident PDFs to EHLQ1 | 0.11000 ± 0.01328 |
| Also apply the historical 20% gluon-fusion prescription | 0.13200 ± 0.01594 |
| GEM’s 1993 electron-fake result, interpolated to 125 GeV | 0.03740 |

With the full modern domain, coherent initial-state radiation gives 0.21577 ± 0.02420 fb, of which the pure production-Z terms give 0.20064 fb, production through a virtual photon gives 0.01141 fb and their interference gives +0.00371 fb. Radiation from the charged leptons gives 0.02950 ± 0.00566 fb. Interference between initial- and final-state emission gives +0.01219 ± 0.00439 fb. The three coherent contributions sum to 0.25746 fb. The two production-Z diagrams embedded in the full amplitude agree with the separately generated decay-chain amplitude on each original collision.

The 15-width interval matches the explicit MadGraph decay-chain integration. It is a choice for this comparison. Leaving the Breit–Wigner tail open and making all the other changes gives 0.13542 fb. The product of the PDF and matrix-element ratios is evaluated separately for each event. The ratio of the original result to the result with the production-Z definition, 10 GeV cut and 15-width range is **1.356 ± 0.076**. Switching from NNPDF to EHLQ1 on that same selection gives a modern-to-historical PDF ratio of **1.726 ± 0.014**. The product, including the historical 1.2 divisor, is **1.950 ± 0.109**. The final difference from GEM after this substitution is **3.53 ± 0.43**. These errors treat GEM’s rounded published rate as fixed. The quoted Monte Carlo errors include the differing window weights and common cross-section integration errors. Some electroweak input and scale choices of the historical generator also differ. [Paired estimates and covariance](hard_projection.json)

The frozen LHE files record their scales in SCALUP. In a separate 4,000-event core integration with identical settings and per-beam PDF bookkeeping enabled, both factorization scales agree with SCALUP to 5 × 10⁻⁷ relative for every event. The EHLQ1 variation uses these scales for both PDF sets. It does not vary the subsequent shower. [Independent scale check](hard_scale.json)

## Reproduction

Export `MG5` as the path to `bin/mg5_aMC`. The integration requires MG5 3.5.15, the project’s Herwig/LHAPDF environment and the central PDF with ID 335900. The standalone reweighting requires a C++ compiler, Python and NumPy. The checked-in [149 event records](hard_window_events.json.gz) include each source index, four-momenta in GeV, flavour, probability weight from the simulated detector and the EHLQ1/NNPDF ratio from the independent PDF calculation. The checked-in [parameter card](hard_parameters.dat) is the original nominal generator input. From the repository root, reproduce the entire set of paired window rates with

```bash
audit=discovery/mass_scan/gem_range/lepton_tag/background_mc/zee_audit
python3 "$audit/hard_matrix.py" "$audit/hard_window_events.json.gz" "$MG5" \
    "$audit/hard_parameters.dat" "$work/standalone" "$work/hard-projection.json"
```

To repeat the integrations and to rebuild the 149-event input from the raw archives and independently derived PDF weights, use

```bash
code=discovery/mass_scan/gem_range/lepton_tag/background_mc
audit="$code/zee_audit"
bash "$audit/hard_legacy_run.sh" "$work/integrations"
python3 "$code/ew_generate.py" --process z_ee_a --work "$work/scale_control" --mg5 "$MG5" \
    --seed 7342853 --events 4000 --cores 6 --set use_syst=True
python3 "$audit/hard_scale_check.py" "$work/scale_control" "$work/scale-check.json"
python3 "$audit/hard_events.py" "$core_campaign" "$complement_campaign" \
    "$code/final/background.json.gz" "$code/zee-confirmation-spectrum.json.gz" \
    "$work/hard_window_events.json.gz" --showers "$replay_trees" --pdf-weights "$work/pdf-events.json"
```

The core and complement campaign paths each point to their separately extracted raw archive or the original campaign. `replay_trees` is the output of the shower replay. The optional `--weights "$work/pdf-events.json"` argument to the PDF reweighting produces the per-collision input to the combined check. [Frozen raw releases](../final/timur-release.json) · [Shower replay](replay_showers.py) · [Incoming PDF reweighting](pdf_reweight.py)
