# Top backgrounds in the tagged diphoton search

The hard events are generated in proton–proton collisions at 40 TeV with MadGraph5_aMC@NLO 3.5.15. The matrix elements are tree level. The hard process and the Herwig 7.3 shower both use the NNPDF4.0 NLO QED PDF with LHAPDF ID 335900. The proton includes the gluon and the five massless quark flavours and their antiquarks. The calculation includes the leading QCD top production amplitudes. The hard process excludes the Higgs signal.

| Independent hard sample | Matrix-element photons used for the two photon candidates | Remaining candidates |
|---|---|---|
| `tt_2a` | Two emitted in production | None |
| `tt_1a` | One emitted in production | An electron or jet misidentified as a photon |
| `tt_0a` | Zero | Two electrons or jets misidentified as photons |
| `tt_1a_1dec` | One emitted in production and one in a top decay | None |
| `tt_0a_1dec` | One emitted in a top decay | An electron or jet misidentified as a photon |
| `tt_0a_2dec_same` | Both emitted in one top decay | None |
| `tt_0a_2dec_split` | One emitted in each top decay | None |

The first three samples also provide an independent control in which Herwig models photons in the decay chains. The standard shower includes QCD and QED splittings. The installed top and W decayers assign the SOPHTY photon generator. The shower-only control and the resolved-decay estimate share events. Use the resolved combination option in the summary script to replace all shower-only decay-photon outcomes before adding the four decay matrix elements.

An optional importance sample, `tt_0a_1dec_ll`, explicitly requires a direct electron or muon and its associated neutrino in each top decay. The photon is emitted in either top decay. Its native hard cross section includes both leptonic decay rates. In the dilepton-enriched combination, this integration supplies only the electron-fake outcome for those two direct leptonic top decays. The inclusive single-photon decay sample supplies all its jet fakes and the electron-fake outcomes for other top decays, including tau and hadronic modes.

For a radiating top, the decay amplitude is `t > b f f a` or `t > b f f a a`, and similarly for the antitop. Charge conservation and the SM couplings select the two light quark pairs, electron, muon and tau channels of the W. Both charges are included. This factorizes production and decay at the top pole, while keeping the contributing charged-current diagrams and their interference within each radiative top decay. The unstable W is internal to that amplitude. Interference between production and decay or between competing assignments of photons to two resonant tops is beyond this pole approximation. Any other on-shell top decays inclusively in Herwig.

The radiative samples use physical charged-lepton masses and allow zero photon–lepton separation. A muon near a photon can pass the calorimeter isolation if a different lepton tags the event. The default hard cuts require photon transverse momentum above 14 GeV and absolute pseudorapidity below 3.5. The two-matrix-element-photon samples require a hard diphoton mass between 60 and 170 GeV to cover the reconstructed 80–140 GeV scan. In the radiative top samples the quark–photon separation is above 0.05. This regulates the collinear massless-quark singularity. Dependence on this cut, quark-to-photon fragmentation, and higher-order QCD and electroweak effects are model uncertainties. The hard cuts and the detector selection are recorded separately.

For each independent integration the event rate in fb is the LHE cross section in pb times 1,000 times the summed detector weight divided by the signed number of hard events attempted. The LHE cross sections for the radiative decay processes already include the corresponding matrix-element decay rates, with the top and W widths recorded in `param_card.dat`. The electron-to-photon transfer is 0.90/33. The showered and smeared object used for this transfer cannot simultaneously supply the lepton tag. The inclusive fake rates and the jet-flavour variation are implemented in the shared detector response.

The first three hard integrations used seeds 721000–721002 and gave the following cross sections before detector selection. The diphoton mass cut applies to the two-photon sample.

| Sample | Hard events | Cross section [pb] |
|---|---:|---:|
| `tt_0a` | 30,000 | 4083.590 ± 3.033376 |
| `tt_1a` | 30,000 | 16.77874 ± 0.03132324 |
| `tt_2a` | 40,000 | 0.03198403 ± 0.00006222933 |

`top_campaign.py` runs the production-photon integrations. `top_radiative.py` writes a private copy of the SM restriction with physical lepton masses and runs the four decay integrations. `top_shower.py` reduces disjoint LHE intervals. Pass `--radiative`, the photon-history plugin, `--photon-matcher response_lhe.py` and the radiative response executable for the resolved-decay samples. An independent random seed is used for each hard integration, shower block and detector realization. The first block of each integration keeps its full ROOT output. Others keep the shell output, cards, photon-matching sidecar if used, detector packet, and hashes after validation.

Example commands for the hard-event stage run from this directory with the Herwig module, LHAPDF grids, and Python environment loaded. `WORK` is the absolute path of the top campaign directory and `MG5` is the MadGraph executable.

```bash
python top_campaign.py --out "$WORK" --mg5 "$MG5" --cores 36 --events 30000 30000 40000
python top_campaign.py --out "$WORK" --mg5 "$MG5" --cores 8 --events 300000 120000 1 --round 1 --samples tt_0a tt_1a
python top_radiative.py --out "$WORK" --generator "$PWD/ew_generate.py" --mg5 "$MG5" --cores 36 --events 10000
python top_radiative.py --out "$WORK" --generator "$PWD/ew_generate.py" --mg5 "$MG5" --cores 6 --events 50000 --samples tt_0a_1dec_ll
```

Run `top_diagnostics.py "$WORK" top-diagnostics.json` after reducing the showers. The summary counts eligible electron-fake pairs with a distinct lepton tag over all masses and in the 125 GeV window. For a zero-event window it reports one-sided 95% binomial bounds conditional on no more than one, three or six eligible pairs per collision, as well as the conservative bound for arbitrary multiplicities. It also bins the selected radiative cross section in the minimum hard photon–quark separation. The 0.05–0.10 bin indicates the contribution closest to the generation cut.
