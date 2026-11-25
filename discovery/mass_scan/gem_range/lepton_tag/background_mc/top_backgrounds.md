# Top backgrounds in the tagged diphoton search

GEM tabulated separate top-pair diphoton and reducible top backgrounds for its lepton-tagged search. The updated calculation generates proton–proton collisions at 40 TeV with MadGraph5_aMC@NLO 3.5.15. The matrix elements are tree level. The hard process and the Herwig 7.3 shower both use the NNPDF4.0 NLO QED PDF with LHAPDF ID 335900. The proton includes the gluon and the five massless quark flavours and their antiquarks. The calculation includes the leading QCD top production amplitudes. [GEM Higgs study, Section 5 and Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

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

For a radiating top, the decay amplitude is `t > b f f a` or `t > b f f a a`, and similarly for the antitop. Charge conservation and the SM couplings select the two light quark pairs, electron, muon and tau channels of the W. Both charges are included. This factorizes production and decay at the top pole, while keeping the contributing charged-current diagrams and their interference within each radiative top decay. The hard-event card uses a top mass and width of 173.0 and 1.4915 GeV. The unstable W is internal to the decay amplitude and has a mass and width of 80.419 and 2.0476 GeV. Interference between production and decay or between competing assignments of photons to two resonant tops is beyond this pole approximation. Any other on-shell top decays inclusively in Herwig.

The radiative samples use physical charged-lepton masses and allow zero photon–lepton separation. A muon near a photon can pass the calorimeter isolation if a different lepton tags the event. The default hard cuts require photon transverse momentum above 14 GeV and absolute pseudorapidity below 3.5. The two-matrix-element-photon samples require a hard diphoton mass between 60 and 170 GeV to cover the reconstructed 80–140 GeV scan. In the radiative top samples the quark–photon separation is above 0.05. This regulates the collinear massless-quark singularity. Dependence on this cut, quark-to-photon fragmentation, and higher-order QCD and electroweak effects are model uncertainties. The hard cuts and the detector selection are recorded separately.

For each independent integration the event rate in fb is the LHE cross section in pb times 1,000 times the summed detector weight divided by the signed number of hard events attempted. The LHE cross sections for the radiative decay processes already include the corresponding matrix-element decay rates, with the top and W widths recorded in `param_card.dat`. The electron-to-photon transfer is 0.90/33. The showered and smeared object used for this transfer cannot simultaneously supply the lepton tag. The inclusive fake rates and the jet-flavour variation are implemented in the shared detector response.

The initial integrations gave the following cross sections before detector selection. The two-photon samples include the hard diphoton mass cut. All samples with an explicit top decay include its branching fraction in the listed cross section. The production-photon integrations used seeds 721000–721002 and the radiative-decay integrations used seeds 751000–751004.

| Sample | Hard events | Cross section [pb] |
|---|---:|---:|
| `tt_0a` | 30,000 | 4083.590 ± 3.033376 |
| `tt_1a` | 30,000 | 16.77874 ± 0.03132324 |
| `tt_2a` | 40,000 | 0.03198403 ± 0.00006222933 |
| `tt_0a_1dec` | 10,000 | 106.9852 ± 0.2819 |
| `tt_0a_1dec_ll` | 50,000 | 11.1646 ± 0.0148 |
| `tt_1a_1dec` | 10,000 | 0.216321 ± 0.000616 |
| `tt_0a_2dec_same` | 10,000 | 0.047715 ± 0.000131 |
| `tt_0a_2dec_split` | 10,000 | 0.363811 ± 0.000738 |

## Top contribution to the discovery fit

The frozen discovery fit uses the following selected cross sections for one experiment in the reconstructed mass interval 123.3–126.7 GeV. The uncertainties include finite hard-event and detector sampling and the hard integration errors. The total is **0.24642 ± 0.02117 fb**. Linear interpolation of the two top rows in GEM Table 19 gives **0.1445 fb** over the same interval, about 59% of the present central value. Differences can come from the hard production and decay models, parton distributions, lepton isolation and photon misidentification probabilities. These model differences are separate from the quoted Monte Carlo uncertainties. [GEM Table 19](https://lss.fnal.gov/archive/other/ssc/ssc-gem-tn-93-373.pdf#page=61)

| Source | Hard collisions in the fit | Selected rate [fb] |
|---|---:|---:|
| `tt_0a`, two fakes | 1,015,000 | 0.01523 ± 0.00063 |
| `tt_1a`, production photon and one fake | 150,000 | 0.05716 ± 0.00924 |
| `tt_2a`, two production photons | 40,000 | 0.02862 ± 0.00380 |
| `tt_0a_1dec`, decay photon and inclusive fake component | 10,000 | 0.02229 ± 0.00789 |
| `tt_0a_1dec_ll`, direct dileptonic electron-fake component | 65,500 | 0.03706 ± 0.01119 |
| `tt_1a_1dec`, production and decay photons | 93,000 | 0.06350 ± 0.00970 |
| `tt_0a_2dec_same`, two photons from one top decay | 40,000 | 0.00292 ± 0.00146 |
| `tt_0a_2dec_split`, one photon from each top decay | 70,000 | 0.01962 ± 0.00802 |

Six independent collisions contributed in the opposite-top double-radiation window. The error of 0.00802 fb meets the source goal of 0.011 fb. The source includes 60,000 fresh collisions from four 15,000-event native grid runs. Their weighted-mean preselection cross section was 0.36589 ± 0.00041 pb, 2.46 combined standard deviations from the parent value. The new and parent distributions had two-sample Kolmogorov–Smirnov p-values of 0.055, 0.179 and 0.362 for partonic centre-of-mass energy, leading-photon transverse momentum and diphoton mass. The parent integration determines the normalization for all 70,000 collisions.

The audit checked 1,641 top detector packets against the discovery spectrum. All 70,000 opposite-top hard collisions had one matrix-element photon descending from each top. In the 40,000 same-top hard collisions, both photons descended from the top in 19,903 events and from the antitop in 20,097. Across the selected radiative sources, the Higgs-window contribution with a hard photon within ΔR = 0.05 of a final-state hard muon was 0.01496 ± 0.00608 fb. The physical muon mass regulates this contribution. Errors here are from finite event sampling. [Frozen rates, checksums, angular diagnostics and independent controls](top_diagnostics.json)

`top_campaign.py` runs the production-photon integrations. `top_radiative.py` writes a private copy of the SM restriction with physical lepton masses and runs the four decay integrations. `top_shower.py` reduces disjoint LHE intervals. Pass `--radiative`, the photon-history plugin, `--photon-matcher response_lhe.py` and the radiative response executable for the resolved-decay samples. An independent random seed is used for each hard integration, shower block and detector realization. The first block of each integration keeps its full ROOT output. Others keep the shell output, cards, photon-matching sidecar if used, detector packet, and hashes after validation.

Example commands for the hard-event stage run from this directory with the Herwig module, LHAPDF grids, and Python environment loaded. `WORK` is the absolute path of the top campaign directory and `MG5` is the MadGraph executable.

```bash
python top_campaign.py --out "$WORK" --mg5 "$MG5" --cores 36 --events 30000 30000 40000
python top_campaign.py --out "$WORK" --mg5 "$MG5" --cores 8 --events 300000 120000 1 --round 1 --samples tt_0a tt_1a
python top_radiative.py --out "$WORK" --generator "$PWD/ew_generate.py" --mg5 "$MG5" --cores 36 --events 10000
python top_radiative.py --out "$WORK" --generator "$PWD/ew_generate.py" --mg5 "$MG5" --cores 6 --events 50000 --samples tt_0a_1dec_ll
```

Run `top_diagnostics.py "$WORK" top-diagnostics.json` after reducing the showers. The summary counts eligible electron-fake pairs with a distinct lepton tag over all masses and in the 125 GeV window. For a zero-event window it reports one-sided 95% binomial bounds conditional on no more than one, three or six eligible pairs per collision, as well as the conservative bound for arbitrary multiplicities. It also bins the selected radiative cross section in the minimum hard photon–quark separation. The 0.05–0.10 bin indicates the contribution closest to the generation cut.

Use `--frozen-spectrum /path/to/background.json` to match an archived fit. This verifies the SHA-256 checksums of the source receipts and detector packets named in the published spectrum. The angular and electron-fake diagnostics then include exactly the hard-event intervals used by the fit. The separately generated control directories can still be supplied via `--ee-control` and `--quark-collar`.

The 1,015,000 inclusive top-pair collisions used in the discovery fit produced 1,014,895 showered events. Two had one eligible dielectron fake pair and a separate tag somewhere in the output mass range. Neither pair was in the 125 GeV window. The corresponding 95% limits on the Higgs-window contribution are 0.00886, 0.02658 and 0.05316 fb under respective caps of one, three and six possible unordered electron pairs per collision. The possible pair multiplicity matters when interpreting an empirical zero. The independently generated direct dielectron control below tests a specified pair without assuming a cap on other pairs in the event.

After the fit was frozen, the inclusive top-pair control reached its planned exposure of 1,230,000 hard collisions and 1,229,874 saved shower events. Three had an eligible dielectron fake pair and a distinct tag across the output mass range. All three were outside the Higgs window. At this exposure the corresponding 95% upper limits are **0.00731, 0.02193 and 0.04386 fb** for caps of one, three and six eligible pairs. The control combines the original 30,000 collisions and three additional rounds of 300,000, 300,000 and 600,000, using hard seeds 721000, 721100, 721200 and 721300. Use all four rounds as input to `top_diagnostics.py` without a frozen-spectrum filter to reproduce its `tt_dielectron_fake` result. The additional 215,000 collisions enter only this control. [Final inclusive dielectron-fake audit](top_auxiliary_dielectron.json)

The optional `--ee-control` reads an independent sample in which both top decays explicitly contain a direct W-decay electron. There is one specified primary electron pair per collision. Its probability of faking two photons is bounded by `(0.90/33)²`. The observed number of collisions with *any* eligible electron pair and a distinct tag bounds the probability of an eligible primary pair. The control generated 15,000 hard collisions (47.035 ± 0.087 pb) and saved 14,995 showered events. One had an eligible pair at any mass and zero had one in the 125 GeV window. The one-sided 95% upper bound on the two primary electrons faking the diphoton pair is 0.00690 fb. Events with alternative W decay modes are covered by the inclusive top sample.

The dielectron control used the following commands, run in a separate work directory to keep its events out of the nominal spectrum. Set `CONTROL` to this directory, `REPO` to the project checkout, and `PLUGIN` and `RESPONSE` to the standard detector plugin and response executable.

```bash
python ew_generate.py --mg5 "$MG5" --process tt_0a --work "$CONTROL/hard" \
  --events 15000 --seed 752941 --cores 8 --pdf 335900 --card top_cards/tt_0a_ee_control.mg5
for start in $(seq 0 500 14500); do
  python ew_shower.py --repo "$REPO" --hard "$CONTROL/hard" \
    --out "$CONTROL/packets/hard-$(printf '%06d' "$start")" --sample tt_0a \
    --start "$start" --events 500 --seed "$((85100000 + start))" \
    --plugin "$PLUGIN" --response "$RESPONSE"
done
```

`top_precision.py` checks each 123.3–126.7 GeV window after the initial showers. The production-plus-decay, same-top and split-top sources have Monte Carlo standard error goals of 0.014, 0.009 and 0.011 fb, respectively. The complementary single-radiation and direct dilepton components have goals of 0.010 and 0.012 fb. For fewer than five contributing collisions, the controller also requires a one-sided 95% binomial upper limit below twice that source's error goal. It generates additional rounds with new hard and detector seeds as needed.

`top_gridpack.py` can continue an expensive phase-space integration with MadGraph's native adaptive grids. It builds a gridpack in a private working copy of the completed process, starts native grid runs with distinct seeds, and checks every event against previous hard-event momenta. The combined LHE file uses a common positive nominal weight. Each native LHE file and its generated cross section are preserved. The detector summary groups the new collisions under the parent integration ID so its cross-section uncertainty enters once. `--fresh-integral` selects a new adaptive integration when varying the generation cuts.

The opposite-top double-radiation gridpack used on Timur was 288 MB compressed and 1.165 GB as an uncompressed tar archive. Each native worker uses its own temporary directory. Thirty-five simultaneous workers expanded about 38 GiB of archive data, plus filesystem metadata. Set `--native-scratch /path/to/private/scratch` to unpack native workers on a faster volume while writing their final LHE files to the campaign directory. This also works across filesystems. A 50-event native test on Timur verified the LHE transfer from a memory-backed volume and cleanup of the temporary directory. The number of workers can be reduced to fit a smaller scratch volume. The final diagnostic checks hard-particle identities and momenta across all rounds of each radiative process and tabulates the signed parent tops of the matrix-element photons.

The separate quark–photon control uses `top_radiative.py --samples tt_0a_1dec --fresh-integral --quark-collar --dr .02 --round 3 --events 10000` in its own work directory. Use the enriched photon-history response to reduce these showers, with the same round and source name. A MadGraph user cut requires at least one photon–quark pair with separation below 0.05. The minimum separations of the 10,000 hard events ranged from 0.020003 to 0.049992. Their hard cross section was 14.374 ± 0.048 pb. Three collisions contributed in the reconstructed Higgs window, giving **0.00113 ± 0.00066 fb** from this additional angular interval. Each included one matrix-element decay photon and one jet fake. The adjacent nominal interval from 0.05 to 0.10 supplied 0.00271 ± 0.00271 fb in the inclusive single-decay sample. Pass the control directory to `top_diagnostics.py --quark-collar` to reproduce the rates and source-event Monte Carlo errors. A prediction down to zero separation would require a consistently defined fragmentation contribution or an infrared-safe photon isolation prescription.
