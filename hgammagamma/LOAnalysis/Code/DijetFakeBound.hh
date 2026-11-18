// Included inside the SSC response implementation so the two analyses use
// identical smearing, acceptance, hard-parton matching and GEM probabilities.
int write_dijet_bound(int argc, char* argv[]) {
  if (argc != 5 && argc != 6)
    throw std::runtime_error("--dijet-bound input.root moments.csv seed [pair-features.csv]");
  std::ofstream pairs;
  if (argc == 6) {
    pairs.open(argv[5]);
    if (!pairs) throw std::runtime_error("could not open pair-feature output");
    pairs << "source,eligible_jets,f1,f2,fq,fg,m_gg,pt_gamma1,eta_gamma1,pt_gamma2,eta_gamma2,deltaR_gg,deltaPhi_gg,pt_gg,y_gg\n";
    pairs << std::setprecision(17);
  }
  TFile input(argv[2], "READ");
  auto* tree = dynamic_cast<TTree*>(input.Get("Data"));
  if (!tree || tree->GetEntries() == 0) throw std::runtime_error("missing HwSim Data entries");
  const unsigned long seed = std::stoul(argv[4]);
  double jets[5][kMaxRecoObjects]{}, bjets[5][kMaxRecoObjects]{}, partons[5][kMaxRecoObjects]{};
  int nj = 0, nb = 0, np = 0;
  double weight = 0;
  tree->SetBranchStatus("*", 0);
  auto bind = [&](const char* name, void* address) {
    if (!tree->GetBranch(name)) throw std::runtime_error(std::string("missing branch ") + name);
    tree->SetBranchStatus(name, 1);
    if (tree->SetBranchAddress(name, address) < 0) throw std::runtime_error(name);
  };
  bind("numJets", &nj); bind("numbJets", &nb); bind("numoutgoing", &np);
  bind("theJets", jets); bind("thebJets", bjets); bind("partons", partons); bind("evweight", &weight);
  const std::array<double, 3> scales = {{0.5, 1., 2.}};
  const std::array<const char*, 5> scenarios = {{"two_eligible_only", "independent_pair_sum",
      "all_quark_independent", "all_quark_any_correlation", "all_gluon_two_eligible"}};
  double sum[scales.size()][scenarios.size()]{}, sum2[scales.size()][scenarios.size()]{};
  long long positive[scales.size()][scenarios.size()]{};
  for (Long64_t event = 0; event < tree->GetEntries(); ++event) {
    if (tree->GetEntry(event) <= 0 || std::abs(weight - 1.) > 1.e-8 ||
        std::min({nj,nb,np}) < 0 || std::max({nj,nb,np}) > kMaxRecoObjects)
      throw std::runtime_error("invalid HwSim entry or non-unit dijet weight");
    std::vector<Candidate> eligible;
    std::vector<int> flavours;
    for (int source = 0; source < 2; ++source) {
      for (int j = 0; j < (source ? nb : nj); ++j) {
        Candidate jet;
        jet.original = array_four_vector(source ? bjets : jets, j);
        if (!(jet.original.e > 0) || !std::isfinite(jet.original.e)) throw std::runtime_error("invalid jet");
        jet.reconstructed = smear_em(jet.original, true, object_seed(seed, event, 30U + source, j));
        if (in_em_fiducial(jet.reconstructed, kPhotonAnalysisEtaMax)) {
          eligible.push_back(jet);
          flavours.push_back(match_parton_id(jet, np, partons));
        }
      }
    }
    double bounds[scales.size()][scenarios.size()]{};
    for (std::size_t i = 0; i < eligible.size(); ++i) {
      for (std::size_t j = i + 1; j < eligible.size(); ++j) {
        const double mass = (eligible[i].reconstructed + eligible[j].reconstructed).mass();
        if (argc == 6) {
          SelectedHypothesis pair;
          pair.first = eligible[i].reconstructed;
          pair.second = eligible[j].reconstructed;
          order_selected_pair(pair);
          const auto diphoton = pair.first + pair.second;
          pairs << event << ',' << eligible.size() << ','
                << jet_to_photon_fake_rate(mass, flavours[i], 1., 1.) << ','
                << jet_to_photon_fake_rate(mass, flavours[j], 1., 1.) << ','
                << jet_to_photon_fake_rate(mass, 1, 1., 1.) << ','
                << jet_to_photon_fake_rate(mass, 21, 1., 1.) << ','
                << mass << ',' << pair.first.pt() << ',' << pair.first.eta() << ','
                << pair.second.pt() << ',' << pair.second.eta() << ','
                << delta_r(pair.first, pair.second) << ',' << delta_phi(pair.first, pair.second) << ','
                << diphoton.pt() << ',' << diphoton.rapidity() << '\n';
        }
        if (mass < 120. || mass > 130.) continue;
        for (std::size_t k = 0; k < scales.size(); ++k) {
          const double a = jet_to_photon_fake_rate(mass, flavours[i], 1., scales[k]);
          const double b = jet_to_photon_fake_rate(mass, flavours[j], 1., scales[k]);
          const double q = jet_to_photon_fake_rate(mass, 1, 1., scales[k]);
          if (eligible.size() == 2) bounds[k][0] += a*b;
          bounds[k][1] += a*b;
          bounds[k][2] += q*q;
          bounds[k][3] += q; // Fréchet upper bound min(q,q) per pair.
          if (eligible.size() == 2) {
            const double g = jet_to_photon_fake_rate(mass, 21, 1., scales[k]);
            bounds[k][4] += g*g;
          }
        }
      }
    }
    for (std::size_t k = 0; k < scales.size(); ++k)
      for (std::size_t m = 0; m < scenarios.size(); ++m) {
        const double value = std::min(1., bounds[k][m]);
        sum[k][m] += value;
        sum2[k][m] += value * value; // square the total for one source event
        positive[k][m] += value > 0.;
      }
  }
  std::ofstream out(argv[3]);
  out << "fake_scale,scenario,events_read,contributing_events,sum_probability,sum_probability_squared\n";
  out << std::setprecision(17);
  for (std::size_t k = 0; k < scales.size(); ++k)
    for (std::size_t m = 0; m < scenarios.size(); ++m)
      out << scales[k] << ',' << scenarios[m] << ',' << tree->GetEntries() << ',' << positive[k][m]
          << ',' << sum[k][m] << ',' << sum2[k][m] << '\n';
  if (!out) throw std::runtime_error("could not write dijet moments");
  return 0;
}
