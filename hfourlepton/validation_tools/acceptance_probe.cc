// Read-only definition-sensitivity diagnostic, not a replacement selection.
// Include the analyzer recorded by the input manifest, using its directory as
// the compiler's -I argument. The Python driver builds in a fresh directory.
//
// Photon-isolation alternatives are deliberately diagnostic: counting every
// photon is not the CMS ancestry-based FSR prescription. The hard-lepton angular
// matching is also only a proxy; shower recoil makes it fail for many events.
// Its rates must not be presented as ancestry-matched fiducial acceptances.
#define main fourlepton_analyzer_program_main
#include "HwSimPostAnalysis_fourlepton.cc"
#undef main
#include <map>
#include <set>

bool sequential_pairing_pass(const std::vector<Lepton>& leptons, const Config& config) {
  std::vector<int> good;
  for(size_t i=0;i<leptons.size();++i) if(passes_baseline(leptons[i],config)) good.push_back(i);
  if(good.size()<4)return false;
  int a=-1,b=-1,c=-1,d=-1;
  double nearest=1.e100;
  for(int i:good)for(int j:good)if(i<j && is_ossf(leptons[i],leptons[j])) {
    double distance=std::abs((leptons[i].reconstructed+leptons[j].reconstructed).mass()-kZMassGeV);
    if(distance<nearest){nearest=distance;a=i;b=j;}
  }
  if(a<0)return false;
  const auto z1=leptons[a].reconstructed+leptons[b].reconstructed;
  if(z1.mass()<=40. || z1.mass()>=120.)return false;
  double sumpt=-1.;
  for(int i:good)for(int j:good)if(i<j && i!=a && i!=b && j!=a && j!=b && is_ossf(leptons[i],leptons[j])) {
    const auto z2=leptons[i].reconstructed+leptons[j].reconstructed;
    const double pt=leptons[i].reconstructed.pt()+leptons[j].reconstructed.pt();
    if(z2.mass()>12. && z2.mass()<120. && pt>sumpt){sumpt=pt;c=i;d=j;}
  }
  if(c<0)return false;
  const auto p4=z1+leptons[c].reconstructed+leptons[d].reconstructed;
  return lhc_topology_pass({{a,b,c,d}},leptons) && p4.mass()>105. && p4.mass()<140.;
}

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  TFile file(argv[1], "READ");
  auto* tree = dynamic_cast<TTree*>(file.Get("Data"));
  if (!tree) return 3;
  static double objects[8][10000], partons[5][100];
  int count = 0, outgoing = 0;
  double weight = 0.;
  tree->SetBranchStatus("*", 0);
  for (const char* name : {"numparticles", "objects", "numoutgoing", "partons", "evweight"}) tree->SetBranchStatus(name, 1);
  tree->SetBranchAddress("numparticles", &count);
  tree->SetBranchAddress("objects", objects);
  tree->SetBranchAddress("numoutgoing", &outgoing);
  tree->SetBranchAddress("partons", partons);
  tree->SetBranchAddress("evweight", &weight);
  Config config;
  config.selection_profile = "cms";
  std::map<std::string, WeightSums> sums;
  Diagnostics diagnostics;
  auto evaluate = [&](const std::vector<Lepton>& leptons) {
    std::vector<int> subset(leptons.size());
    std::iota(subset.begin(), subset.end(), 0);
    return analyze_subset(leptons, subset, config);
  };
  auto passed = [](const SubsetResult& r) {
    return (r.cut_mask & kSelectedMask) == kSelectedMask && r.four_lepton.mass() > 105. && r.four_lepton.mass() < 140.;
  };
  for (Long64_t e = 0; e < tree->GetEntries(); ++e) {
    tree->GetEntry(e);
    if (count < 0 || count > 10000 || outgoing < 0 || outgoing > 100) return 4;
    std::vector<Particle> particles, hard;
    for (int i=0; i<count; ++i) particles.push_back({from_array(objects,i), int(std::lround(objects[4][i])),i});
    for (int i=0; i<outgoing; ++i) {
      const int id=int(std::lround(partons[4][i]));
      if (std::abs(id)==11 || std::abs(id)==13) hard.push_back({{partons[1][i],partons[2][i],partons[3][i],partons[0][i]},id,i});
    }
    sums["all"].add(weight);
    auto leptons=response_leptons(particles, config,e,diagnostics);
    auto chosen=evaluate(leptons);
    bool nominal=passed(chosen);
    if(sequential_pairing_pass(leptons,config))sums["sequential_pairing"].add(weight);
    if (nominal) {
      sums["nominal"].add(weight);
      if (chosen.four_lepton.mass()>120. && chosen.four_lepton.mass()<130.) sums["nominal_120_130"].add(weight);
      int accepted=std::count_if(leptons.begin(),leptons.end(),[&](const Lepton& l){return passes_baseline(l,config);});
      sums[accepted>4 ? "nominal_extra_accepted" : "nominal_four_accepted"].add(weight);
    }
    auto own=leptons, full=leptons, bare=leptons;
    for (size_t i=0; i<leptons.size(); ++i) {
      own[i].isolation_et=isolation_et(own[i],particles,config);
      auto temp=full[i]; temp.dressed_photon_indices.clear();
      full[i].isolation_et=isolation_et(temp,particles,config);
      bare[i].reconstructed=bare[i].bare;
      bare[i].dressed_photon_indices.clear();
      bare[i].isolation_et=isolation_et(bare[i],particles,config);
    }
    for (const auto& item : {std::pair<const char*,const std::vector<Lepton>*>{"only_own_photons_removed_from_iso",&own},
                            {"all_photons_in_iso",&full},{"bare",&bare}}) {
      const bool p=passed(evaluate(*item.second));
      if (p) sums[item.first].add(weight);
      if (nominal && !p) sums[std::string(item.first)+"_lost"].add(weight);
      if (!nominal && p) sums[std::string(item.first)+"_gained"].add(weight);
      if(sequential_pairing_pass(*item.second,config))sums[std::string(item.first)+"_sequential"].add(weight);
    }
    if (hard.size()!=4) {sums["hard_not_four"].add(weight);continue;}
    FourVector hard_p4;
    for(const auto& p:hard)hard_p4=hard_p4+p.p4;
    const double hard_mass=hard_p4.mass();
    const std::string bin=hard_mass<100.?"below100":hard_mass<105.?"100_105":hard_mass<115.?"105_115":hard_mass<135.?"115_135":hard_mass<140.?"135_140":hard_mass<145.?"140_145":hard_mass<160.?"145_160":"above160";
    sums["hard_mass_"+bin].add(weight);
    if(nominal)sums["selected_hard_mass_"+bin].add(weight);
    if(nominal && chosen.four_lepton.mass()>120. && chosen.four_lepton.mass()<130.)sums["window_hard_mass_"+bin].add(weight);
    auto born=response_leptons(hard,config,e,diagnostics);
    const auto born_result=evaluate(born);
    if (passed(born_result)) sums["born_no_isolation"].add(weight);
    std::vector<int> match;
    std::set<int> used;
    bool matched=true;
    for (const auto& h:hard) {
      int best=-1; double dr=.1;
      for(size_t i=0;i<leptons.size();++i) if(leptons[i].pdg_id==h.pdg_id && delta_r(leptons[i].bare,h.p4)<dr) {
        dr=delta_r(leptons[i].bare,h.p4);best=i;
      }
      if(best<0 || used.count(best)) {matched=false;break;}
      used.insert(best);match.push_back(best);
    }
    if(!matched) {sums["hard_match_failed"].add(weight);if(nominal)sums["nominal_hard_match_failed"].add(weight);continue;}
    sums["hard_match_complete"].add(weight);
    std::vector<Lepton> hleptons;
    for(int i:match)hleptons.push_back(leptons[i]);
    auto matched_result=evaluate(hleptons);
    if(passed(matched_result))sums["matched_four_only"].add(weight);
    if(nominal) {
      bool all_h=true;
      for(int i:chosen.output_indices) all_h &= used.count(i)>0;
      sums[all_h?"nominal_all_four_matched":"nominal_includes_extra"].add(weight);
    }
    for(size_t i=0;i<born.size();++i) born[i].isolation_et=hleptons[i].isolation_et;
    if(passed(evaluate(born)))sums["born_with_reco_isolation"].add(weight);
  }
  std::cout << std::setprecision(17) << "{";
  bool first=true;
  for(const auto& item:sums) {
    if(!first)std::cout << ",";first=false;
    const auto& s=item.second;
    std::cout << '"' << item.first << "\":{\"entries\":" << s.entries << ",\"sumw\":" << s.sumw << ",\"sumw2\":" << s.sumw2 << ",\"sumabsw\":" << s.sumabsw << "}";
  }
  std::cout << "}\n";
}
