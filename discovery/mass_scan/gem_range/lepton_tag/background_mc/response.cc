// Build with -I<repository-root> and root-config --cflags --libs.
#define TAGGED_BACKGROUND_RESPONSE
#include "../signal/measure.cc"
#include <map>
#include <sstream>

namespace {
enum class Kind { HardPhoton, DecayPhoton, Electron, Jet };
struct Object {
  FourVector p;
  Kind kind;
  int stable = -1, hard = -1;
  bool bjet = false;
};
struct Tag { int stable; FourVector p; double efficiency; };

struct Ancestry : Graph {
  std::vector<int>* hard = nullptr;
  std::vector<int>* photon_index = nullptr;
  std::vector<double>* photon_p4 = nullptr;
  // Require charge-consistent e/mu lines for the top and antitop separately.
  // A tau/hadron/photon intermediate cannot qualify. An off-shell W in an
  // ME decay need not appear as a separate particle in the event record.
  int direct_lepton_top(int i) const {
    const int sign=pdg->at(object->at(i))<0 ? 1 : -1;
    for(int p : first_parents(object->at(i))) {
      if(pdg->at(p)==sign*6) return sign;
      if(pdg->at(p)!=sign*24) continue;
      for(int t : first_parents(p)) if(pdg->at(t)==sign*6) return sign;
    }
    return 0;
  }
  // A hard outgoing particle and all its Herwig copies have the same PDG
  // identity. A photon with a hard parton ancestor alone is not an ME photon.
  int hard_copy(int i) const {
    const int start = object->at(i), id = pdg->at(start);
    std::set<int> visited;
    std::vector<int> todo{start};
    while (!todo.empty()) {
      const int j = todo.back();
      todo.pop_back();
      if (j < 0 || pdg->at(j) != id || !visited.insert(j).second) continue;
      if (hard->at(j) == 1) return j;
      todo.push_back(previous->at(j));
      for (int k = offsets->at(j); k < offsets->at(j+1); ++k) todo.push_back(parents->at(k));
    }
    return -1;
  }
  // Follow the emitter's line through partons/leptons alone. A photon from a
  // pi0 (or another hadron) in a top jet belongs to the jet fake response.
  bool decay_photon(int i) const {
    const auto direct = first_parents(object->at(i));
    for (int j : direct) if (std::abs(pdg->at(j)) >= 100) return false;
    std::set<int> visited;
    std::vector<int> todo(direct.begin(), direct.end());
    while (!todo.empty()) {
      const int j=todo.back();
      todo.pop_back();
      if (j < 0 || !visited.insert(j).second) continue;
      const int id=std::abs(pdg->at(j));
      if (id==24 || id==6) return true;
      // Stop at a hard particle produced with the top pair, MPI, or hadronization.
      if (id>=100 || hard->at(j)!=0) continue;
      todo.push_back(previous->at(j));
      for (int k=offsets->at(j); k<offsets->at(j+1); ++k) todo.push_back(parents->at(k));
    }
    return false;
  }
};

#include "response_me.hh"

int photon_count(const std::string& sample) {
  static const std::map<std::string,int> counts{
    {"tt_0a",0},{"tt_1a",1},{"tt_2a",2},{"w_lnu_aa",2},{"z_ll_aa",2},{"z_ee_a",1},
    {"tt_0a_1dec",0},{"tt_0a_1dec_ll",0},{"tt_1a_1dec",1},{"tt_0a_2dec_same",0},{"tt_0a_2dec_split",0}};
  return counts.at(sample);
}
double id_probability(const Object& o, bool gluon_jets) {
  switch(o.kind) {
    case Kind::HardPhoton:
    case Kind::DecayPhoton: return .90 * .96;
    case Kind::Electron: return .90 / 33.;
    case Kind::Jet:
      // Fixed 125 GeV working point for a single coherent detector-ID
      // outcome across all tested mass windows. Keep b-tagged jets quark-like.
      return jet_to_photon_fake_rate(125., o.bjet ? 5 : 0, gluon_jets ? 0. : 1., 1.);
  }
  throw std::logic_error("Unknown photon candidate");
}
std::string label(Kind k) {
  switch(k) {
    case Kind::HardPhoton: return "hard";
    case Kind::DecayPhoton: return "decay";
    case Kind::Electron: return "electron";
    case Kind::Jet: return "jet";
  }
  throw std::logic_error("Unknown photon candidate");
}

// Two leading successfully identified photons define one candidate. Every
// higher-ranked object outside the chosen pair must have failed its photon ID.
// An electron's tag and photon-fake outcomes are mutually exclusive.
double event_probability(const std::vector<Object>& objects, const std::vector<Tag>& tags,
                         size_t a, size_t b, bool gluon_jets) {
  double prob=.988*id_probability(objects[a],gluon_jets)*id_probability(objects[b],gluon_jets);
  std::map<int,double> electron_miss;
  for(size_t k=0; k<b; ++k) {
    if(k==a) continue;
    const double miss=1.-id_probability(objects[k],gluon_jets);
    prob*=miss;
    if(objects[k].kind==Kind::Electron) electron_miss[objects[k].stable]=miss;
  }
  double miss_tag=1.;
  for(const auto& t : tags) {
    if(t.stable==objects[a].stable || t.stable==objects[b].stable) continue;
    auto it=electron_miss.find(t.stable);
    miss_tag*=1.-t.efficiency/(it==electron_miss.end() ? 1. : it->second);
  }
  return prob*(1.-miss_tag);
}
}

int main(int argc, char** argv) {
  if(argc!=5 && argc!=6) throw std::runtime_error("Use response sample shower.root packet.csv seed [hard-photons.csv]");
  const std::string sample=argv[1];
  const int n_hard=photon_count(sample);
  const bool top=sample.rfind("tt_",0)==0;
  const bool electron_fake=sample=="z_ee_a";
  const int n_decay=sample.find("2dec")!=std::string::npos ? 2 : int(sample.find("1dec")!=std::string::npos);
  const bool record_dilepton=n_decay==1 && n_hard==0;
  if(bool(n_decay)!=(argc==6)) throw std::runtime_error("Dedicated radiative sample requires an LHE photon sidecar");
  const auto seed=std::stoull(argv[4]);
  TFile input(argv[2],"READ");
  auto* tree=dynamic_cast<TTree*>(input.Get("Data"));
  if(!tree) throw std::runtime_error("Missing Data tree");
  const auto me=n_decay ? load_photons(argv[5],tree->GetEntries(),n_hard,n_decay) : MEEvents{};
  tree->SetBranchStatus("*",0);
  static double particles[8][10000],jets[5][100],bjets[5][100],photons[4][100];
  int n=0, nj=0, nb=0, ng=0, schema=0;
  double weight=0.;
  Ancestry graph;
  read(tree,"evweight",&weight);
  read(tree,"numJets",&nj);read(tree,"theJets",jets);
  read(tree,"numbJets",&nb);read(tree,"thebJets",bjets);
  read(tree,"numPhotons",&ng);read(tree,"thePhotons",photons);
  read(tree,"numparticles",&n);read(tree,"objects",particles);
  read(tree,"truthSchemaVersion",&schema);read(tree,"truthHard",&graph.hard);
  read(tree,"truthObject",&graph.object);read(tree,"truthPdg",&graph.pdg);
  read(tree,"truthParents",&graph.parents);read(tree,"truthParentOffsets",&graph.offsets);
  read(tree,"truthPrevious",&graph.previous);
  if(n_decay) {
    read(tree,"truthPhotonIndex",&graph.photon_index);
    read(tree,"truthPhotonP4",&graph.photon_p4);
  }
  std::ofstream out(argv[3]);
  out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(16)<<"source_event,evweight,component,mass_gev,probability,probability_gluon_jets";
  if(record_dilepton) out<<",direct_top_dilepton";
  out<<'\n';
  int matched_all=0, n_dilepton=0;
  for(Long64_t e=0;e<tree->GetEntries();++e) {
    if(tree->GetEntry(e)<=0 || schema!=1 || n<0 || n>10000 || n!=int(graph.object->size()) ||
       !std::isfinite(weight)) throw std::runtime_error("Invalid shower event or ancestry schema");
    std::vector<RawParticle> stable;
    stable.reserve(n);
    for(int i=0;i<n;++i) {
      auto p=FourVector{particles[1][i],particles[2][i],particles[3][i],particles[0][i]};
      stable.push_back({p,std::abs(int(particles[4][i])),p.pt(),p.eta(),p.phi()});
    }
    std::set<int> hard_ids;
    const auto matched=n_decay ? match_photons(stable,graph,me[e]) : std::map<int,int>{};
    std::vector<Object> candidates;
    std::vector<Tag> tags;
    std::vector<FourVector> jet_cleaning_objects;
    bool top_plus=false, top_minus=false;
    for(int i=0;i<n;++i) {
      const auto& s=stable[i];
      if(record_dilepton && (s.id==11 || s.id==13)) {
        const int sign=graph.direct_lepton_top(i);
        if(sign==1) top_plus=true;
        if(sign==-1) top_minus=true;
      }
      if(s.id==22) {
        if(n_decay && !matched.count(i)) continue;
        const int hard=n_decay ? (me[e][matched.at(i)].decay ? -1 : matched.at(i)) : graph.hard_copy(i);
        if(hard>=0) hard_ids.insert(hard);
        if(s.pt<10. || (!n_decay && hard<0 && (!top || !graph.decay_photon(i)))) continue;
        auto smear_seed=object_seed(seed,e,90U,i);
        // Same detector realization if this photon is also in HwSim's photon table.
        for(int k=0;k<std::min(ng,100);++k) {
          const auto g=array_four_vector(photons,k);
          if(delta_r(s.p,g)<1.e-7 && std::abs(s.p.e-g.e)<1.e-6*std::max(1.,g.e)) {
            smear_seed=object_seed(seed,e,1U,k);break;
          }
        }
        const auto p=smear_em(s.p,true,smear_seed);
        if(in_em_fiducial(p,2.5) && isolated(stable,i,.45)) {
          candidates.push_back({p,hard<0 ? Kind::DecayPhoton : Kind::HardPhoton,i,hard});
          jet_cleaning_objects.push_back(p);
        }
      } else if((s.id==11 || s.id==13) && s.pt>10.) {
        const bool electron=s.id==11;
        const auto p=electron ? smear_em(s.p,true,object_seed(seed,e,91U,i)) : s.p;
        const bool fid=electron ? in_em_fiducial(p,2.5) : (in_muon_fiducial(p) && p.pt()>20.);
        if(!fid) continue;
        if(electron_fake && (!electron || (graph.hard_copy(i)<0 &&
            !graph.parent_is(graph.first_parents(graph.object->at(i)),23)))) continue;
        if(isolated(stable,i,.30)) {
          tags.push_back({i,p,electron ? kDefaultElectronEfficiency : kDefaultMuonEfficiency});
          jet_cleaning_objects.push_back(p);
        }
        if((top || electron_fake) && electron && isolated(stable,i,.45))
          candidates.push_back({p,Kind::Electron,i});
      }
    }
    const bool dilepton=top_plus && top_minus;
    if(dilepton) ++n_dilepton;
    if(n_decay ? int(matched.size())==n_hard+n_decay : int(hard_ids.size())==n_hard) ++matched_all;
    if(top && n_hard<2 && !tags.empty()) {
      for(int c=0;c<2;++c) {
        const auto arr=c ? bjets : jets;
        for(int i=0;i<std::min(100,c?nb:nj);++i) {
          auto orig=array_four_vector(arr,i);
          if(!(orig.e>0.)) continue;
          const auto p=smear_em(orig,true,object_seed(seed,e,30U+c,i));
          if(!in_em_fiducial(p,2.5)) continue;
          bool overlaps=false;
          for(const auto& other:jet_cleaning_objects) if(delta_r(other,p)<.45) overlaps=true;
          if(!overlaps) candidates.push_back({p,Kind::Jet,-1,-1,c==1});
        }
      }
    }
    std::stable_sort(candidates.begin(),candidates.end(),[](const Object& a,const Object& b){
      return a.p.pt()>b.p.pt();
    });
    bool any=false;
    for(size_t a=0;a<candidates.size() && !tags.empty();++a) for(size_t b=a+1;b<candidates.size();++b) {
      const auto& one=candidates[a]; const auto& two=candidates[b];
      const int hard=(one.kind==Kind::HardPhoton)+(two.kind==Kind::HardPhoton);
      if(hard!=n_hard || (hard==2 && one.hard==two.hard)) continue;
      const int decay=(one.kind==Kind::DecayPhoton)+(two.kind==Kind::DecayPhoton);
      if(n_decay && decay!=n_decay) continue;
      if(electron_fake && one.kind!=Kind::Electron && two.kind!=Kind::Electron) continue;
      const auto pair=one.p+two.p;
      if(pair.pt()<=40. || ((one.kind==Kind::Jet || two.kind==Kind::Jet) && delta_r(one.p,two.p)<.45)) continue;
      const double prob=event_probability(candidates,tags,a,b,false);
      if(prob<=0.) continue;
      std::string first=label(one.kind),second=label(two.kind);
      if(first>second) std::swap(first,second);
      std::string component=first+"_"+second;
      if(record_dilepton && dilepton && component=="decay_electron") component+="_ll";
      out<<e<<','<<weight<<','<<component<<','<<pair.mass()<<','<<prob<<','
         <<event_probability(candidates,tags,a,b,true);
      if(record_dilepton) out<<','<<dilepton;
      out<<'\n';
      any=true;
    }
    if(!any) {
      out<<e<<','<<weight<<",none,0,0,0";
      if(record_dilepton) out<<','<<dilepton;
      out<<'\n';
    }
  }
  std::cerr<<"Entries with all expected hard photons matched: "<<matched_all<<"/"<<tree->GetEntries()<<"\n";
  // A systematic loss usually indicates that primarySubProcess exposed only
  // an intermediate resonance. Never interpret such a shower as a zero rate.
  if(n_hard+n_decay>0 && tree->GetEntries()>=20 && matched_all<.80*tree->GetEntries())
    throw std::runtime_error("Less than 80% of shower events have the expected hard photon ancestry");
  if(record_dilepton) std::cerr<<"Direct e/mu decays of both tops: "<<n_dilepton<<"/"<<tree->GetEntries()<<"\n";
  if(sample=="tt_0a_1dec_ll" && n_dilepton!=tree->GetEntries())
    throw std::runtime_error("LHE dilepton sample disagrees with signed top truth classification");
  return 0;
}
