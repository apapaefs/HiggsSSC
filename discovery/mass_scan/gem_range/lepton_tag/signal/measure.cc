// Apply the inclusive SSC photon response and the additional GEM lepton tag
// to the same showered events. Compile with the repository root on the include path.
#define main inclusive_response_cli
#include "hgammagamma/LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC.cc"
#undef main
#include <set>
#include <iomanip>

namespace {
struct Graph {
  std::vector<int> *object=nullptr, *pdg=nullptr, *parents=nullptr, *offsets=nullptr, *previous=nullptr;
  // Follow a particle's copies backwards. Leave hadronization and shower
  // recoil siblings out of the direct H->gamma and W->lepton ancestry test.
  std::set<int> first_parents(int start) const {
    std::set<int> end, seen;
    std::vector<int> todo{start};
    const int identity=std::abs(pdg->at(start));
    while (!todo.empty()) {
      const int i=todo.back(); todo.pop_back();
      if (i<0 || !seen.insert(i).second) continue;
      const int id=std::abs(pdg->at(i));
      if (id!=identity) { end.insert(i); continue; }
      todo.push_back(previous->at(i));
      for (int j=offsets->at(i);j<offsets->at(i+1);++j) todo.push_back(parents->at(j));
    }
    return end;
  }
  bool parent_is(const std::set<int>& nodes,int id) const {
    for (int p:nodes) if (std::abs(pdg->at(p))==id) return true;
    return false;
  }
  // Exclusive categories in priority order.
  int lepton_origin(int i) const {
    const auto parents1=first_parents(object->at(i));
    if (parent_is(parents1,24)) return 0;
    for (int p:parents1) if (std::abs(pdg->at(p))==15 && parent_is(first_parents(p),24)) return 1;
    // Follow ancestry through hadron decay chains to distinguish b/c.
    std::vector<int> todo(parents1.begin(),parents1.end());
    std::set<int> seen;
    while (!todo.empty()) {
      int j=todo.back();todo.pop_back();
      if (j<0 || !seen.insert(j).second) continue;
      int id=std::abs(pdg->at(j));
      if (id<=100 || id>=1000000000) continue;
      for (int pos:{10,100,1000}) if ((id/pos)%10==4 || (id/pos)%10==5) return 2;
      todo.push_back(previous->at(j));
      for(int k=offsets->at(j);k<offsets->at(j+1);++k) todo.push_back(parents->at(k));
    }
    return 3;
  }
};
struct RawParticle { FourVector p; int id; double pt,eta,phi; };
bool isolated(const std::vector<RawParticle>& stable,int index,double radius) {
  double et=0.;
  const auto &candidate=stable.at(index);
  for (size_t j=0;j<stable.size();++j) {
    if (j==(size_t)index) continue;
    const auto &other=stable[j];
    if (other.id==12 || other.id==13 || other.id==14 || other.id==16 || other.pt<=0) continue;
    double deta=other.eta-candidate.eta;
    if (std::abs(deta)>=radius) continue;
    double dphi=std::remainder(other.phi-candidate.phi,2*M_PI);
    if (deta*deta+dphi*dphi<radius*radius) et+=other.pt;
    if (et>=5.) return false;
  }
  return true;
}
template <class T> void read(TTree *t,const char *name,T* address) {
  if(!t->GetBranch(name)) throw std::runtime_error(std::string("Missing ")+name);
  t->SetBranchStatus(name,1);t->SetBranchAddress(name,address);
}
}
#ifndef TAGGED_BACKGROUND_RESPONSE
int main(int argc,char** argv) {
  if(argc!=4) throw std::runtime_error("Use measure shower.root output.csv response-seed");
  TFile input(argv[1],"READ");
  TTree* tree=dynamic_cast<TTree*>(input.Get("Data"));
  if(!tree) throw std::runtime_error("Missing Data tree");
  tree->SetBranchStatus("*",0);
  Config cfg;cfg.seed=std::stoul(argv[3]);
  double jets[5][100],bjets[5][100],photons[4][100];
  static double particles[8][10000];
  int n=0,ng=0,nj=0,nb=0;double event_weight=0;
  Graph graph;
  read(tree,"numJets",&nj);read(tree,"theJets",jets);
  read(tree,"numbJets",&nb);read(tree,"thebJets",bjets);
  read(tree,"numPhotons",&ng);read(tree,"thePhotons",photons);
  read(tree,"numparticles",&n);read(tree,"objects",particles);read(tree,"evweight",&event_weight);
  read(tree,"truthObject",&graph.object);read(tree,"truthPdg",&graph.pdg);
  read(tree,"truthParents",&graph.parents);read(tree,"truthParentOffsets",&graph.offsets);
  read(tree,"truthPrevious",&graph.previous);
  std::ofstream out(argv[2]);out<<std::setprecision(12)
      <<"evweight,mass_inclusive,mass_tag,p_inclusive,p_tag,p_overlap,p_direct,p_tau,p_hf,p_other,n_higgs_photons\n";
  for (Long64_t e=0;e<tree->GetEntries();++e) {
    tree->GetEntry(e);
    if(n>10000 || n!=(int)graph.object->size()) throw std::runtime_error("Corrupt stable object array");
    // Match the existing inclusive implementation including truncation to 100
    // reconstructed objects, its two leading isolated photons and exact seed.
    std::vector<Candidate> isolation_jets, inclusive_photons;
    for(int c=0;c<2;++c) {
      const double (*arr)[100]=c?bjets:jets;
      for(int i=0;i<std::min(100,c?nb:nj);++i) {
        const auto p=array_four_vector(arr,i);
        if(!(p.e>0.)) continue;
        Candidate jet;jet.reconstructed=smear_jet(p,object_seed(cfg.seed,e,20U+c,i));
        if(in_isolation_jet_fiducial(jet.reconstructed)) isolation_jets.push_back(jet);
      }
    }
    std::vector<Candidate> all;
    for(int i=0;i<std::min(100,ng);++i) {
      Candidate g;g.input_index=i;g.original=array_four_vector(photons,i);
      if(!(g.original.e>0.)) continue;
      g.reconstructed=smear_em(g.original,true,object_seed(cfg.seed,e,1U,i));
      all.push_back(g);
      if(in_em_fiducial(g.reconstructed,kPhotonAnalysisEtaMax) && photon_isolated_from_jets(g.reconstructed,isolation_jets))
        inclusive_photons.push_back(g);
    }
    sort_by_pt(inclusive_photons);
    double p_inc=0,m_inc=0,eff_inc=0,m_tag=0,p_tag=0,p_both=0;
    if(inclusive_photons.size()>=2) {
      m_inc=(inclusive_photons[0].reconstructed+inclusive_photons[1].reconstructed).mass();
      eff_inc=photon_efficiency_for_mass(m_inc,cfg);p_inc=eff_inc*eff_inc;
    }
    std::vector<RawParticle> stable;std::vector<int> genuine,leptons;
    for(int i=0;i<n;++i) {
      int id=std::abs(int(particles[4][i]));
      FourVector p{particles[1][i],particles[2][i],particles[3][i],particles[0][i]};
      stable.push_back({p,id,p.pt(),p.eta(),p.phi()});
      if(id==22 && p.pt()>2. && graph.parent_is(graph.first_parents(graph.object->at(i)),25)) genuine.push_back(i);
      if((id==11 || id==13) && p.pt()>10. && std::abs(p.eta())<2.5) leptons.push_back(i);
    }
    double components[4]={};
    if(genuine.size()==2) {
      std::vector<Candidate> tag;
      for(int idx:genuine) {
        auto p=stable[idx].p;
        Candidate g;g.input_index=-idx-1;g.original=p;
        g.reconstructed=smear_em(p,true,object_seed(cfg.seed,e,90U,idx));
        for(auto &match:all)
          if(delta_r(p,match.original)<1.e-7 && std::abs(p.e-match.original.e)<1.e-6*std::max(p.e,1.)) {
            g=match;break;
          }
        if(in_em_fiducial(g.reconstructed,kPhotonAnalysisEtaMax) && isolated(stable,idx,.45)) tag.push_back(g);
      }
      if(tag.size()==2 && (tag[0].reconstructed+tag[1].reconstructed).pt()>40.) {
        m_tag=(tag[0].reconstructed+tag[1].reconstructed).mass();
        double miss[4]={1.,1.,1.,1.};
        for(int i:leptons) {
          bool electron=stable[i].id==11;
          const auto lepton=electron ? smear_em(stable[i].p,true,object_seed(cfg.seed,e,91U,i)) : stable[i].p;
          bool geom=electron ? in_em_fiducial(lepton,2.5) : (in_muon_fiducial(lepton) && lepton.pt()>20.);
          if(geom && isolated(stable,i,.30))
            miss[graph.lepton_origin(i)]*=1-(electron?cfg.electron_efficiency:cfg.muon_efficiency);
        }
        const double eff_tag=photon_efficiency_for_mass(m_tag,cfg);
        double remaining=.988*eff_tag*eff_tag;
        for(int j=0;j<4;++j) {components[j]=remaining*(1-miss[j]);remaining*=miss[j];p_tag+=components[j];}
        if(p_inc && p_tag) {
          // Common detector-ID outcome for shared photons; independent ID
          // for a genuinely different photon hypothesis in the same collision.
          p_both=p_tag;
          for(int i=0;i<2;++i) {
            int idx=inclusive_photons[i].input_index;
            bool shared=idx==tag[0].input_index || idx==tag[1].input_index;
            p_both*=shared?std::min(eff_inc,eff_tag)/eff_tag:eff_inc;
          }
        }
      }
    }
    out<<event_weight<<','<<m_inc<<','<<m_tag<<','<<p_inc<<','<<p_tag<<','<<p_both;
    for(double c:components) out<<','<<c;
    out<<','<<genuine.size()<<'\n';
  }
}
#endif
