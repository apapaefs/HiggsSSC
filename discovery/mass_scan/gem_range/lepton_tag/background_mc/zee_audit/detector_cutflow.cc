// Include the exact production identification, smearing, isolation and ancestry helpers.
#define TAGGED_BACKGROUND_HELPERS_ONLY
#include "../response.cc"

namespace {
struct Em {
  int stable;
  int pdg;
  FourVector reco, truth;
  bool iso_photon, iso_tag, loose_photon, loose_tag;
};

bool basic(const Em& o, bool smear) {
  const auto& p=smear ? o.reco : o.truth;
  // The generation and post-shower pT floors used by the original response
  // differ for electrons and photons. They matter only before analysis cuts.
  const bool input=o.pdg==22 ? o.truth.pt()>=10. : o.truth.pt()>10.;
  return input && p.pt()>=20. && std::abs(p.eta())<2.5;
}

// The 1992 report used a perfect EMC geometry and a relative term in the
// isolation threshold. Use this as a selection variation only; its calorimeter
// shower-shape discriminator cannot be reconstructed from HwSim particles.
bool loose_iso(const std::vector<RawParticle>& particles, int idx, double radius) {
  const auto& a=particles[idx];
  const double limit=5.+.10*a.pt;
  double sum=0.;
  for(size_t i=0;i<particles.size();++i) {
    if(i==size_t(idx)) continue;
    const auto& b=particles[i];
    if(b.id==12 || b.id==13 || b.id==14 || b.id==16 || b.pt<=0.) continue;
    const double deta=a.eta-b.eta;
    if(std::abs(deta)>=radius) continue;
    const double dphi=std::remainder(a.phi-b.phi,2*M_PI);
    if(deta*deta+dphi*dphi<radius*radius) sum+=b.pt;
    if(sum>=limit) return false;
  }
  return true;
}

void emit(std::ostream& out, Long64_t event, const std::vector<Em>& all, int mode) {
  const bool smear=mode==1, precursor=mode==2;
  auto photon_iso=[precursor](const Em& o){return precursor ? o.loose_photon : o.iso_photon;};
  auto tag_iso=[precursor](const Em& o){return precursor ? o.loose_tag : o.iso_tag;};
  auto p4=[smear](const Em& o) -> const FourVector& {return smear ? o.reco : o.truth;};
  auto fiducial=[&](const Em& o){return basic(o,smear) && (precursor || in_em_fiducial(p4(o),2.5));};
  std::vector<Object> candidates;
  std::vector<Tag> tags;
  for(const auto& o:all) {
    if(!fiducial(o)) continue;
    if(photon_iso(o)) candidates.push_back({p4(o),o.pdg==22 ? Kind::HardPhoton : Kind::Electron,o.stable});
    if(o.pdg!=22 && tag_iso(o)) tags.push_back({o.stable,p4(o),kDefaultElectronEfficiency});
  }
  std::stable_sort(candidates.begin(),candidates.end(), [](const Object& a,const Object& b){return a.p.pt()>b.p.pt();});
  std::map<int,size_t> rank;
  for(size_t i=0;i<candidates.size();++i) rank.emplace(candidates[i].stable,i);
  for(const auto& g:all) {
    if(g.pdg!=22) continue;
    for(const auto& e:all) {
      if(e.pdg==22) continue;
      const auto pair=p4(g)+p4(e);
      const double m=pair.mass();
      if(m<75. || m>=165.) continue;
      int raw=0, kin=0, fid=0, iso=0;
      for(const auto& other:all) {
        if(other.pdg==22 || other.stable==e.stable) continue;
        ++raw;
        if(!basic(other,smear)) continue;
        ++kin;
        if(!fiducial(other)) continue;
        ++fid;
        if(tag_iso(other)) ++iso;
      }
      if(!raw) continue;
      // Each row is a possible (real photon, fake electron), with a separate
      // identified electron for the tag. Count the row once even if several tags exist.
      const bool cuts[]={basic(g,smear), basic(e,smear), kin>0,
                         fiducial(g)&&fiducial(e)&&fid>0,
                         photon_iso(g), photon_iso(e), iso>0, pair.pt()>40.};
      int ncut=0;
      for(bool pass:cuts) {if(!pass) break; ++ncut;}
      double probability=0.;
      if(ncut>=7) {
        const auto a=rank.at(g.stable),b=rank.at(e.stable);
        // Factor at the fully identified/triggered step, with or without the
        // pair-pT requirement. Includes the mutual exclusion of tag and fake.
        probability=event_probability(candidates,tags,std::min(a,b),std::max(a,b),false);
      }
      out<<event<<','<<mode<<','<<m<<','<<ncut<<','<<e.pdg<<','<<iso<<','<<probability<<'\n';
    }
  }
}
}

int main(int argc,char** argv) {
  if(argc!=4) throw std::runtime_error("Use detector_cutflow shower.root output.csv response-seed");
  const auto seed=std::stoull(argv[3]);
  TFile input(argv[1],"READ");
  auto* tree=dynamic_cast<TTree*>(input.Get("Data"));
  if(!tree) throw std::runtime_error("Missing Data tree");
  tree->SetBranchStatus("*",0);
  static double particles[8][10000],photons[4][100];
  int n=0,ng=0,schema=0;
  double weight=0.;
  Ancestry graph;
  read(tree,"evweight",&weight);read(tree,"numPhotons",&ng);read(tree,"thePhotons",photons);
  read(tree,"numparticles",&n);read(tree,"objects",particles);read(tree,"truthSchemaVersion",&schema);
  read(tree,"truthHard",&graph.hard);read(tree,"truthObject",&graph.object);read(tree,"truthPdg",&graph.pdg);
  read(tree,"truthParents",&graph.parents);read(tree,"truthParentOffsets",&graph.offsets);
  read(tree,"truthPrevious",&graph.previous);
  std::ofstream out(argv[2]);
  out.exceptions(std::ios::badbit|std::ios::failbit);
  out<<std::setprecision(16)<<"event,mode,mass,passed,pdg_fake,eligible_tags,probability\n";
  for(Long64_t ev=0;ev<tree->GetEntries();++ev) {
    if(tree->GetEntry(ev)<=0 || schema!=1 || n<0 || n>10000 || n!=int(graph.object->size()) ||
       std::abs(weight-1.)>1e-8) throw std::runtime_error("Malformed event or unexpected event weight");
    std::vector<RawParticle> stable;
    for(int i=0;i<n;++i) {
      const FourVector p{particles[1][i],particles[2][i],particles[3][i],particles[0][i]};
      stable.push_back({p,std::abs(int(particles[4][i])),p.pt(),p.eta(),p.phi()});
    }
    std::vector<Em> em;
    for(int i=0;i<n;++i) {
      const auto& s=stable[i];
      if(s.id!=11 && s.id!=22) continue;
      const int hard=graph.hard_copy(i);
      if(hard<0 && (s.id==22 || !graph.parent_is(graph.first_parents(graph.object->at(i)),23))) continue;
      auto smear_seed=object_seed(seed,ev,s.id==22 ? 90U : 91U,i);
      if(s.id==22) {
        for(int k=0;k<std::min(ng,100);++k) {
          const auto g=array_four_vector(photons,k);
          if(delta_r(s.p,g)<1.e-7 && std::abs(s.p.e-g.e)<1.e-6*std::max(1.,g.e)) {
            smear_seed=object_seed(seed,ev,1U,k);break;
          }
        }
      }
      em.push_back({i,graph.pdg->at(graph.object->at(i)),smear_em(s.p,true,smear_seed),s.p,
                    isolated(stable,i,.45),s.id==11 && isolated(stable,i,.30),
                    loose_iso(stable,i,.45),s.id==11 && loose_iso(stable,i,.30)});
    }
    for(int mode=0;mode<3;++mode) emit(out,ev,em,mode);
  }
  std::cerr<<tree->GetEntries()<<" saved collisions\n";
}
