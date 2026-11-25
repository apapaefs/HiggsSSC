// Pointlike calorimeter-deposit proxies for the GEM isolation and preselection.
// This counter leaves the original reconstructed momenta, response seed and
// mutually exclusive tag/fake probabilities unchanged.
#define TAGGED_BACKGROUND_HELPERS_ONLY
#include "../../response.cc"
#include <array>
#include <iomanip>

namespace {
constexpr double pitch=.026;
struct GridCell {
  int eta_bin,phi_bin;
  double eta,phi,dphi;
};
GridCell locate(double eta,double phi,bool stepped,bool hadronic=false,double uniform_pitch=pitch) {
  // Approximate the η transitions and relative HC/EM pitch in Fig. 3 of
  // TN-93-376. Round bin counts to tessellate each range and the full azimuth.
  constexpr std::array<double,6> edges{0.,1.56,1.90,2.20,2.66,4.0};
  constexpr std::array<double,5> pitches{.026,.037,.052,.075,.11};
  double angle=uniform_pitch, lo=0., width=uniform_pitch;
  int offset=0;
  if(stepped) {
    size_t band=0;
    for(;band+1<pitches.size() && std::abs(eta)>=edges[band+1];++band) {
      const double target=pitches[band]*(hadronic ? 3. : 1.);
      offset+=std::lround((edges[band+1]-edges[band])/target);
    }
    angle=pitches[band]*(hadronic ? 3. : 1.);
    lo=edges[band];
    const double span=edges[band+1]-lo;
    width=span/std::lround(span/angle);
  }
  const int eb=static_cast<int>(std::floor((std::abs(eta)-lo)/width));
  const double center=lo+(eb+.5)*width;
  const int np=std::lround(2.*M_PI/angle);
  const double dp=2.*M_PI/np;
  const int pb=(static_cast<int>(std::floor((phi+M_PI)/dp))+np)%np;
  return {eta>=0. ? offset+eb : -1-offset-eb,pb,eta>=0. ? center : -center,
          (pb+.5)*dp-M_PI,dp};
}
int cell_distance(const GridCell& cell,const GridCell& seed) {
  const double dphi=std::abs(std::remainder(cell.phi-seed.phi,2.*M_PI));
  return std::max(std::abs(cell.eta_bin-seed.eta_bin),int(std::floor(dphi/seed.dphi+.5)));
}
struct TaggedEM {int index,pdg; FourVector reco;};
struct Cell {double energy=0.; GridCell pos{}; bool em=false;};
// Energy-weighted mean three-dimensional opening angle in GEM TN-92-126,
// Eqs. 8-9, using the original particle directions before showering in matter.
double mean_opening(const std::vector<FourVector>& particles) {
  if(particles.size()<2) return 0.;
  std::array<double,3> axis{};
  double energy=0.;
  for(const auto& p:particles) {
    const double scale=p.e/p.pabs();
    axis[0]+=scale*p.px;axis[1]+=scale*p.py;axis[2]+=scale*p.pz;
    energy+=p.e;
  }
  const double norm=std::hypot(std::hypot(axis[0],axis[1]),axis[2]);
  if(norm<=0. || energy<=0.) throw std::runtime_error("Undefined EM centroid");
  double weighted=0.;
  for(const auto& p:particles) {
    const double dot=(axis[0]*p.px+axis[1]*p.py+axis[2]*p.pz)/(norm*p.pabs());
    weighted+=p.e*std::acos(std::clamp(dot,-1.,1.));
  }
  return weighted/energy;
}
struct Features {
  double mean_em_angle_rad=0.,mean92_em_angle_rad=0.;
  double cone=0.,core_iso=0.,cell_iso_et=0.,cell_iso_e=0.,core_em_et=0.,core3_em_et=0.;
  double core_had_et=0.,core_had_e=0.,extra_track_pt=0.,extra_track3_pt=0.,extra_em_pt=0.,extra_em_dr=-1.;
  int core_merged=0;
};
Features compute(const std::vector<RawParticle>& stable, const std::vector<double>& charge,
                 const TaggedEM& c, double radius,bool stepped) {
  const auto& a=stable[c.index];
  const auto seed=locate(a.eta,a.phi,stepped);
  const double primary=a.p.e/std::cosh(seed.eta);
  const auto seed92=locate(a.eta,a.phi,false,false,.04);
  std::vector<FourVector> em_directions,em92_directions;
  std::map<std::tuple<bool,int,int>,Cell> cells;
  Features f;
  for (size_t i=0;i<stable.size();++i) {
    const auto& o=stable[i];
    if (o.id==12 || o.id==14 || o.id==16 || o.pt<=0.) continue;
    const double det=o.eta-a.eta,dph=std::remainder(o.phi-a.phi,2.*M_PI);
    // Include HC towers whose centers enter the cone from outside it.
    if (std::abs(det)>radius+.2 || std::abs(dph)>radius+.2) continue;
    const double dr=std::hypot(det,dph);
    if (int(i)!=c.index && o.id!=13 && dr<radius) f.cone+=o.pt;
    const auto empos=locate(o.eta,o.phi,stepped);
    const int step=cell_distance(empos,seed);
    const bool is_em=o.id==11 || o.id==22;
    if(is_em && step<=2) em_directions.push_back(o.p);
    if(is_em && cell_distance(locate(o.eta,o.phi,false,false,.04),seed92)<=1)
      em92_directions.push_back(o.p);
    if (step<=2 && int(i)!=c.index) {
      if (charge[i]!=0.) {
        f.extra_track_pt=std::max(f.extra_track_pt,o.pt);
        if (step<=1) f.extra_track3_pt=std::max(f.extra_track3_pt,o.pt);
      }
      if (is_em) {
        f.core_merged+=1;
        if (o.pt>f.extra_em_pt) {f.extra_em_pt=o.pt;f.extra_em_dr=dr;}
      }
    }
    if (o.id==13) continue;
    const auto pos=is_em ? empos : locate(o.eta,o.phi,stepped,true);
    auto& bin=cells[{is_em,pos.eta_bin,pos.phi_bin}];
    bin.pos=pos;bin.energy+=o.p.e;bin.em=is_em;
  }
  for(const auto& entry:cells) {
    const auto& cell=entry.second;
    const double ds=std::hypot(cell.pos.eta-a.eta,std::remainder(cell.pos.phi-a.phi,2*M_PI));
    if(ds>=radius) continue;
    const auto empos=cell.em ? cell.pos : locate(cell.pos.eta,cell.pos.phi,stepped);
    const int step=cell_distance(empos,seed);
    const double et=cell.energy/std::cosh(cell.pos.eta);
    if (step<=2 && cell.em) {
      f.core_em_et+=et;
      if (step<=1) f.core3_em_et+=et;
      continue;
    }
    if (step<=2 && !cell.em) {f.core_had_et+=et; f.core_had_e+=cell.energy;}
    f.core_iso+=et;
    if(et>.5) f.cell_iso_et+=et;
    if(cell.energy>.5) f.cell_iso_e+=et;
  }
  f.core_em_et=std::max(0.,f.core_em_et-primary);
  f.core3_em_et=std::max(0.,f.core3_em_et-primary);
  f.mean_em_angle_rad=mean_opening(em_directions);
  f.mean92_em_angle_rad=mean_opening(em92_directions);
  return f;
}
void emit_features(std::ostream& o,const Features& f,const RawParticle& s,bool stepped) {
  const auto seed=locate(s.eta,s.phi,stepped);
  o<<"\"ieta\":"<<seed.eta_bin<<",\"iphi\":"<<seed.phi_bin<<",\"primary_et\":"<<s.p.e/std::cosh(seed.eta)
   <<",\"cone\":"<<f.cone<<",\"core_iso\":"<<f.core_iso<<",\"cell_iso_et\":"<<f.cell_iso_et
   <<",\"cell_iso_e\":"<<f.cell_iso_e<<",\"core_em_et\":"<<f.core_em_et<<",\"core3_em_et\":"<<f.core3_em_et
   <<",\"core_had_et\":"<<f.core_had_et<<",\"core_had_e\":"<<f.core_had_e
   <<",\"extra_track_pt\":"<<f.extra_track_pt<<",\"extra_track3_pt\":"<<f.extra_track3_pt
   <<",\"extra_em_pt\":"<<f.extra_em_pt<<",\"extra_em_dr\":"<<f.extra_em_dr<<",\"core_merged\":"<<f.core_merged
   <<",\"mean_em_angle_rad\":"<<f.mean_em_angle_rad<<",\"mean92_em_angle_rad\":"<<f.mean92_em_angle_rad;
}
void emit_object(std::ostream& o,const char* role,const std::vector<RawParticle>& stable,
                 const std::vector<double>& charge,const TaggedEM& c,double radius) {
  const auto& s=stable[c.index];
  o<<",\""<<role<<"\":{\"pdg\":"<<c.pdg<<",\"pt\":"<<c.reco.pt()<<",\"eta\":"<<c.reco.eta()
   <<",\"phi\":"<<c.reco.phi()<<",\"true_pt\":"<<s.pt<<",\"true_eta\":"<<s.eta<<",\"true_phi\":"<<s.phi<<",";
  emit_features(o,compute(stable,charge,c,radius,true),s,true);
  o<<",\"uniform\":{";
  emit_features(o,compute(stable,charge,c,radius,false),s,false);
  o<<"}}";
}
void process(std::ostream& o,int ev,const std::vector<TaggedEM>& accepted,
             const std::vector<RawParticle>& stable,const std::vector<double>& charge,
             double mass_low=123.3,double mass_high=126.7) {
  std::vector<TaggedEM> photons,electrons;
  for(const auto& c:accepted) (c.pdg==22 ? photons : electrons).push_back(c);
  if(photons.empty() || electrons.size()<2) return;
  for(const auto& g:photons) for(const auto& e:electrons) {
    const auto pair=g.reco+e.reco;
    if (pair.mass()<mass_low || pair.mass()>=mass_high || pair.pt()<=40.) continue;
    if (photons.size()!=1 || electrons.size()!=2) throw std::runtime_error("Ambiguous hard ee-gamma final state");
    const auto& t=electrons[electrons[0].index==e.index ? 1 : 0];
    // With precisely three EM objects, the identified tag and its fake-photon
    // outcome are exclusive. The exact event probability reduces to this product.
    const double prob=.988*(.90*.96)*(.90/33.)*kDefaultElectronEfficiency;
    const bool pass=isolated(stable,g.index,.45) && isolated(stable,e.index,.45) && isolated(stable,t.index,.30);
    o<<"{\"event\":"<<ev<<",\"mass\":"<<pair.mass()<<",\"probability\":"<<prob<<",\"nominal_probability\":"<<(pass ? prob : 0.);
    emit_object(o,"real",stable,charge,g,.45);emit_object(o,"fake",stable,charge,e,.45);emit_object(o,"tag",stable,charge,t,.30);
    o<<"}\n";
  }
}
}
#ifndef GEM_SELECTION_HELPERS_ONLY
int main(int argc,char** argv) {
  if(argc!=4) throw std::runtime_error("Use herwig_cells shower.root output.jsonl response-seed");
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
  out<<std::setprecision(16);
  for(Long64_t ev=0;ev<tree->GetEntries();++ev) {
    if(tree->GetEntry(ev)<=0 || schema!=1 || n<0 || n>10000 || n!=int(graph.object->size()) ||
       std::abs(weight-1.)>1e-8) throw std::runtime_error("Malformed event or unexpected event weight");
    std::vector<RawParticle> stable;
    std::vector<double> charge;
    for(int i=0;i<n;++i) {
      const FourVector p{particles[1][i],particles[2][i],particles[3][i],particles[0][i]};
      charge.push_back(particles[5][i]);
      stable.push_back({p,std::abs(int(particles[4][i])),p.pt(),p.eta(),p.phi()});
    }
    std::vector<TaggedEM> accepted;
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
      const auto reco=smear_em(s.p,true,smear_seed);
      const bool input=s.id==22 ? s.pt>=10. : s.pt>10.;
      if(input && reco.pt()>=20. && in_em_fiducial(reco,2.5))
        accepted.push_back({i,graph.pdg->at(graph.object->at(i)),reco});
    }
    process(out,ev,accepted,stable,charge);
  }
  std::cerr<<tree->GetEntries()<<" saved collisions\n";
}
#endif
