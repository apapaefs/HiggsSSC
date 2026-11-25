// Apply the same tagged EM smearing, cracks, 20 GeV threshold, and isolation
// as the main analysis to a native PYTHIA 5.6 / JETSET 7.3 LUJETS record.
#define TAGGED_BACKGROUND_RESPONSE
#include "../../signal/measure.cc"
#include <array>
#include <iomanip>
#include <iostream>

extern "C" struct {int n, k[5][4000]; float p[5][4000], v[5][4000];} lujets_;
namespace {
constexpr int NWindows = 3, NStages = 6;
constexpr double bounds[NWindows][2] = {{75.,165.},{80.,140.},{123.3,126.7}};
const char* stages[NStages] = {"kinematics","fiducial","isolated","pair_pt","ideal92_iso","ideal92_pair_pt"};
double sum[NWindows][NStages]{}, sum2[NWindows][NStages]{};
long matched = 0, total = 0;
constexpr double probability = .988 * .90 * (.90 / 33.) * (.90 * .96);

// Separate 1992 Table 15 diagnostic. This changes the geometrical and
// isolation choices only. The GEANT cluster shape and its correlations
// cannot be recovered from a stable particle event record.
bool isolated92(const std::vector<RawParticle>& stable, int index, double radius) {
  const auto& candidate = stable.at(index);
  const double limit = 5. + .10*candidate.pt;
  double sum_pt = 0.;
  for (size_t j=0;j<stable.size();++j) {
    const auto& other=stable[j];
    if (j==size_t(index) || other.id==12 || other.id==13 || other.id==14 || other.id==16) continue;
    if (std::hypot(other.eta-candidate.eta,std::remainder(other.phi-candidate.phi,2*M_PI))<radius)
      sum_pt+=other.pt;
  }
  return sum_pt<limit;
}

// In subprocess 19 entries 7/8 of the hard record are Z and gamma,
// entries 9/10 are the Z daughters. Trace stable particles back along
// same-PDG copies. Photons radiated by charged leptons do not qualify.
int hard_ancestor(int i) {
  const auto& k=lujets_.k;
  const int pdg=k[1][i];
  while (i>=0) {
    if (k[0][i]==21 && ((i==6 || i==7) && pdg==22)) return 0;
    // After QED radiation JETSET can assign *both* physical charged
    // daughters to the first hard Z decay slot. Charge of the final
    // electron distinguishes the two. The hard record must identify
    // a charged lepton from this Z before either slot can be accepted.
    if (k[0][i]==21 && (i==8 || i==9) && std::abs(k[1][i])==11 && std::abs(pdg)==11
        && (k[2][i]==7 || k[2][i]==8) && k[1][k[2][i]-1]==23) return pdg>0 ? 1 : 2;
    if (k[1][i]!=pdg) break;
    const int next=k[2][i]-1;
    if (next>=i) break;
    i=next;
  }
  return -1;
}
}

extern "C" void measure_(const int* seed, const int* event) {
  const auto& rec=lujets_;
  if (rec.n<10 || rec.n>4000) throw std::runtime_error("Unexpected LUJETS length");
  ++total;
  std::array<int,3> record = {-1,-1,-1};
  for (int i=0;i<rec.n;++i) {
    if (rec.k[0][i]<=0 || rec.k[0][i]>10) continue;
    const int id=std::abs(rec.k[1][i]);
    if (id!=22 && id!=11) continue;
    const int idx=hard_ancestor(i);
    if (idx<0) continue;
    if (record[idx]>=0) throw std::runtime_error("Ambiguous hard-object copy");
    record[idx]=i;
  }
  if (*std::min_element(record.begin(),record.end())<0) return;
  ++matched;
  std::array<FourVector,3> reco;
  for (int a=0;a<3;++a) {
    const int i=record[a];
    FourVector p{rec.p[0][i],rec.p[1][i],rec.p[2][i],rec.p[3][i]};
    if (p.pt()<=10.) return;
    reco[a]=smear_em(p,true,object_seed(*seed,*event,a==0 ? 90U : 91U,i));
    if (reco[a].pt()<20. || std::abs(reco[a].eta())>=2.5) return;
  }
  std::vector<RawParticle> stable;
  std::array<int,3> index{};
  for (int i=0;i<rec.n;++i) {
    if (rec.k[0][i]<=0 || rec.k[0][i]>10) continue;
    for (int a=0;a<3;++a) if (i==record[a]) index[a]=stable.size();
    FourVector p{rec.p[0][i],rec.p[1][i],rec.p[2][i],rec.p[3][i]};
    stable.push_back({p,std::abs(rec.k[1][i]),p.pt(),p.eta(),p.phi()});
  }
  double in_event[NWindows][NStages]{};
  const bool fiducial=std::all_of(reco.begin(),reco.end(),[](const FourVector& p){return in_em_fiducial(p,2.5);});
  for (int fake=1;fake<=2;++fake) {
    const auto pair=reco[0]+reco[fake];
    const bool iso=fiducial && isolated(stable,index[0],.45) && isolated(stable,index[fake],.45)
                            && isolated(stable,index[3-fake],.30);
    const bool iso92=isolated92(stable,index[0],.45) && isolated92(stable,index[fake],.45)
                      && isolated92(stable,index[3-fake],.30);
    const bool pass[NStages]={true,fiducial,iso,iso && pair.pt()>40.,iso92,iso92 && pair.pt()>40.};
    for (int w=0;w<NWindows;++w) if (pair.mass()>=bounds[w][0] && pair.mass()<bounds[w][1]) {
      for (int j=0;j<NStages;++j) if(pass[j]) ++in_event[w][j];
    }
  }
  for(int w=0;w<NWindows;++w) for(int j=0;j<NStages;++j) {
    const double v=in_event[w][j];
    sum[w][j]+=v;
    sum2[w][j]+=v*v;
  }
}

extern "C" void finish_() {
  std::cout<<"MATCH "<<matched<<" "<<total<<"\n";
  std::cout<<"PROB "<<std::setprecision(14)<<probability<<"\n";
  for(int w=0;w<NWindows;++w) for(int j=0;j<NStages;++j)
    std::cout<<"BIN "<<w<<" "<<stages[j]<<" "<<sum[w][j]<<" "<<sum2[w][j]<<"\n";
}
