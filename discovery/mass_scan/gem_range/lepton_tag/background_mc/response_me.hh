// Incoming beam identity is validated by response_lhe.py before these rows
// are written. The per-photon four-vectors are the particle-level LHE inputs.
struct MEPhoton { FourVector p; bool decay; };
using MEEvents = std::vector<std::vector<MEPhoton>>;
MEEvents load_photons(const char* path, Long64_t entries, int hard, int decay) {
  std::ifstream in(path);
  std::string line;
  if(!std::getline(in,line) || line!="saved_entry,lhe_entry,decay,px,py,pz,energy")
    throw std::runtime_error("Invalid LHE photon sidecar");
  MEEvents events(entries);
  while(std::getline(in,line)) {
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream row(line);
    int saved,source,origin;
    FourVector p;
    if(!(row>>saved>>source>>origin>>p.px>>p.py>>p.pz>>p.e) || saved<0 || saved>=entries ||
       source<0 || (origin!=0 && origin!=1) || p.e<=0.) throw std::runtime_error("Invalid LHE photon record");
    events[saved].push_back({p,origin==1});
  }
  for(const auto& photons:events) {
    int nd=std::count_if(photons.begin(),photons.end(),[](const MEPhoton& p){return p.decay;});
    if(nd!=decay || int(photons.size())!=hard+decay) throw std::runtime_error("Unexpected LHE photon origins");
  }
  return events;
}

// Compare pre-recoil four-momenta anywhere along each surviving photon's own
// copy chain. Stable photon directions alone can move substantially under ISR.
// The LHE mother indices supply the production vs decay classification.
std::map<int,int> match_photons(const std::vector<RawParticle>& stable,const Ancestry& graph,
                               const std::vector<MEPhoton>& me) {
  if(graph.photon_index->size()*4!=graph.photon_p4->size())
    throw std::runtime_error("Malformed photon history");
  std::vector<int> index(graph.pdg->size(),-1);
  for(size_t i=0;i<graph.photon_index->size();++i) index.at(graph.photon_index->at(i))=i;
  // The comparison uses the maximum absolute four-component residual,
  // expressed in units of the original energy.
  constexpr double tolerance=1.e-5;
  std::vector<std::vector<std::pair<int,double>>> options(me.size(),{{-1,3.*tolerance}});
  for(size_t i=0;i<stable.size();++i) {
    if(stable[i].id!=22) continue;
    std::vector<double> distance(me.size(),1.);
    std::vector<int> todo{graph.object->at(i)};
    std::set<int> seen;
    while(!todo.empty()) {
      const int k=todo.back();
      todo.pop_back();
      if(k<0 || graph.pdg->at(k)!=22 || !seen.insert(k).second) continue;
      const int offset=index.at(k)*4;
      if(offset<0) throw std::runtime_error("Missing momentum for photon-copy node");
      const auto q=graph.photon_p4->data()+offset;
      for(size_t j=0;j<me.size();++j) {
        const auto& p=me[j].p;
        const double diff=std::max({std::abs(q[0]-p.px),std::abs(q[1]-p.py),
                                    std::abs(q[2]-p.pz),std::abs(q[3]-p.e)})/p.e;
        distance[j]=std::min(distance[j],diff);
      }
      todo.push_back(graph.previous->at(k));
      for(int t=graph.offsets->at(k);t<graph.offsets->at(k+1);++t) todo.push_back(graph.parents->at(t));
    }
    for(size_t j=0;j<me.size();++j) if(distance[j]<tolerance) options[j].push_back({i,distance[j]});
  }
  if(options.size()==1) options.push_back({{-1,0.}});
  std::map<int,int> result;
  double best=1.e6;
  for(const auto& a:options[0]) for(const auto& b:options[1]) {
    if((a.first>=0 && a.first==b.first) || a.second+b.second>=best) continue;
    best=a.second+b.second;
    result.clear();
    if(a.first>=0) result[a.first]=0;
    if(b.first>=0) result[b.first]=1;
  }
  return result;
}
