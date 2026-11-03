#ifndef HWSIM_ANCESTRY_H
#define HWSIM_ANCESTRY_H
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <vector>

// Schema 1 records actual parent links and ThePEG copy links separately.
// Do not walk through beam ancestors to infer promptness: every pp particle
// ultimately has a proton ancestor. Classify the first non-copy parent.
namespace hwsim {
enum Origin { Hard = 1, Boson = 2, Higgs = 4, Tau = 8, Hadron = 16, Other = 32 };
struct Ancestry {
  std::vector<int> lineage, origins;
  std::vector<std::vector<int>> emitters;

  Ancestry(const std::vector<int>& pdg, const std::vector<int>& parents,
           const std::vector<int>& offsets, const std::vector<int>& previous,
           const std::vector<int>& hard) {
    const int n = pdg.size();
    if (previous.size() != pdg.size() || hard.size() != pdg.size() ||
        offsets.size() != pdg.size() + 1 || offsets.front() != 0 ||
        offsets.back() != static_cast<int>(parents.size()) ||
        !std::is_sorted(offsets.begin(), offsets.end()))
      throw std::runtime_error("invalid HwSim ancestry graph dimensions");
    for (int p : parents) if (p < 0 || p >= n)
      throw std::runtime_error("invalid HwSim ancestry parent index");
    for (int p : previous) if (p < -1 || p >= n)
      throw std::runtime_error("invalid HwSim previous-copy index");
    lineage.resize(n);
    std::iota(lineage.begin(), lineage.end(), 0);
    const auto root = [&](int a) {
      while (lineage[a] != a) { lineage[a] = lineage[lineage[a]]; a = lineage[a]; }
      return a;
    };
    const auto join = [&](int a, int b) {
      if (b < 0 || pdg[a] != pdg[b]) return;
      a = root(a); b = root(b);
      lineage[std::max(a, b)] = std::min(a, b);
    };
    for (int i = 0; i < n; ++i) {
      join(i, previous[i]);
      for (int j = offsets[i]; j < offsets[i + 1]; ++j) join(i, parents[j]);
    }
    for (int i = 0; i < n; ++i) lineage[i] = root(i);
    origins.assign(n, 0); emitters.resize(n);
    for (int i = 0; i < n; ++i) {
      const int group = lineage[i];
      if (hard[i] == 1) origins[group] |= Hard;
      for (int j = offsets[i]; j < offsets[i + 1]; ++j) {
        const int parent = parents[j], id = std::abs(pdg[parent]);
        if (lineage[parent] == group) continue;
        origins[group] |= id == 25 ? Higgs : (id == 23 || id == 24) ? Boson
                         : id == 15 ? Tau : id > 100 ? Hadron : Other;
        if (pdg[i] == 22 && (id == 11 || id == 13)) emitters[group].push_back(lineage[parent]);
      }
    }
    for (auto& e : emitters) {
      std::sort(e.begin(), e.end()); e.erase(std::unique(e.begin(), e.end()), e.end());
    }
  }
  bool prompt(int node) const {
    const int flags = origins.at(lineage.at(node));
    return (flags & (Hard | Boson | Higgs)) && !(flags & (Tau | Hadron));
  }
};
}
#endif
