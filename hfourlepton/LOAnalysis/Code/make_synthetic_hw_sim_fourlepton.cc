#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>

namespace {

constexpr int kMaxParticles = 10000;
constexpr int kMaxRecoObjects = 100;
constexpr double kPi = 3.14159265358979323846;

void add_particle(double objects[][kMaxParticles], int& count, int pdg_id,
                  double pt, double eta, double phi, double mass = 0.0) {
  if (count >= kMaxParticles) {
    throw std::runtime_error("synthetic particle buffer overflow");
  }
  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = pt * std::sinh(eta);
  const double energy =
      std::sqrt(px * px + py * py + pz * pz + mass * mass);
  objects[0][count] = energy;
  objects[1][count] = px;
  objects[2][count] = py;
  objects[3][count] = pz;
  objects[4][count] = pdg_id;
  objects[5][count] =
      std::abs(pdg_id) == 11 || std::abs(pdg_id) == 13
          ? (pdg_id > 0 ? -1.0 : 1.0)
          : 0.0;
  ++count;
}

void set_reco_object(double objects[][kMaxRecoObjects], int index, double pt,
                     double eta, double phi) {
  objects[0][index] = pt * std::cosh(eta);
  objects[1][index] = pt * std::cos(phi);
  objects[2][index] = pt * std::sin(phi);
  objects[3][index] = pt * std::sinh(eta);
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    if (argc < 2 || argc > 3) {
      std::cerr
          << "Usage: make_synthetic_hw_sim_fourlepton OUTPUT.root "
             "[minimal|unnamed-weights]\n";
      return 1;
    }
    const std::string output_name = argv[1];
    const bool minimal = argc == 3 && std::string(argv[2]) == "minimal";
    const bool unnamed_weights =
        argc == 3 && std::string(argv[2]) == "unnamed-weights";
    if (argc == 3 && !minimal && !unnamed_weights) {
      throw std::runtime_error(
          "the supported modes are 'minimal' and 'unnamed-weights'");
    }

    TFile output(output_name.c_str(), "RECREATE");
    if (output.IsZombie()) {
      throw std::runtime_error("failed to create " + output_name);
    }
    TTree tree("Data", "Synthetic HwSim Data tree for four-lepton tests");

    int numparticles = 0;
    double objects[8][kMaxParticles] = {};
    double evweight = 1.0;
    tree.Branch("numparticles", &numparticles, "numparticles/I");
    tree.Branch("objects", objects, "objects[8][10000]/D");
    tree.Branch("evweight", &evweight, "evweight/D");

    int numJets = 0;
    double theJets[5][kMaxRecoObjects] = {};
    int numbJets = 0;
    double thebJets[5][kMaxRecoObjects] = {};
    double theETmiss[4] = {};
    std::vector<double> theOptWeights;
    std::vector<std::string> theOptWeightsNames;
    if (!minimal) {
      tree.Branch("numJets", &numJets, "numJets/I");
      tree.Branch("theJets", theJets, "theJets[5][100]/D");
      tree.Branch("numbJets", &numbJets, "numbJets/I");
      tree.Branch("thebJets", thebJets, "thebJets[5][100]/D");
      tree.Branch("theETmiss", theETmiss, "theETmiss[4]/D");
      tree.Branch("theOptWeights", &theOptWeights);
      tree.Branch("theOptWeightsNames", &theOptWeightsNames);
    }

    for (int event = 0; event < 8; ++event) {
      numparticles = 0;
      std::fill(&objects[0][0], &objects[0][0] + 8 * kMaxParticles, 0.0);
      std::fill(&theJets[0][0],
                &theJets[0][0] + 5 * kMaxRecoObjects, 0.0);
      std::fill(&thebJets[0][0],
                &thebJets[0][0] + 5 * kMaxRecoObjects, 0.0);
      std::fill(std::begin(theETmiss), std::end(theETmiss), 0.0);

      if (event == 0) {
        // 2e2mu with a photon dressed into the negative electron.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, -13, 20.0, 0.0, 1.0 + kPi,
                     0.105658);
        add_particle(objects, numparticles, 22, 2.0, 0.01, 0.03);
        evweight = 1.0;
      } else if (event == 1) {
        // 4e plus a fifth isolated electron, exercising best-quadruplet logic.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 11, 20.0, 0.0, 1.0);
        add_particle(objects, numparticles, -11, 20.0, 0.0, 1.0 + kPi);
        add_particle(objects, numparticles, 11, 11.0, 2.0, 2.0);
        evweight = 2.0;
      } else if (event == 2) {
        // Negative-weight 4mu source event.
        add_particle(objects, numparticles, 13, 45.0, 0.0, 0.0, 0.105658);
        add_particle(objects, numparticles, -13, 45.0, 0.0, kPi, 0.105658);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, -13, 20.0, 0.0, 1.0 + kPi,
                     0.105658);
        evweight = -1.0;
      } else if (event == 3) {
        // Exactly pT=10 fails the strict analysis threshold.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, -13, 10.0, 0.0, 1.0 + kPi,
                     0.105658);
        evweight = 1.0;
      } else if (event == 4) {
        // One electron in the open 1.01<|eta|<1.16 crack.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 11, 20.0, 1.10, 1.0);
        add_particle(objects, numparticles, -11, 20.0, -1.10, 1.0 + kPi);
        evweight = 1.0;
      } else if (event == 5) {
        // A 6 GeV charged constituent inside one electron's isolation cone.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, -13, 20.0, 0.0, 1.0 + kPi,
                     0.105658);
        add_particle(objects, numparticles, 211, 6.0, 0.05, 0.15, 0.13957);
        evweight = 1.0;
      } else if (event == 6) {
        // Four reconstructed leptons but no opposite-sign pair.
        add_particle(objects, numparticles, 11, 45.0, 0.0, 0.0);
        add_particle(objects, numparticles, 11, 45.0, 0.0, kPi);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, 13, 20.0, 0.0, 1.0 + kPi,
                     0.105658);
        evweight = 1.0;
      } else {
        // Valid OSSF pairs with Z1 outside its analysis window.
        add_particle(objects, numparticles, 11, 30.0, 0.0, 0.0);
        add_particle(objects, numparticles, -11, 30.0, 0.0, kPi);
        add_particle(objects, numparticles, 13, 15.0, 0.0, 1.0, 0.105658);
        add_particle(objects, numparticles, -13, 15.0, 0.0, 1.0 + kPi,
                     0.105658);
        evweight = 1.0;
      }

      if (!minimal) {
        numJets = 1;
        set_reco_object(theJets, 0, 30.0 + event, 0.5, 2.5);
        numbJets = event == 1 ? 1 : 0;
        if (numbJets != 0) {
          set_reco_object(thebJets, 0, 25.0, 0.8, -2.0);
        }
        const double met_pt = 15.0 + event;
        const double met_phi = -0.7;
        theETmiss[0] = met_pt;
        theETmiss[1] = met_pt * std::cos(met_phi);
        theETmiss[2] = met_pt * std::sin(met_phi);
        theOptWeights = {evweight * 0.9, evweight * 1.1};
        theOptWeightsNames =
            unnamed_weights ? std::vector<std::string>{}
                            : std::vector<std::string>{"scale_down", "scale_up"};
      }
      tree.Fill();
    }

    output.cd();
    tree.Write();
    output.Close();
    std::cout << "Wrote " << output_name
              << (minimal ? " (required branches only)" : " (full optional schema)")
              << std::endl;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return 1;
  }
}
