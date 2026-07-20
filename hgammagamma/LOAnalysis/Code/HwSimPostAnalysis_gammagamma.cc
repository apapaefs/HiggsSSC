#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include "TopHist.h"

using namespace std;

namespace {

struct FourVector {
  double px;
  double py;
  double pz;
  double e;

  double perp() const {
    return sqrt(px * px + py * py);
  }

  double pabs() const {
    return sqrt(px * px + py * py + pz * pz);
  }

  double eta() const {
    const double p = pabs();
    const double numerator = p + pz;
    const double denominator = p - pz;
    if (numerator <= 0.0 || denominator <= 0.0) {
      return (pz >= 0.0) ? 1.0e9 : -1.0e9;
    }
    return 0.5 * log(numerator / denominator);
  }

  double rap() const {
    const double numerator = e + pz;
    const double denominator = e - pz;
    if (numerator <= 0.0 || denominator <= 0.0) {
      return (pz >= 0.0) ? 1.0e9 : -1.0e9;
    }
    return 0.5 * log(numerator / denominator);
  }

  double phi() const {
    return atan2(py, px);
  }

  double m() const {
    const double m2 = e * e - px * px - py * py - pz * pz;
    return sqrt(max(0.0, m2));
  }
};

FourVector operator+(const FourVector& p1, const FourVector& p2) {
  return {p1.px + p2.px, p1.py + p2.py, p1.pz + p2.pz, p1.e + p2.e};
}

char* getCmdOption(char** begin, char** end, const string& option) {
  char** itr = find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const string& option) {
  return find(begin, end, option) != end;
}

string replace_extension(string path, const string& replacement) {
  const string root_ext = ".root";
  const string input_ext = ".input";

  size_t pos = path.rfind(root_ext);
  if (pos != string::npos && pos + root_ext.size() == path.size()) {
    path.replace(pos, root_ext.size(), replacement);
    return path;
  }

  pos = path.rfind(input_ext);
  if (pos != string::npos && pos + input_ext.size() == path.size()) {
    path.replace(pos, input_ext.size(), replacement);
    return path;
  }

  return path + replacement;
}

double delta_phi(const FourVector& p1, const FourVector& p2) {
  double dphi = p2.phi() - p1.phi();
  while (dphi > M_PI) {
    dphi -= 2.0 * M_PI;
  }
  while (dphi <= -M_PI) {
    dphi += 2.0 * M_PI;
  }
  return dphi;
}

double delta_r(const FourVector& p1, const FourVector& p2) {
  const double dphi = delta_phi(p1, p2);
  const double dy = p1.rap() - p2.rap();
  return sqrt(dy * dy + dphi * dphi);
}

void reset_variables(double variables[10]) {
  for (int ii = 0; ii < 10; ++ii) {
    variables[ii] = -999.0;
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (!argv[1]) {
    cout << "Use: ./HwSimPostAnalysis_gammagamma [input.root|input.input] [options]" << endl;
    cout << "Options: -t tag -n maxevents -nmin first -nmax last -w weight_scale" << endl;
    return 1;
  }

  char* infile = argv[1];

  string tag;
  if (cmdOptionExists(argv, argv + argc, "-t")) {
    tag = getCmdOption(argv, argv + argc, "-t");
    cout << "Adding tag: " << tag << endl;
    tag = "-" + tag;
  }

  double weight_scale = 1.0;
  if (cmdOptionExists(argv, argv + argc, "-w")) {
    weight_scale = atof(getCmdOption(argv, argv + argc, "-w"));
  } else if (cmdOptionExists(argv, argv + argc, "--weight-scale")) {
    weight_scale = atof(getCmdOption(argv, argv + argc, "--weight-scale"));
  }
  cout << "Using event weight scale: " << weight_scale << endl;

  TChain t("Data");

  int numPhotons = 0;
  double thePhotons[4][100];
  double evweight = 0.0;

  string stringin;
  ifstream inputlist;
  if (string(infile).find(".input") != string::npos) {
    inputlist.open(infile);
    if (!inputlist) {
      cerr << "Error: Failed to open input file " << infile << endl;
      return 1;
    }
    while (inputlist >> stringin) {
      if (stringin.empty() || stringin[0] == '#') {
        continue;
      }
      t.Add(stringin.c_str());
      cout << "Adding " << stringin << endl;
    }
    inputlist.close();
  } else if (string(infile).find(".root") != string::npos) {
    cout << "Adding " << infile << endl;
    t.Add(infile);
  } else {
    cerr << "Error: input must end in .root or .input" << endl;
    return 1;
  }

  t.SetBranchAddress("thePhotons", &thePhotons);
  t.SetBranchAddress("numPhotons", &numPhotons);
  t.SetBranchAddress("evweight", &evweight);

  const int EventNumber = int(t.GetEntries());
  cout << "Total number of events in " << infile << " : " << EventNumber << endl;
  if (EventNumber == 0) {
    cerr << "Error: no events found" << endl;
    return 1;
  }

  int minevents = 0;
  int maxevents = EventNumber;

  if (cmdOptionExists(argv, argv + argc, "-n")) {
    maxevents = atoi(getCmdOption(argv, argv + argc, "-n"));
    if (maxevents > EventNumber) {
      maxevents = EventNumber;
    }
    cout << "Analyzing up to " << maxevents << endl;
  }

  if (cmdOptionExists(argv, argv + argc, "-nmax") && !cmdOptionExists(argv, argv + argc, "-n")) {
    maxevents = atoi(getCmdOption(argv, argv + argc, "-nmax"));
    if (maxevents > EventNumber) {
      maxevents = EventNumber;
    }
    cout << "Analyzing up to " << maxevents << endl;
  }

  if (cmdOptionExists(argv, argv + argc, "-nmin")) {
    minevents = atoi(getCmdOption(argv, argv + argc, "-nmin"));
    if (minevents < 0) {
      minevents = 0;
    }
    if (minevents > maxevents) {
      minevents = 0;
    }
    cout << "Analyzing from " << minevents << endl;
  }

  if (maxevents <= minevents) {
    cerr << "Error: requested event range is empty" << endl;
    return 1;
  }

  const string output = replace_extension(infile, tag + ".top");
  const string output_dat = replace_extension(infile, tag + ".dat");
  const string output_evp = replace_extension(infile, tag + ".evp");
  const string output_root = replace_extension(infile, tag + "_var.root");

  ofstream outdat(output_dat.c_str(), ios::out);
  ofstream outevp(output_evp.c_str(), ios::out);

  TFile* dat2 = new TFile(output_root.c_str(), "RECREATE");
  TTree* Data2 = new TTree("Data2", "Gamma Gamma Data Tree");

  double variables[10];
  double eventweight[1];
  Data2->Branch("variables", variables, "variables[10]/D");
  Data2->Branch("eventweight", eventweight, "eventweight[1]/D");

  const double cut_pt_photon = 10.0;
  const double cut_eta_photon = 6.0;

  TopHist h_dummy(10, output, "dummy histo", 0, 1);
  TopHist h_nphotons(11, output, "number of selected photons", -0.5, 10.5);
  TopHist h_pT_photons(80, output, "pT of selected photons", 0, 2000);
  TopHist h_eta_photons(60, output, "eta of selected photons", -6, 6);
  TopHist h_pT_gamma1(80, output, "pT of leading photon", 0, 2000);
  TopHist h_pT_gamma2(80, output, "pT of subleading photon", 0, 2000);
  TopHist h_eta_gamma1(60, output, "eta of leading photon", -6, 6);
  TopHist h_eta_gamma2(60, output, "eta of subleading photon", -6, 6);
  TopHist h_mgg(100, output, "diphoton invariant mass", 0, 1000);
  TopHist h_deltaR_gg(80, output, "DeltaR of two leading photons", 0, 12);
  TopHist h_deltaPhi_gg(64, output, "DeltaPhi of two leading photons", -3.2, 3.2);
  TopHist h_pT_gg(80, output, "pT of diphoton system", 0, 2000);
  TopHist h_y_gg(60, output, "rapidity of diphoton system", -6, 6);

  double sum_weight = 0.0;
  double sum_diphoton_weight = 0.0;
  int selected_events = 0;

  for (int ii = minevents; ii < maxevents; ++ii) {
    t.GetEntry(ii);
    if (ii % 1000 == 0) {
      cout << "Event number: " << ii << "\r" << flush;
    }

    const double scaled_weight = evweight * weight_scale;
    sum_weight += scaled_weight;

    vector<FourVector> photons;
    for (int jj = 0; jj < numPhotons; ++jj) {
      FourVector photon{thePhotons[1][jj], thePhotons[2][jj], thePhotons[3][jj], thePhotons[0][jj]};
      if (photon.perp() < cut_pt_photon) {
        continue;
      }
      if (fabs(photon.eta()) > cut_eta_photon) {
        continue;
      }
      photons.push_back(photon);
    }

    sort(photons.begin(), photons.end(),
         [](const FourVector& left, const FourVector& right) {
           return left.perp() > right.perp();
         });

    reset_variables(variables);
    eventweight[0] = scaled_weight;
    variables[9] = photons.size();

    h_nphotons.thfill(photons.size(), scaled_weight);
    for (size_t jj = 0; jj < photons.size(); ++jj) {
      h_pT_photons.thfill(photons[jj].perp(), scaled_weight);
      h_eta_photons.thfill(photons[jj].eta(), scaled_weight);
    }

    if (!photons.empty()) {
      h_pT_gamma1.thfill(photons[0].perp(), scaled_weight);
      h_eta_gamma1.thfill(photons[0].eta(), scaled_weight);
      variables[1] = photons[0].perp();
      variables[2] = photons[0].eta();
    }

    if (photons.size() > 1) {
      h_pT_gamma2.thfill(photons[1].perp(), scaled_weight);
      h_eta_gamma2.thfill(photons[1].eta(), scaled_weight);
      variables[3] = photons[1].perp();
      variables[4] = photons[1].eta();

      const FourVector diphoton = photons[0] + photons[1];
      const double mgg = diphoton.m();
      const double dRgg = delta_r(photons[0], photons[1]);
      const double dPhigg = delta_phi(photons[0], photons[1]);

      h_mgg.thfill(mgg, scaled_weight);
      h_deltaR_gg.thfill(dRgg, scaled_weight);
      h_deltaPhi_gg.thfill(dPhigg, scaled_weight);
      h_pT_gg.thfill(diphoton.perp(), scaled_weight);
      h_y_gg.thfill(diphoton.rap(), scaled_weight);

      variables[0] = mgg;
      variables[5] = dRgg;
      variables[6] = dPhigg;
      variables[7] = diphoton.perp();
      variables[8] = diphoton.rap();

      sum_diphoton_weight += scaled_weight;
      selected_events++;
      outevp << ii << endl;
    }

    Data2->Fill();
  }

  Data2->Write();
  dat2->Close();
  cout << endl;
  cout << "A root tree has been written to the file: " << output_root << endl;

  h_dummy.plot(1, 0);
  h_nphotons.add(output, 0, 0);
  h_pT_photons.add(output, 0, 0);
  h_eta_photons.add(output, 0, 0);
  h_pT_gamma1.add(output, 0, 0);
  h_eta_gamma1.add(output, 0, 0);
  h_pT_gamma2.add(output, 0, 0);
  h_eta_gamma2.add(output, 0, 0);
  h_mgg.add(output, 0, 0);
  h_deltaR_gg.add(output, 0, 0);
  h_deltaPhi_gg.add(output, 0, 0);
  h_pT_gg.add(output, 0, 0);
  h_y_gg.add(output, 0, 0);

  outdat << "# HwSimPostAnalysis_gammagamma summary" << endl;
  outdat << "input " << infile << endl;
  outdat << "analysis legacy_direct_photons" << endl;
  outdat << "detector_response none" << endl;
  outdat << "response_mode genuine" << endl;
  outdat << "weighted_hypotheses 0" << endl;
  outdat << "events_read " << (maxevents - minevents) << endl;
  outdat << "events_with_two_selected_photons " << selected_events << endl;
  outdat << "weight_scale " << setprecision(12) << weight_scale << endl;
  outdat << "sum_weight " << setprecision(12) << sum_weight << endl;
  outdat << "sum_diphoton_weight " << setprecision(12) << sum_diphoton_weight << endl;

  cout << "------------------" << endl;
  cout << "events with >= 2 selected photons = " << selected_events << endl;
  cout << "sum weight = " << sum_weight << endl;
  cout << "sum diphoton weight = " << sum_diphoton_weight << endl;
  cout << "------------------" << endl;

  return 0;
}
