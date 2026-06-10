#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

//ROOT include files
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

//Fastjet headers
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include <fastjet/tools/JHTopTagger.hh>
#include <fastjet/Selector.hh>

//Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>

//custom headers
#include "TopHist.h"
#include "complex_d.h"

using namespace std;
using namespace fastjet;

//----------------------------------------------------------------------
// Some four-vector operators
//----------------------------------------------------------------------
double dot(fastjet::PseudoJet p1, fastjet::PseudoJet p2);
double deltaR(fastjet::PseudoJet p1, fastjet::PseudoJet p2);

/* jet to lepton mistag */
double Pjet_to_lepton(double pt);

/* jet to photon mistag */
double Pjet_to_photon(double pt);

//----------------------------------------------------------------------// forward declaration for printing out info about a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream &, const PseudoJet &);

//----------------------------------------------------------------------
// command line parameters
//----------------------------------------------------------------------
char* getCmdOption(char ** begin, char ** end, const std::string & option);
bool cmdOptionExists(char** begin, char** end, const std::string& option);

//----------------------------------------------------------------------
// Analysis functions
//----------------------------------------------------------------------

// smearing of jets, leptons and photons.
fastjet::PseudoJet smear_jet(fastjet::PseudoJet jet_in);
fastjet::PseudoJet smear_lepton(fastjet::PseudoJet lepton_in, int lepton_id);
fastjet::PseudoJet smear_photon(fastjet::PseudoJet photon_in);

// acceptance efficiency for leptons, photons, jets
bool lepton_efficiency_accept(fastjet::PseudoJet lepton_in, int lepton_id);
bool photon_efficiency_accept(fastjet::PseudoJet photon_in);
bool jet_efficiency_accept(fastjet::PseudoJet jet_in);

// IDs of B-hadrons used by the btag_hadrons function
int bhadronid[105] = {5122, -5122, 15122, -15122, 5124, -5124, 5334, -5334, 5114, -5114, 5214, -5214, 5224, -5224, 5112, -5112, 5212, -5212, 5222, -5222, 15322, -15322, 15312, -15312, 15324, -15324, 15314, -15314, 5314, -5314, 5324, -5324, 5132, -5132, 5232, -5232, 5312, -5312, 5322, -5322, 551, 10555, 100551, 200551, 553, 557, 555, 100555, 200555, 20523, -20523, 20513, -20513, 20543, -20543, 20533, -20533, 511, 521, -511, -521, 531, -531, 541, -541, 513, 523, -513, -523, 533, -533, 543, -543, 10513, 10523, -10513, -10523, 10533, -10533, 10543, -10543, 10511, 10521, -10511, -10521, 10531, -10531, 10541, -10541, 20513, 20523, -20513, -20523, 20533, -20533, 20543, -20543, 515, 525, -515, -525, 535, -535, 545, -545};

/*
 * CREATE ROOT CHAIN TO READ IN THE FILES
 */
TChain t("Data");


/*
 * DECLARE RANDOM NUMBERS
 */
TRandom3 rnd;
TRandom3 rndint;

/*****
 ***** SWITCHES FOR SMEARING/EFFICIENCIES
 *****/
bool donotsmear_jets = 0;
bool donotsmear_leptons = 0;
bool donot_apply_efficiency = 0;
bool donotsmear_photons = 0;

int main(int argc, char *argv[]) {

  //take command line options
  char* output;
  char* infile = "";
  if(argv[1]) { infile = argv[1]; } else { cout << "Use: ./HwSimAnalysis [input] [options]" << endl; exit(1); }

  //set the variables and addresses to be read from root file
  //total number of particles in an event
  int numparticles;

  /** particle information in the order:
   * 4 momenta (E,x,y,z), id, other info
   **/
  double objects[8][10000];
  /* the event weight */
  double evweight;
  /* The missing energy four-vector */
  double theETmiss[4];
  /* The incoming partons */
  double incoming[5][2];
  /* the outgoing partons */
  double partons[5][100];



  /* the optional weight values */
  std::vector<double> *theOptWeights;

  /* the optional weight names */
  std::vector<string> *theOptWeightsNames;


  /*
   * SET THE ROOT BRANCH ADDRESSES
   */
  t.SetBranchAddress("numparticles",&numparticles);
  t.SetBranchAddress("objects",&objects);
  t.SetBranchAddress("evweight",&evweight);
  t.SetBranchAddress("partons",&partons);
  t.SetBranchAddress("incoming",&incoming);


  t.SetBranchAddress("theETmiss", &theETmiss);

  // uncomment if you wish to use optional weights:
  /* t.SetBranchAddress("theOptWeightsNames", &theOptWeightsNames);
     t.SetBranchAddress("theOptWeights", &theOptWeights);*/

  /* Set up random number
   * generator
   */
  rnd.SetSeed(14101983);

  /* Add up all the input
   * files to the chain
   */
  string stringin = "";
  ifstream inputlist;
  if (std::string(infile).find(".input") != std::string::npos) {
    inputlist.open(infile);
    if(!inputlist) {  cerr << "Error: Failed to open input file " << infile << endl; exit(1); }
    while(inputlist) {
      inputlist >> stringin;
      if(stringin!="") { t.Add(stringin.c_str());
        cout << "Adding " << stringin.c_str() << endl;
      }
      stringin = "";
    }
    inputlist.close();
  } else if (std::string(infile).find(".root") != std::string::npos) {
    cout << "Adding " << infile << endl;
    t.Add(infile);
  }

  /* Get Number of events
   * and print
   */
  int EventNumber(int(t.GetEntries()));
  cout << "Total number of events in " << infile << " : " << EventNumber << endl;

  /*
   * -b: USED TO REANALYZE PREVIOUSLY PASSED EVENTS ONLY, DEFAULT IS ALL EVENTS
   */

  //whether the analysis performed is level-2 or level-3
  bool basic = true;
  if(cmdOptionExists(argv, argv+argc, "-b")) {
    cout << "Looking for .evp2 file, running over all events" << endl;
    basic = false;
  }


  /*
   * -t: ADD AN EXTENSION TAG TO YOUR OUTPUT FILES
   */
  string tag;
  tag = "";
  if(cmdOptionExists(argv, argv+argc, "-t")) {
    tag = getCmdOption(argv, argv + argc, "-t");
    cout << "Adding tag: " << tag << endl;
    tag = "-" + tag;
  }

  /*
   * -n: RUN FROM START OF FILE UP TO A GIVEN NUMBER OF EVENTS
   */
  char * switch_maxevents;
  char * switch_minevents;
  int maxevents(0), minevents(0);
  if(cmdOptionExists(argv, argv+argc, "-n")) {
    switch_maxevents = getCmdOption(argv, argv + argc, "-n");
    maxevents=(atoi(switch_maxevents));
    if(maxevents > EventNumber) { maxevents = EventNumber; }
    cout << "Analyzing up to " << maxevents << endl;
    if(maxevents < 1 || maxevents > 1E10) { cout << "Error: maxevents must be in the range [1,1E10]" << endl; exit(1); }
  }

  /*
   * -nmax: RUN FROM START OF FILE UP TO A GIVEN NUMBER OF EVENTS, TO BE USED IN CONJUNCTION WITH -nmin
   */
  //maximum number of events to analyze
  if(cmdOptionExists(argv, argv+argc, "-nmax") && !cmdOptionExists(argv, argv+argc, "-n")) {
   switch_maxevents = getCmdOption(argv, argv + argc, "-nmax");
    maxevents=(atoi(switch_maxevents));
    if(maxevents > EventNumber) { maxevents = EventNumber; }
    cout << "Analyzing up to " << maxevents << endl;
    if(maxevents < 1 || maxevents > 1E10) { cout << "Error: maxevents must be in the range [1,1E10]" << endl; exit(1); }
  }
  if(!cmdOptionExists(argv, argv+argc, "-nmax") && !cmdOptionExists(argv, argv+argc, "-n")) { maxevents = EventNumber; }

  /*
   * -nmin: RUN FROM POINT nmin OF FILE UP TO A GIVEN NUMBER OF EVENTS SPECIFIED BY -nmax
   */
  //starting number of events to analyse
  if(cmdOptionExists(argv, argv+argc, "-nmin")) {
    switch_minevents = getCmdOption(argv, argv + argc, "-nmin");
    minevents=(atoi(switch_minevents));
    if(minevents > maxevents) { minevents = 0; }
    cout << "Analyzing from " << minevents << endl;
    if(minevents < 1 || minevents > 1E10) { cout << "Error: minevents must be in the range [1,1E10]" << endl; exit(1); }
  }


  /*
   * CREATE THE OUTPUT FILE STRINGS
   */
  string outnew = "";
  outnew = std::string(infile);
  string replacement = tag + ".top";
  boost::replace_all(outnew, ".root", replacement);
  boost::replace_all(outnew, ".input", replacement);
  char* output2 = new char[outnew.length() + 1];
  //  cout << outnew.c_str() << endl;
  strcpy (output2, outnew.c_str());
  output = output2;

  char* output_dat;
  string outnew2 = "";
  outnew2 = std::string(infile);
  replacement = tag + ".dat";
  boost::replace_all(outnew2, ".root", replacement);
  boost::replace_all(outnew2, ".input", replacement);
  char* output3 = new char[outnew2.length() + 1];
  strcpy (output3, outnew2.c_str());
  output_dat = output3;
  ofstream outdat(output_dat, ios::out);

  //load events that have passed the second stage of analysis
  //if basic = false;
  string ineventpass;
  ifstream inevt;
  string inevt_curr;
  int passed_event[20000];
  int npassed_previous(0);
  if(basic == false) {
    ineventpass = std::string(infile);
    replacement = tag + ".evp";
    boost::replace_all(ineventpass,".input", replacement);
    boost::replace_all(ineventpass,".root", replacement);
    inevt.open(ineventpass.c_str());
    if(!inevt) { cerr << "Error: Cannot open "<< ineventpass.c_str() << endl; exit(1); }
    for(int ii = 0; ii < 1000; ii++) { passed_event[ii] = -1; }
    while(inevt) {
      inevt >> inevt_curr;
      // cout << inevt_curr.c_str() << endl;
      passed_event[npassed_previous] = atoi(inevt_curr.c_str());
      npassed_previous++;
    }
  }
  //for(int pp = 0; pp < npassed_previous; pp++) { coust << passed_event[pp] << endl; }

  string outeventpass = "";
  ofstream outevp;

  if(basic == false) {
    outeventpass = std::string(infile);
    replacement = tag + ".evp2";
    boost::replace_all(outeventpass,".root", replacement);
    boost::replace_all(outeventpass,".input", replacement);
    boost::replace_all(outeventpass,".top", replacement);
    outevp.open(outeventpass.c_str());
  } else if(basic == true) {
    outeventpass = std::string(infile);
    replacement = tag + ".evp";
    boost::replace_all(outeventpass,".root", replacement);
    boost::replace_all(outeventpass,".input", replacement);
    boost::replace_all(outeventpass,".top", replacement);
    outevp.open(outeventpass.c_str());
  }

  /*
   * PREPARES THE OUTPUT ARRAY FOR *_var.root: USED FOR FURTHER ANALYSIS
   */
  std::cout << "Preparing Root Tree for event variables" << endl;
  TTree* Data2;
  TFile* dat2;
  string fnameroot = std::string(infile);
  replacement = tag + "_var.root";
  boost::replace_all(fnameroot,".root", replacement);
  boost::replace_all(fnameroot,".input", replacement);
  dat2 = new TFile(fnameroot.c_str(), "RECREATE");
  Data2 = new TTree ("Data2", "Data Tree");
  //variables to fill in the .root file
  double variables[10];
  double eventweight[1];
  double muonevent[1];
  Data2->Branch("variables", &variables, "variables[10]/D");
  Data2->Branch("eventweight", &eventweight, "eventweight[1]/D");

  /*
   * CUTS DEFINED HERE IN GEV
   */

  /* ALL PARTICLE CUTS */
  double eta_cut(5.0); //global pseudorapidity cut of particles
  double pt_cut_part(0.1); //global pt cut for particles

  /* JET CUTS */
  double cut_pt_jet(20.0); //pt cut for jets
  double cut_eta_jet(5.0); //pseudo-rapidity cut for jets

  /* LEPTON CUTS */
  double cut_eta_electron(3.0); // pseudo-rapidity cut for electrons
  double cut_pt_electron(10.0); // pt cut for electrons

  double cut_eta_muon(3.0); // pseudo-rapidity cut for muons
  double cut_pt_muon(10.0); // pt cut for muons

  /* PHOTON CUTS */
  double cut_pt_photon(20.0); //pt cut for jets
  double cut_eta_photon(5.0); //pseudo-rapidity cut for jets

  /*
   * COUNTERS FOR NUMBER OF EVENTS THAT PASS CUTS
   */
  double passcuts(0); //passed all cuts

  /*
   * PARAMETERS AND SWITCHES
   */

  /*
   * HISTOGRAMS DEFINED HERE
   */
  TopHist h_dummy(10,output,"dummy histo", 0,1);
  TopHist h_pT_leptons(50,output,"pT of leptons",0, 1000);
  TopHist h_m4l(50,output,"invariant mass of 4 leptons",0,1000);
  TopHist h_ETmiss(50,output,"missing transverse energy",0,1000);


  /*
   *
   * LOOP OVER EVENTS
   * AND
   * PERFORM ANALYSIS
   *
   */
  bool perform_analysis_on_event = false;
  for(int ii = minevents; ii < maxevents; ii++) {

    /* IF LEVEL 3 ANALYSIS THEN
     * CHECK IF EVENT IS IN .evp FILE
     */
    perform_analysis_on_event = false;
    if(basic == false) {
       for(int pp = 0; pp < npassed_previous; pp++) { if(ii == passed_event[pp]) { perform_analysis_on_event = true; } }
    }
    if(!perform_analysis_on_event && basic == false) { continue; }

    /* GRAB EVENT ENTRY
     * FROM ROOT FILE
     * AND PRINT EVENT NUMBER
     */
    t.GetEntry(ii);

    if(ii%1 == 0) { cout << "Event number: " << ii << "\r" << flush; }

    /*
     * Print optional weights for debugging
     */
    /*for(size_t ww = 0; ww < theOptWeightsNames->size(); ww++) {
      std::string name = (*theOptWeightsNames)[ww];
      std::cout << "theOptWeightsNames = " << name.c_str() << ":\t " << (*theOptWeights)[ww] << endl;
      }*/

    // Print IDs of incoming partons
    /*cout << endl;
      cout << incoming[4][0] << "\t" << incoming[4][1] << endl;*/


     fastjet::PseudoJet ETmiss = fastjet::PseudoJet(theETmiss[1], theETmiss[2], theETmiss[3], theETmiss[0]);

    /*
     * FILL IN THE HISTOGRAMS
     */
     //find the leptons first:
     std::vector<fastjet::PseudoJet> leptons_array;
     for(int jj = 0; jj < numparticles; jj++) {
       fastjet::PseudoJet fourvector = fastjet::PseudoJet(objects[1][jj], objects[2][jj], objects[3][jj], objects[0][jj]);
       if(fabs(objects[4][jj]) == 11 || fabs(objects[4][jj]) == 13) { // 11 or 13 for electrons or muons
	 leptons_array.push_back(fourvector);
	 h_pT_leptons.thfill(fourvector.perp());
       }
     }
     h_ETmiss.thfill(ETmiss.perp());
     if(leptons_array.size() >= 4) {
       double m4l = (leptons_array[0] + leptons_array[1] + leptons_array[2] + leptons_array[3]).m();
       h_m4l.thfill(m4l);
     }

     h_ETmiss.thfill(ETmiss.perp());

    /*
     * APPLY CUTS
     */

    /*
     * DOES THE EVENT PASS ALL THE CUTS?
     * IF SO INCREMENT THE WEIGHT
     */
    passcuts+=evweight;

    /*
     * Fill in the _var.root file for further analysis.
     */

    //    variables[0] = ;
    eventweight[0] = evweight;

    Data2->Fill();

    /* IF EVENT HAS PASSED CUTS
     * PRINT TO .evp or .evp2 FILE
     * INCREMENT AND CONTINUE
     */
    outevp << ii << endl;

  } /* LOOP OVER EVENTS ENDS HERE
     * ENDS HERE
     */
  Data2->GetCurrentFile();
  Data2->Write();
  dat2->Close();
  cout << "A root tree has been written to the file: " << fnameroot << endl;

  /* OUTPUT HISTOGRAMS
   * HERE AND
   * FINISH
   */
  h_dummy.plot(output,1,0);
  h_pT_leptons.add(output,1,0);
  h_m4l.add(output,1,0);
  h_ETmiss.add(output,1,0);


  cout << "------------------" << endl;
  cout << "passed = " << passcuts << endl;
  cout << "------------------" << endl;

  return 0;
}

double dot(fastjet::PseudoJet p1, fastjet::PseudoJet p2) {
  return (p1.e() * p2.e() - p1.px() * p2.px() - p1.py() * p2.py() - p1.pz() * p2.pz() );
}


double deltaR(fastjet::PseudoJet p1, fastjet::PseudoJet p2) {
  double dphi_tmp;

  dphi_tmp = p2.phi() - p1.phi();
  if(dphi_tmp > M_PI)
    dphi_tmp = 2 * M_PI - dphi_tmp;
  else if( dphi_tmp < - M_PI)
    dphi_tmp = 2 * M_PI + dphi_tmp;

  //  return sqrt(sqr(p1.eta() - p2.eta()) + sqr(dphi_tmp));
  return sqrt(sqr(p1.rap() - p2.rap()) + sqr(dphi_tmp));
}

//----------------------------------------------------------------------
// does the actual work for printing out a jet
//----------------------------------------------------------------------
ostream & operator<<(ostream & ostr, const PseudoJet & jet) {
  ostr << "e, pt, y, phi ="
       << " " << setw(10) <<  jet.e()
       << " " << setw(10) << jet.perp()
       << " " << setw(6) <<  jet.rap()
       << " " << setw(6) <<  jet.phi()
       << ", mass = " << setw(10) << jet.m()
       << ", btag = " << jet.user_index();
  return ostr;
}
char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


bool btag_hadrons(fastjet::PseudoJet jet) {
  bool btagged(false);
  /* search constintuents of jets for b-mesons */
  for(int cc = 0; cc < jet.constituents().size(); cc++) {
    for(int bb = 0; bb < 105; bb++) {
      if(jet.constituents()[cc].user_index() == bhadronid[bb]) {
	btagged = true;
	//	cout << "Jet B-tagged!" << endl;
	//	cout << jet << endl;
      }
    }
  }
  return btagged;
}


fastjet::PseudoJet smear_jet(fastjet::PseudoJet jet_in) {
  if(donotsmear_jets) { return jet_in; }

  fastjet::PseudoJet smeared;
  double smearing = 20, smeared_pt(0);

  double pt = jet_in.perp();
  double eta = fabs(jet_in.eta());
  double sigma(0);

  double a, b, S, C;
  if(eta < 0.8) { a = 3.2; b = 0.07; S = 0.74; C = 0.05; }
  if(eta > 0.8 && eta < 1.2) { a = 3.0; b = 0.07; S = 0.81; C = 0.05; }
  if(eta > 1.2 && eta < 2.8) { a = 3.3; b = 0.08; S = 0.54; C = 0.05; }
  if(eta > 2.8 /*&& eta < 3.6*/) { a = 2.8; b = 0.11; S = 0.83; C = 0.05; }

  double mu_pileup = 40;
  double N = a + b * mu_pileup;

  sigma = pt * sqrt( sqr(N)/sqr(pt) + sqr(S) / pt + sqr(C) );

  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;

  double deltaE = - jet_in.e() + sqrt( sqr(jet_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*jet_in.px() + smeared_pt*sin(theta)*sin(phi)*jet_in.py() + smeared_pt*cos(theta)*jet_in.pz()));

  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);


  smeared = jet_in + smearing_vector;

  return smeared;
}

fastjet::PseudoJet smear_photon(fastjet::PseudoJet photon_in) {
  if(donotsmear_photons) { return photon_in; }

  fastjet::PseudoJet smeared;
  double smeared_pt = 0;
  double smear_frac = 0.1E-2;
  double smear_sampling = 0.15;
  double sigma(smear_sampling * sqrt(photon_in.perp()) + smear_frac*photon_in.perp());

  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;

  double deltaE = - photon_in.e() + sqrt( sqr(photon_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*photon_in.px() + smeared_pt*sin(theta)*sin(phi)*photon_in.py() + smeared_pt*cos(theta)*photon_in.pz()));

  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);
  smeared = photon_in + smearing_vector;
  //cout << "smeared mass = " << smeared.m() << endl;
  return smeared;

}



fastjet::PseudoJet smear_lepton(fastjet::PseudoJet lepton_in, int lepton_id) {

  if(donotsmear_leptons) { return lepton_in; }


  fastjet::PseudoJet smeared;
  double smeared_pt = 0;
  double smearing = 20.;

  double pt = lepton_in.perp();
  double lepton_energy = lepton_in.e();
  double eta = fabs(lepton_in.eta());
  double sigma(0);

  //see ATL-PHYS-PUB-2013-009
  if(lepton_id == 13) {
    double sigma_id = 0;
    double sigma_ms = 0;
    double sigma_cb = 0;
    double a1, a2, b0, b1, b2;

    if(eta < 0.18) { a1 = 0.01061; a2 = 0.000157; }
    if(eta > 0.18 && eta < 0.36) { a1 = 0.01084; a2 = 0.000153; }
    if(eta > 0.36 && eta < 0.54) { a1 = 0.01124; a2 = 0.000150; }
    if(eta > 0.54 && eta < 0.72) { a1 = 0.01173; a2 = 0.000149; }
    if(eta > 0.72 && eta < 0.90) { a1 = 0.01269; a2 = 0.000148; }
    if(eta > 0.90 && eta < 1.08) { a1 = 0.01406; a2 = 0.000161; }
    if(eta > 1.08 && eta < 1.26) { a1 = 0.01623; a2 = 0.000192; }
    if(eta > 1.26 && eta < 1.44) { a1 = 0.01755; a2 = 0.000199; }
    if(eta > 1.44 && eta < 1.62) { a1 = 0.01997; a2 = 0.000232; }
    if(eta > 1.62 && eta < 1.80) { a1 = 0.02453; a2 = 0.000261; }
    if(eta > 1.80 && eta < 1.98) { a1 = 0.03121; a2 = 0.000297; }
    if(eta > 1.98 && eta < 2.16) { a1 = 0.03858; a2 = 0.000375; }
    if(eta > 2.16 && eta < 2.34) { a1 = 0.05273; a2 = 0.000465; }
    if(eta > 2.34 && eta < 2.52) { a1 = 0.05329; a2 = 0.000642; }
    if(eta > 2.52 /*&& eta < 2.70*/) { a1 = 0.05683; a2 = 0.000746; }

    if(eta < 1.05) { b1 = 0.02676; b2 = 0.00012; }
    if(eta > 1.05) { b1 = 0.03880; b2 = 0.00016; }

    sigma_id = pt * sqrt( a1 + sqr(a2 * pt) );
    sigma_ms = pt * sqrt( sqr(b0/pt) + sqr(b1) + sqr(b2*pt) );
    sigma = (sigma_id * sigma_ms)/sqrt( sqr(sigma_id) + sqr(sigma_ms) ); //sigma_cb

  }


  if(lepton_id == 11) {
    double sigma = 0;
    if(eta < 1.4) { sigma = sqrt( sqr(0.3) + sqr(0.10 * sqrt(lepton_energy)) + sqr( 0.010 * lepton_energy ) ); }
    if(eta > 1.4 /* && eta < 2.47 */) { sigma = sqrt( sqr(0.3) + sqr(0.15 * sqrt(lepton_energy)) + sqr( 0.015 * lepton_energy ) ); }
  }

  smeared_pt = fabs(rnd.Gaus(0,sigma));
  double theta = rnd.Rndm()*M_PI;
  double phi = rnd.Rndm()*2.*M_PI;


  double deltaE = - lepton_in.e() + sqrt( sqr(lepton_in.e()) + sqr(smeared_pt) + 2 * (smeared_pt*sin(theta)*cos(phi)*lepton_in.px() + smeared_pt*sin(theta)*sin(phi)*lepton_in.py() + smeared_pt*cos(theta)*lepton_in.pz()));

  fastjet::PseudoJet smearing_vector(smeared_pt*sin(theta)*cos(phi),smeared_pt*sin(theta)*sin(phi), smeared_pt*cos(theta), deltaE);

  smeared = lepton_in + smearing_vector;

  //smeared.reset(smeared.px(), smeared.py(), smeared.pz(), eprime);

  return smeared;
}

bool lepton_efficiency_accept(fastjet::PseudoJet lepton_in, int lepton_id) {
  bool accepted(1);
  if(donot_apply_efficiency) { return accepted; }

  double pt = lepton_in.perp();
  double eta = fabs(lepton_in.eta());

  double epsilon = 0;
  if(lepton_id == 11) {
    epsilon = 0.85 - 0.191 * exp(1 - pt/20);
  }
  if(lepton_id == 13) {
    if(eta<0.1) { epsilon = 0.54; }
    if(eta>0.1) { epsilon = 0.97; }
  }
  double random_num = rnd.Rndm();
  //  cout << lepton_id << " " << pt << " " << eta << " " << random_num << " " << epsilon << endl;
  if(random_num > epsilon) { accepted = 0; }
  return accepted;
}
bool photon_efficiency_accept(fastjet::PseudoJet photon_in) {
  bool accepted(1);
  if(donot_apply_efficiency) { return accepted; }

  double pt = photon_in.perp();
  double eta = fabs(photon_in.eta());

  double epsilon = 0;

  epsilon = 0.76 - 1.98 * exp(-pt/16.1);

  double random_num = rnd.Rndm();
  if(random_num > epsilon) { accepted = 0; }
  return accepted;
}
bool jet_efficiency_accept(fastjet::PseudoJet jet_in) {
    bool accepted(1);
    if(donot_apply_efficiency) { return accepted; }

    double pt = jet_in.perp();
    double epsilon = 0;

    epsilon = 0.75 + (0.95 - 0.75) * pt / (50. - 20.);
    if(epsilon < 0) { epsilon = 0; }
    if(epsilon > 1.0) { epsilon = 1.0; }


    double random_num = rnd.Rndm();
    if(random_num > epsilon) { accepted = 0; }
    return accepted;

}
