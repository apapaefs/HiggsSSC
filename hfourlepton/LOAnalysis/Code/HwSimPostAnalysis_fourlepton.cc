#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

namespace {

constexpr int kMaxParticles = 10000;
constexpr int kMaxRecoObjects = 100;
constexpr int kMaxEnumeratedLeptons = 10;
constexpr int kLeptonsPerQuadruplet = 4;
constexpr int kCutStageCount = 10;

constexpr double kPi = 3.14159265358979323846;
constexpr double kMuonMassGeV = 0.1056583755;
constexpr double kZMassGeV = 91.1876;

constexpr double kLeptonPtMinGeV = 10.0;
constexpr double kLeptonEtaMax = 2.5;
constexpr double kElectronCrackEtaMin = 1.01;
constexpr double kElectronCrackEtaMax = 1.16;
constexpr double kDressingDeltaR = 0.10;
constexpr double kIsolationDeltaR = 0.35;
constexpr double kIsolationEtMaxGeV = 5.0;
constexpr double kIsolationConstituentEtMinGeV = 0.5;
constexpr double kAllOssfMassMinGeV = 4.0;
constexpr double kZ1MassMinGeV = 70.0;
constexpr double kZ1MassMaxGeV = 100.0;
constexpr double kZ2MassMinGeV = 10.0;
constexpr double kZ2MassMaxGeV = 100.0;
constexpr double kJetPtMinGeV = 20.0;
constexpr double kJetEtaMax = 5.5;

// GEM electromagnetic calorimeter response. Energies are in GeV.
constexpr double kEmBarrelSampling = 0.060;
constexpr double kEmEndcapSampling = 0.085;
constexpr double kEmConstant = 0.004;
constexpr double kEmBarrelThermalNoiseEt = 0.100;
constexpr double kEmEndcapThermalNoiseEt = 0.175;
constexpr double kEmCentralPileupNoiseEt = 0.120;

constexpr double kElectronEfficiency = 0.90;
constexpr double kMuonEfficiency = 0.85 * 0.95;
constexpr double kFourElectronTriggerEfficiency = 0.98;
constexpr double kOtherFourLeptonTriggerEfficiency = 0.99;

enum CutBit : std::uint32_t {
  kAtLeastFour = 1U << 0U,
  kKinematicAcceptance = 1U << 1U,
  kElectronCrackAcceptance = 1U << 2U,
  kIsolation = 1U << 3U,
  kOssfPairing = 1U << 4U,
  kAllOssfMasses = 1U << 5U,
  kZ1Window = 1U << 6U,
  kZ2Window = 1U << 7U,
  kTrigger = 1U << 8U,
  kSelected = 1U << 9U,
};

constexpr std::uint32_t kSelectedMask =
    kAtLeastFour | kKinematicAcceptance | kElectronCrackAcceptance | kIsolation |
    kOssfPairing | kAllOssfMasses | kZ1Window | kZ2Window | kTrigger | kSelected;

enum class ResponseProfile {
  Perfect,
  Ssc,
};

struct FourVector {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e = 0.0;

  double pt() const { return std::hypot(px, py); }

  double pabs() const { return std::sqrt(px * px + py * py + pz * pz); }

  double transverse_energy() const {
    const double momentum = pabs();
    return momentum > 0.0 ? e * pt() / momentum : 0.0;
  }

  double eta() const {
    const double momentum = pabs();
    const double plus = momentum + pz;
    const double minus = momentum - pz;
    if (plus <= 0.0 || minus <= 0.0) {
      return pz >= 0.0 ? 1.0e9 : -1.0e9;
    }
    return 0.5 * std::log(plus / minus);
  }

  double rapidity() const {
    const double plus = e + pz;
    const double minus = e - pz;
    if (plus <= 0.0 || minus <= 0.0) {
      return pz >= 0.0 ? 1.0e9 : -1.0e9;
    }
    return 0.5 * std::log(plus / minus);
  }

  double phi() const { return std::atan2(py, px); }

  double mass() const {
    const double mass_squared = e * e - px * px - py * py - pz * pz;
    return std::sqrt(std::max(0.0, mass_squared));
  }
};

FourVector operator+(const FourVector& left, const FourVector& right) {
  return {left.px + right.px, left.py + right.py, left.pz + right.pz, left.e + right.e};
}

FourVector operator-(const FourVector& left, const FourVector& right) {
  return {left.px - right.px, left.py - right.py, left.pz - right.pz,
          left.e - right.e};
}

struct Particle {
  FourVector p4;
  int pdg_id = 0;
  int source_index = -1;
};

struct Lepton {
  FourVector bare;
  FourVector dressed;
  FourVector reconstructed;
  int pdg_id = 0;
  int flavour = 0;
  int charge = 0;
  int source_index = -1;
  std::vector<int> dressed_photon_indices;
  double isolation_et = 0.0;
  double efficiency = 1.0;
};

struct RecoContext {
  int n_jets = -1;
  int n_bjets = -1;
  double leading_jet_pt = -1.0;
  double met_pt = -1.0;
  double met_phi = 0.0;
};

struct PairingChoice {
  bool has_pairing = false;
  std::array<int, 4> lepton_indices = {{-1, -1, -1, -1}};
  FourVector z1;
  FourVector z2;
  FourVector four_lepton;
  double scalar_pt_sum = 0.0;
  std::array<int, 4> sorted_source_indices = {{-1, -1, -1, -1}};
  std::array<int, 4> pairing_source_indices = {{-1, -1, -1, -1}};
};

struct SubsetResult {
  std::uint32_t cut_mask = 0U;
  int channel = 0;
  std::array<int, 4> output_indices = {{-1, -1, -1, -1}};
  FourVector four_lepton;
  FourVector z1;
  FourVector z2;
  double min_delta_r = -1.0;
  double z1_delta_r = -1.0;
  double z2_delta_r = -1.0;
};

struct Config {
  std::string input_name;
  std::string output_dir;
  std::string tag;
  std::string sample;
  std::string category;
  std::string requested_channel;
  std::string optional_weight_names_file;
  std::string response_profile_name = "perfect";
  ResponseProfile response_profile = ResponseProfile::Perfect;
  std::uint64_t seed = 14101983ULL;
  double weight_scale = 1.0;
  double muon_resolution_scale = 1.0;
  double electron_efficiency = kElectronEfficiency;
  double muon_efficiency = kMuonEfficiency;
  double four_electron_trigger_efficiency = kFourElectronTriggerEfficiency;
  double other_trigger_efficiency = kOtherFourLeptonTriggerEfficiency;
  bool include_pileup_noise = true;
  long long first_event = 0;
  long long last_event = -1;
};

struct Diagnostics {
  long long invalid_particle_records = 0;
  long long lepton_overflow_events = 0;
  long long marginalized_nonbaseline_leptons = 0;
  long long muon_pt_underflow = 0;
  long long muon_pt_overflow = 0;
  long long optional_weight_name_syntheses = 0;
  long long optional_weight_names_from_input = 0;
  long long optional_weight_names_from_file = 0;
  long long optional_weight_name_file_validations = 0;
  long long optional_weight_name_mismatches = 0;
  long long closure_failures = 0;
  double maximum_closure_delta = 0.0;
};

struct WeightSums {
  long long entries = 0;
  long long positive_entries = 0;
  long long negative_entries = 0;
  long long zero_entries = 0;
  double sumw = 0.0;
  double sumabsw = 0.0;
  double sumw2 = 0.0;
  double sumw_positive = 0.0;
  double sumw_negative = 0.0;

  void add(double weight) {
    ++entries;
    sumw += weight;
    sumabsw += std::fabs(weight);
    sumw2 += weight * weight;
    if (weight > 0.0) {
      ++positive_entries;
      sumw_positive += weight;
    } else if (weight < 0.0) {
      ++negative_entries;
      sumw_negative += weight;
    } else {
      ++zero_entries;
    }
  }

  double negative_fraction() const {
    return entries == 0 ? 0.0
                        : static_cast<double>(negative_entries) / static_cast<double>(entries);
  }
};

struct SourceRecord {
  Long64_t source_index = -1;
  double generator_weight = 0.0;
  double sample_weight = 1.0;
  double base_event_weight = 0.0;
  double response_failure_probability = 1.0;
  double response_probability_in_tree = 0.0;
  double selected_probability = 0.0;
  double selected_event_weight = 0.0;
  double response_closure_delta = 0.0;
  int raw_lepton_count = 0;
  int accepted_lepton_count = 0;
  int subset_hypothesis_count = 0;
  double cut_probabilities[kCutStageCount] = {};
  std::vector<std::string> optional_weight_names;
  std::vector<double> optional_weights;
  std::vector<double> optional_selected_event_weights;
};

struct FourLeptonRecord {
  Long64_t source_index = -1;
  ULong64_t subset_mask = 0;
  UInt_t cut_mask = 0U;
  int channel = 0;
  int candidate_count = 0;
  double hypothesis_probability = 0.0;
  double generator_weight = 0.0;
  double sample_weight = 1.0;
  double response_weight = 0.0;
  double trigger_weight = 1.0;
  double pretrigger_event_weight = 0.0;
  double event_weight = 0.0;
  std::vector<double> optional_event_weights;
  double m4l = -1.0;
  double mZ1 = -1.0;
  double mZ2 = -1.0;
  double pt4l = -1.0;
  double y4l = 0.0;
  double min_delta_r = -1.0;
  double z1_delta_r = -1.0;
  double z2_delta_r = -1.0;
  int n_jets = -1;
  int n_bjets = -1;
  double leading_jet_pt = -1.0;
  double met_pt = -1.0;
  double met_phi = 0.0;
  double lepton_pt[kLeptonsPerQuadruplet] = {};
  double lepton_eta[kLeptonsPerQuadruplet] = {};
  double lepton_phi[kLeptonsPerQuadruplet] = {};
  double lepton_isolation[kLeptonsPerQuadruplet] = {};
  int lepton_flavour[kLeptonsPerQuadruplet] = {};
  int lepton_charge[kLeptonsPerQuadruplet] = {};
  int lepton_source_index[kLeptonsPerQuadruplet] = {};
  int lepton_dressed_photon_count[kLeptonsPerQuadruplet] = {};
};

std::string cut_mask_definition() {
  return "bit0:at_least_four_reconstructed;"
         "bit1:pt_eta_acceptance;"
         "bit2:electron_crack_acceptance;"
         "bit3:isolation;"
         "bit4:disjoint_ossf_pairing;"
         "bit5:all_ossf_masses_gt_4;"
         "bit6:z1_70_100;"
         "bit7:z2_10_100;"
         "bit8:trigger_weight_applied;"
         "bit9:selected";
}

std::string cut_bit_zero_semantics() {
  return "at least four reconstructed source electrons or muons before "
         "pT, eta, electron-crack, and isolation requirements; this is not "
         "the baseline-accepted multiplicity";
}

void print_usage(std::ostream& output) {
  output
      << "Usage:\n"
      << "  HwSimPostAnalysis_fourlepton INPUT.root|INPUT.input [options]\n"
      << "  HwSimPostAnalysis_fourlepton --input-list INPUT.input [options]\n\n"
      << "Options:\n"
      << "  --response-profile perfect|ssc   Detector response (default: perfect)\n"
      << "  --detector-response perfect|ssc  Alias for --response-profile\n"
      << "  --seed INTEGER                   Deterministic response seed\n"
      << "  --weight-scale FLOAT             Sample normalization factor (default: 1)\n"
      << "  --muon-resolution-scale FLOAT    Multiplier for sigma(pT)/pT (default: 1)\n"
      << "  --optional-weight-names-file FILE\n"
      << "                                    One semantic optional-weight name per line\n"
      << "  --tag TAG                        Output tag\n"
      << "  --output-dir DIRECTORY           Output directory\n"
      << "  --sample NAME                    Sample metadata and output stem\n"
      << "  --category NAME                  Category metadata\n"
      << "  --channel NAME                   Generated-channel metadata\n"
      << "  -n N | -nmin I -nmax J           Optional source-entry range\n"
      << "  --no-pileup-noise                Disable EM thermal/pileup noise term\n"
      << "  --help                            Show this message\n\n"
      << "With --output-dir/--sample/--tag, outputs are\n"
      << "  <output-dir>/<sample>_<tag>_<profile>.root\n"
      << "  <output-dir>/<sample>_<tag>_<profile>.summary.json\n";
}

std::string require_option_value(int& index, int argc, char* argv[],
                                 const std::string& option) {
  if (index + 1 >= argc) {
    throw std::runtime_error("missing value for " + option);
  }
  ++index;
  return argv[index];
}

double parse_double(const std::string& text, const std::string& option) {
  std::size_t parsed = 0;
  const double value = std::stod(text, &parsed);
  if (parsed != text.size() || !std::isfinite(value)) {
    throw std::runtime_error("invalid numeric value for " + option + ": " + text);
  }
  return value;
}

long long parse_int64(const std::string& text, const std::string& option) {
  std::size_t parsed = 0;
  const long long value = std::stoll(text, &parsed);
  if (parsed != text.size()) {
    throw std::runtime_error("invalid integer value for " + option + ": " + text);
  }
  return value;
}

std::uint64_t parse_uint64(const std::string& text, const std::string& option) {
  std::size_t parsed = 0;
  const unsigned long long value = std::stoull(text, &parsed);
  if (parsed != text.size()) {
    throw std::runtime_error("invalid integer value for " + option + ": " + text);
  }
  return static_cast<std::uint64_t>(value);
}

void require_probability(double value, const std::string& name) {
  if (!(value >= 0.0 && value <= 1.0)) {
    throw std::runtime_error(name + " must lie in [0,1]");
  }
}

ResponseProfile parse_response_profile(const std::string& text) {
  if (text == "perfect") {
    return ResponseProfile::Perfect;
  }
  if (text == "ssc") {
    return ResponseProfile::Ssc;
  }
  throw std::runtime_error("unknown response profile '" + text +
                           "' (expected perfect or ssc)");
}

Config parse_config(int argc, char* argv[]) {
  Config config;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    if (argument == "--help" || argument == "-h") {
      print_usage(std::cout);
      std::exit(0);
    } else if (argument == "--input-list" || argument == "--input") {
      config.input_name = require_option_value(index, argc, argv, argument);
    } else if (argument == "--response-profile" || argument == "--detector-response") {
      config.response_profile_name = require_option_value(index, argc, argv, argument);
    } else if (argument == "--seed") {
      config.seed = parse_uint64(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--weight-scale" || argument == "-w") {
      config.weight_scale =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--muon-resolution-scale") {
      config.muon_resolution_scale =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--optional-weight-names-file") {
      config.optional_weight_names_file =
          require_option_value(index, argc, argv, argument);
    } else if (argument == "--tag" || argument == "--run-tag" || argument == "-t") {
      config.tag = require_option_value(index, argc, argv, argument);
    } else if (argument == "--output-dir") {
      config.output_dir = require_option_value(index, argc, argv, argument);
    } else if (argument == "--sample") {
      config.sample = require_option_value(index, argc, argv, argument);
    } else if (argument == "--category") {
      config.category = require_option_value(index, argc, argv, argument);
    } else if (argument == "--channel") {
      config.requested_channel = require_option_value(index, argc, argv, argument);
    } else if (argument == "--electron-efficiency") {
      config.electron_efficiency =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--muon-efficiency") {
      config.muon_efficiency =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--trigger-4e-efficiency") {
      config.four_electron_trigger_efficiency =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--trigger-other-efficiency") {
      config.other_trigger_efficiency =
          parse_double(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "--no-pileup-noise") {
      config.include_pileup_noise = false;
    } else if (argument == "-n") {
      config.last_event =
          parse_int64(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "-nmin") {
      config.first_event =
          parse_int64(require_option_value(index, argc, argv, argument), argument);
    } else if (argument == "-nmax") {
      config.last_event =
          parse_int64(require_option_value(index, argc, argv, argument), argument);
    } else if (!argument.empty() && argument[0] == '-') {
      throw std::runtime_error("unknown option " + argument);
    } else if (config.input_name.empty()) {
      config.input_name = argument;
    } else {
      throw std::runtime_error("unexpected positional argument " + argument);
    }
  }

  if (config.input_name.empty()) {
    throw std::runtime_error("no input ROOT file or input list was specified");
  }
  config.response_profile = parse_response_profile(config.response_profile_name);
  if (config.muon_resolution_scale < 0.0) {
    throw std::runtime_error("--muon-resolution-scale must be non-negative");
  }
  require_probability(config.electron_efficiency, "electron efficiency");
  require_probability(config.muon_efficiency, "muon efficiency");
  require_probability(config.four_electron_trigger_efficiency, "4e trigger efficiency");
  require_probability(config.other_trigger_efficiency, "other-channel trigger efficiency");
  if (config.first_event < 0 || config.last_event == 0) {
    throw std::runtime_error("event-range values must be positive, except -nmin may be zero");
  }
  return config;
}

std::string trim(const std::string& text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::vector<std::string> load_optional_weight_names(
    const std::string& path) {
  if (path.empty()) {
    return {};
  }
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error(
        "failed to open optional-weight names file " + path);
  }
  std::vector<std::string> names;
  std::string line;
  while (std::getline(input, line)) {
    const std::string name = trim(line);
    if (name.empty() || name.front() == '#') {
      continue;
    }
    if (std::find(names.begin(), names.end(), name) != names.end()) {
      throw std::runtime_error(
          "duplicate optional-weight name in " + path + ": " + name);
    }
    names.push_back(name);
  }
  if (names.empty()) {
    throw std::runtime_error(
        "optional-weight names file contains no names: " + path);
  }
  return names;
}

std::string optional_weight_names_source(const Diagnostics& diagnostics,
                                         bool names_file_configured,
                                         bool optional_weights_available) {
  std::vector<std::string> sources;
  if (diagnostics.optional_weight_names_from_input != 0) {
    sources.push_back(
        names_file_configured ? "HwSim_branch_validated_by_file"
                              : "HwSim_branch");
  }
  if (diagnostics.optional_weight_names_from_file != 0) {
    sources.push_back("configured_file");
  }
  if (diagnostics.optional_weight_name_syntheses != 0) {
    sources.push_back("synthesized_positional");
  }
  if (sources.empty()) {
    return optional_weights_available ? "empty_optional_weight_vectors"
                                      : "not_available";
  }
  std::ostringstream output;
  for (std::size_t index = 0; index < sources.size(); ++index) {
    if (index != 0) {
      output << '+';
    }
    output << sources[index];
  }
  return output.str();
}

bool has_suffix(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void add_input_files(TChain& chain, const std::string& input_name) {
  if (has_suffix(input_name, ".root")) {
    if (chain.Add(input_name.c_str()) == 0) {
      throw std::runtime_error("failed to add ROOT input " + input_name);
    }
    return;
  }
  if (!has_suffix(input_name, ".input")) {
    throw std::runtime_error("input must end in .root or .input");
  }

  std::ifstream input(input_name);
  if (!input) {
    throw std::runtime_error("failed to open input list " + input_name);
  }
  std::string line;
  int files_added = 0;
  while (std::getline(input, line)) {
    const std::string path = trim(line);
    if (path.empty() || path[0] == '#') {
      continue;
    }
    if (chain.Add(path.c_str()) == 0) {
      throw std::runtime_error("failed to add ROOT input " + path);
    }
    ++files_added;
  }
  if (files_added == 0) {
    throw std::runtime_error("input list contains no ROOT files: " + input_name);
  }
}

std::string safe_filename_component(const std::string& input) {
  std::string output;
  output.reserve(input.size());
  for (const unsigned char character : input) {
    if (std::isalnum(character) || character == '_' || character == '-' ||
        character == '.') {
      output.push_back(static_cast<char>(character));
    } else {
      output.push_back('_');
    }
  }
  while (!output.empty() && output.front() == '.') {
    output.erase(output.begin());
  }
  return output.empty() ? "unnamed" : output;
}

std::pair<std::filesystem::path, std::filesystem::path> output_paths(const Config& config) {
  const std::filesystem::path input_path(config.input_name);
  const std::filesystem::path directory =
      config.output_dir.empty()
          ? (input_path.has_parent_path() ? input_path.parent_path() : std::filesystem::path("."))
          : std::filesystem::path(config.output_dir);
  std::filesystem::create_directories(directory);

  std::string stem;
  if (!config.sample.empty()) {
    stem = safe_filename_component(config.sample);
  } else {
    stem = safe_filename_component(input_path.stem().string());
  }
  if (!config.tag.empty()) {
    stem += "_" + safe_filename_component(config.tag);
  }
  stem += "_" + safe_filename_component(config.response_profile_name);
  return {directory / (stem + ".root"), directory / (stem + ".summary.json")};
}

double delta_phi(const FourVector& first, const FourVector& second) {
  double difference = second.phi() - first.phi();
  while (difference > kPi) {
    difference -= 2.0 * kPi;
  }
  while (difference <= -kPi) {
    difference += 2.0 * kPi;
  }
  return difference;
}

double delta_r(const FourVector& first, const FourVector& second) {
  return std::hypot(first.eta() - second.eta(), delta_phi(first, second));
}

FourVector from_array(const double objects[][kMaxParticles], int index) {
  return {objects[1][index], objects[2][index], objects[3][index], objects[0][index]};
}

FourVector from_reco_array(const double objects[][kMaxRecoObjects], int index) {
  return {objects[1][index], objects[2][index], objects[3][index], objects[0][index]};
}

FourVector with_pt_eta_phi_mass(double pt, double eta, double phi, double mass) {
  const double safe_pt = std::max(1.0e-9, pt);
  const double px = safe_pt * std::cos(phi);
  const double py = safe_pt * std::sin(phi);
  const double pz = safe_pt * std::sinh(eta);
  const double momentum_squared = px * px + py * py + pz * pz;
  return {px, py, pz, std::sqrt(std::max(0.0, momentum_squared + mass * mass))};
}

FourVector massless_with_energy_direction(double energy, double eta, double phi) {
  const double safe_energy = std::max(1.0e-9, energy);
  return with_pt_eta_phi_mass(safe_energy / std::cosh(eta), eta, phi, 0.0);
}

std::uint64_t splitmix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

unsigned int object_seed(std::uint64_t base_seed, std::uint64_t event_index,
                         std::uint64_t object_kind, std::uint64_t object_index) {
  std::uint64_t value = splitmix64(base_seed);
  value ^= splitmix64(event_index + 0x100000001b3ULL);
  value ^= splitmix64((object_kind << 32U) ^ object_index);
  unsigned int seed = static_cast<unsigned int>(splitmix64(value) & 0xffffffffULL);
  return seed == 0U ? 1U : seed;
}

double gaussian(double mean, double sigma, unsigned int seed) {
  if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma < 0.0) {
    throw std::runtime_error("invalid Gaussian response parameters");
  }
  if (sigma == 0.0) {
    return mean;
  }
  TRandom3 random(seed);
  return random.Gaus(mean, sigma);
}

double em_noise_et(double abs_eta, bool include_pileup_noise) {
  const double thermal = abs_eta < kElectronCrackEtaMin ? kEmBarrelThermalNoiseEt
                                                        : kEmEndcapThermalNoiseEt;
  if (!include_pileup_noise) {
    return thermal;
  }
  const double pileup =
      abs_eta < 1.4 ? kEmCentralPileupNoiseEt
                    : kEmCentralPileupNoiseEt + 0.366 * (abs_eta - 1.4);
  return std::hypot(thermal, pileup);
}

FourVector smear_electron(const FourVector& input, bool include_pileup_noise,
                          unsigned int seed) {
  const double abs_eta = std::fabs(input.eta());
  const double sampling =
      abs_eta < kElectronCrackEtaMin ? kEmBarrelSampling : kEmEndcapSampling;
  double relative_variance =
      sampling * sampling / std::max(input.e, 1.0e-9) + kEmConstant * kEmConstant;
  if (input.pt() > 0.0) {
    const double noise_fraction =
        em_noise_et(abs_eta, include_pileup_noise) / input.pt();
    relative_variance += noise_fraction * noise_fraction;
  }
  const double energy =
      std::max(1.0e-9,
               gaussian(input.e, input.e * std::sqrt(relative_variance), seed));
  return massless_with_energy_direction(energy, input.eta(), input.phi());
}

template <std::size_t N>
double interpolate_linear_clamped(double value, const std::array<double, N>& knots,
                                  const std::array<double, N>& table) {
  if (value <= knots.front()) {
    return table.front();
  }
  if (value >= knots.back()) {
    return table.back();
  }
  for (std::size_t index = 0; index + 1 < N; ++index) {
    if (value <= knots[index + 1]) {
      const double fraction =
          (value - knots[index]) / (knots[index + 1] - knots[index]);
      return table[index] + fraction * (table[index + 1] - table[index]);
    }
  }
  return table.back();
}

double muon_fractional_resolution(double pt, double abs_eta, Diagnostics& diagnostics) {
  // Approximate digitization of GEM TDR Fig. 4-20. These are fractional
  // sigma(pT)/pT values, not an exact numerical table from the TDR.
  constexpr std::array<double, 9> eta_knots =
      {{0.1, 0.6, 1.0, 1.3, 1.6, 1.85, 2.0, 2.2, 2.5}};
  constexpr std::array<double, 9> pt10 =
      {{0.028, 0.0265, 0.0255, 0.0215, 0.0190, 0.0195, 0.0200, 0.0205, 0.0220}};
  constexpr std::array<double, 9> pt25 =
      {{0.0167, 0.0155, 0.0150, 0.0175, 0.0182, 0.0195, 0.0140, 0.0135, 0.0125}};
  constexpr std::array<double, 9> pt50 =
      {{0.0148, 0.0138, 0.0138, 0.0185, 0.0198, 0.0200, 0.0145, 0.0145, 0.0135}};
  constexpr std::array<double, 9> pt100 =
      {{0.0135, 0.0125, 0.0125, 0.0195, 0.0230, 0.0250, 0.0205, 0.0235, 0.0245}};
  constexpr std::array<double, 5> pt_knots = {{10.0, 25.0, 50.0, 100.0, 500.0}};

  if (pt < pt_knots.front()) {
    ++diagnostics.muon_pt_underflow;
  } else if (pt > pt_knots.back()) {
    ++diagnostics.muon_pt_overflow;
  }
  const double clamped_pt = std::max(pt_knots.front(), std::min(pt, pt_knots.back()));
  const double clamped_eta = std::max(0.0, std::min(abs_eta, 2.5));
  const std::array<double, 5> eta_values = {{
      interpolate_linear_clamped(clamped_eta, eta_knots, pt10),
      interpolate_linear_clamped(clamped_eta, eta_knots, pt25),
      interpolate_linear_clamped(clamped_eta, eta_knots, pt50),
      interpolate_linear_clamped(clamped_eta, eta_knots, pt100),
      0.05 + (0.12 - 0.05) * clamped_eta / 2.5,
  }};

  if (clamped_pt <= pt_knots.front()) {
    return eta_values.front();
  }
  if (clamped_pt >= pt_knots.back()) {
    return eta_values.back();
  }
  const double log_pt = std::log(clamped_pt);
  for (std::size_t index = 0; index + 1 < pt_knots.size(); ++index) {
    if (clamped_pt <= pt_knots[index + 1]) {
      const double fraction =
          (log_pt - std::log(pt_knots[index])) /
          (std::log(pt_knots[index + 1]) - std::log(pt_knots[index]));
      return eta_values[index] +
             fraction * (eta_values[index + 1] - eta_values[index]);
    }
  }
  return eta_values.back();
}

FourVector smear_muon(const FourVector& input, int charge, double resolution_scale,
                      unsigned int seed, Diagnostics& diagnostics) {
  const double pt = std::max(input.pt(), 1.0e-9);
  const double fractional_resolution =
      muon_fractional_resolution(pt, std::fabs(input.eta()), diagnostics) *
      resolution_scale;
  const double nominal_curvature = static_cast<double>(charge) / pt;
  const double smeared_curvature =
      gaussian(nominal_curvature, fractional_resolution / pt, seed);
  // Charge misidentification is deliberately outside the detector model.
  const double smeared_pt =
      std::min(1.0e7, 1.0 / std::max(std::fabs(smeared_curvature), 1.0e-7));
  return with_pt_eta_phi_mass(smeared_pt, input.eta(), input.phi(),
                              kMuonMassGeV);
}

bool finite_four_vector(const FourVector& vector) {
  return std::isfinite(vector.px) && std::isfinite(vector.py) &&
         std::isfinite(vector.pz) && std::isfinite(vector.e) && vector.e > 0.0;
}

bool is_neutrino(int pdg_id) {
  const int absolute_id = std::abs(pdg_id);
  return absolute_id == 12 || absolute_id == 14 || absolute_id == 16 ||
         absolute_id == 18;
}

bool is_ossf(const Lepton& first, const Lepton& second) {
  return first.flavour == second.flavour && first.charge * second.charge == -1;
}

std::vector<Lepton> build_dressed_leptons(const std::vector<Particle>& particles) {
  std::vector<Lepton> leptons;
  std::vector<const Particle*> photons;
  for (const Particle& particle : particles) {
    const int absolute_id = std::abs(particle.pdg_id);
    if (absolute_id == 11 || absolute_id == 13) {
      Lepton lepton;
      lepton.bare = particle.p4;
      lepton.dressed = particle.p4;
      lepton.reconstructed = particle.p4;
      lepton.pdg_id = particle.pdg_id;
      lepton.flavour = absolute_id;
      lepton.charge = particle.pdg_id > 0 ? -1 : 1;
      lepton.source_index = particle.source_index;
      leptons.push_back(lepton);
    } else if (particle.pdg_id == 22) {
      photons.push_back(&particle);
    }
  }

  for (const Particle* photon : photons) {
    int nearest = -1;
    double nearest_delta_r = kDressingDeltaR;
    for (std::size_t index = 0; index < leptons.size(); ++index) {
      const double distance = delta_r(leptons[index].bare, photon->p4);
      if (distance < nearest_delta_r ||
          (distance == nearest_delta_r && nearest >= 0 &&
           leptons[index].source_index < leptons[nearest].source_index)) {
        nearest = static_cast<int>(index);
        nearest_delta_r = distance;
      }
    }
    if (nearest >= 0) {
      leptons[nearest].dressed = leptons[nearest].dressed + photon->p4;
      leptons[nearest].dressed_photon_indices.push_back(photon->source_index);
    }
  }
  return leptons;
}

double isolation_et(const Lepton& candidate, const std::vector<Particle>& particles) {
  double scalar_et = 0.0;
  for (const Particle& particle : particles) {
    if (particle.source_index == candidate.source_index || is_neutrino(particle.pdg_id) ||
        particle.p4.transverse_energy() <= kIsolationConstituentEtMinGeV) {
      continue;
    }
    if (std::find(candidate.dressed_photon_indices.begin(),
                  candidate.dressed_photon_indices.end(),
                  particle.source_index) != candidate.dressed_photon_indices.end()) {
      continue;
    }
    if (delta_r(candidate.reconstructed, particle.p4) < kIsolationDeltaR) {
      scalar_et += particle.p4.transverse_energy();
    }
  }
  return scalar_et;
}

bool passes_kinematic_acceptance(const Lepton& lepton) {
  return lepton.reconstructed.pt() > kLeptonPtMinGeV &&
         std::fabs(lepton.reconstructed.eta()) < kLeptonEtaMax;
}

bool passes_electron_crack_acceptance(const Lepton& lepton) {
  const double abs_eta = std::fabs(lepton.reconstructed.eta());
  return lepton.flavour != 11 ||
         !(abs_eta > kElectronCrackEtaMin &&
           abs_eta < kElectronCrackEtaMax);
}

bool passes_isolation(const Lepton& lepton) {
  return lepton.isolation_et < kIsolationEtMaxGeV;
}

bool passes_baseline(const Lepton& lepton) {
  return passes_kinematic_acceptance(lepton) &&
         passes_electron_crack_acceptance(lepton) &&
         passes_isolation(lepton);
}

std::vector<Lepton> response_leptons(const std::vector<Particle>& particles,
                                     const Config& config, Long64_t event_index,
                                     Diagnostics& diagnostics) {
  std::vector<Lepton> leptons = build_dressed_leptons(particles);
  for (Lepton& lepton : leptons) {
    if (config.response_profile == ResponseProfile::Ssc) {
      if (lepton.flavour == 11) {
        lepton.reconstructed =
            smear_electron(lepton.dressed, config.include_pileup_noise,
                           object_seed(config.seed, event_index, 11U,
                                       static_cast<std::uint64_t>(lepton.source_index)));
        lepton.efficiency = config.electron_efficiency;
      } else {
        const FourVector dressing_photons = lepton.dressed - lepton.bare;
        lepton.reconstructed =
            smear_muon(
                lepton.bare, lepton.charge, config.muon_resolution_scale,
                object_seed(
                    config.seed, event_index, 13U,
                    static_cast<std::uint64_t>(lepton.source_index)),
                diagnostics) +
            dressing_photons;
        lepton.efficiency = config.muon_efficiency;
      }
    } else {
      lepton.reconstructed = lepton.dressed;
      lepton.efficiency = 1.0;
    }

    lepton.isolation_et = isolation_et(lepton, particles);
  }
  std::sort(leptons.begin(), leptons.end(), [](const Lepton& left, const Lepton& right) {
    return left.source_index < right.source_index;
  });
  return leptons;
}

std::array<int, 2> ordered_pair(int first, int second,
                                const std::vector<Lepton>& leptons) {
  if (leptons[first].charge < leptons[second].charge) {
    return {{first, second}};
  }
  if (leptons[second].charge < leptons[first].charge) {
    return {{second, first}};
  }
  return leptons[first].source_index < leptons[second].source_index
             ? std::array<int, 2>{{first, second}}
             : std::array<int, 2>{{second, first}};
}

bool pairing_is_better(const PairingChoice& candidate, const PairingChoice& best) {
  if (!best.has_pairing) {
    return true;
  }
  const double candidate_z_distance = std::fabs(candidate.z1.mass() - kZMassGeV);
  const double best_z_distance = std::fabs(best.z1.mass() - kZMassGeV);
  if (std::fabs(candidate_z_distance - best_z_distance) > 1.0e-12) {
    return candidate_z_distance < best_z_distance;
  }
  if (std::fabs(candidate.scalar_pt_sum - best.scalar_pt_sum) > 1.0e-12) {
    return candidate.scalar_pt_sum > best.scalar_pt_sum;
  }
  if (candidate.sorted_source_indices != best.sorted_source_indices) {
    return candidate.sorted_source_indices < best.sorted_source_indices;
  }
  return candidate.pairing_source_indices < best.pairing_source_indices;
}

PairingChoice best_pairing(const std::vector<Lepton>& leptons,
                           const std::vector<int>& subset) {
  PairingChoice best;
  const int subset_size = static_cast<int>(subset.size());
  for (int a = 0; a < subset_size; ++a) {
    for (int b = a + 1; b < subset_size; ++b) {
      for (int c = b + 1; c < subset_size; ++c) {
        for (int d = c + 1; d < subset_size; ++d) {
          const std::array<int, 4> quad =
              {{subset[a], subset[b], subset[c], subset[d]}};
          constexpr std::array<std::array<int, 4>, 3> pairings = {{
              {{0, 1, 2, 3}},
              {{0, 2, 1, 3}},
              {{0, 3, 1, 2}},
          }};
          for (const std::array<int, 4>& pattern : pairings) {
            const int first_a = quad[pattern[0]];
            const int first_b = quad[pattern[1]];
            const int second_a = quad[pattern[2]];
            const int second_b = quad[pattern[3]];
            if (!is_ossf(leptons[first_a], leptons[first_b]) ||
                !is_ossf(leptons[second_a], leptons[second_b])) {
              continue;
            }

            std::array<int, 2> pair_one = ordered_pair(first_a, first_b, leptons);
            std::array<int, 2> pair_two = ordered_pair(second_a, second_b, leptons);
            FourVector z_one = leptons[pair_one[0]].reconstructed +
                               leptons[pair_one[1]].reconstructed;
            FourVector z_two = leptons[pair_two[0]].reconstructed +
                               leptons[pair_two[1]].reconstructed;
            const double distance_one = std::fabs(z_one.mass() - kZMassGeV);
            const double distance_two = std::fabs(z_two.mass() - kZMassGeV);
            if (distance_two < distance_one ||
                (std::fabs(distance_one - distance_two) <= 1.0e-12 &&
                 std::array<int, 2>{{leptons[pair_two[0]].source_index,
                                     leptons[pair_two[1]].source_index}} <
                     std::array<int, 2>{{leptons[pair_one[0]].source_index,
                                         leptons[pair_one[1]].source_index}})) {
              std::swap(pair_one, pair_two);
              std::swap(z_one, z_two);
            }

            PairingChoice candidate;
            candidate.has_pairing = true;
            candidate.lepton_indices =
                {{pair_one[0], pair_one[1], pair_two[0], pair_two[1]}};
            candidate.z1 = z_one;
            candidate.z2 = z_two;
            candidate.four_lepton = z_one + z_two;
            candidate.scalar_pt_sum =
                std::accumulate(quad.begin(), quad.end(), 0.0,
                                [&](double sum, int index) {
                                  return sum + leptons[index].reconstructed.pt();
                                });
            for (int index = 0; index < 4; ++index) {
              candidate.sorted_source_indices[index] =
                  leptons[quad[index]].source_index;
              candidate.pairing_source_indices[index] =
                  leptons[candidate.lepton_indices[index]].source_index;
            }
            std::sort(candidate.sorted_source_indices.begin(),
                      candidate.sorted_source_indices.end());
            if (pairing_is_better(candidate, best)) {
              best = candidate;
            }
          }
        }
      }
    }
  }
  return best;
}

bool all_ossf_masses_pass(const std::array<int, 4>& indices,
                          const std::vector<Lepton>& leptons) {
  for (int first = 0; first < 4; ++first) {
    for (int second = first + 1; second < 4; ++second) {
      const Lepton& left = leptons[indices[first]];
      const Lepton& right = leptons[indices[second]];
      if (is_ossf(left, right) &&
          !((left.reconstructed + right.reconstructed).mass() >
            kAllOssfMassMinGeV)) {
        return false;
      }
    }
  }
  return true;
}

int classify_channel(const std::array<int, 4>& indices,
                     const std::vector<Lepton>& leptons) {
  int electron_count = 0;
  for (const int index : indices) {
    electron_count += leptons[index].flavour == 11 ? 1 : 0;
  }
  if (electron_count == 4) {
    return 1;
  }
  if (electron_count == 0) {
    return 2;
  }
  if (electron_count == 2) {
    return 3;
  }
  return 0;
}

SubsetResult analyze_subset(const std::vector<Lepton>& leptons,
                            const std::vector<int>& subset) {
  SubsetResult result;
  if (subset.size() < 4) {
    return result;
  }

  const auto fill_leading_quadruplet = [&](const std::vector<int>& candidates) {
    std::vector<int> ordered = candidates;
    std::stable_sort(ordered.begin(), ordered.end(),
                     [&](int left, int right) {
                       if (leptons[left].reconstructed.pt() !=
                           leptons[right].reconstructed.pt()) {
                         return leptons[left].reconstructed.pt() >
                                leptons[right].reconstructed.pt();
                       }
                       return leptons[left].source_index <
                              leptons[right].source_index;
                     });
    for (int index = 0; index < 4; ++index) {
      result.output_indices[index] = ordered[index];
      result.four_lepton =
          result.four_lepton + leptons[ordered[index]].reconstructed;
    }
    result.channel = classify_channel(result.output_indices, leptons);
  };

  result.cut_mask = kAtLeastFour;
  std::vector<int> kinematic;
  std::copy_if(subset.begin(), subset.end(), std::back_inserter(kinematic),
               [&](int index) {
                 return passes_kinematic_acceptance(leptons[index]);
               });
  if (kinematic.size() < 4) {
    fill_leading_quadruplet(subset);
  } else {
    result.cut_mask |= kKinematicAcceptance;
    std::vector<int> crack_accepted;
    std::copy_if(kinematic.begin(), kinematic.end(),
                 std::back_inserter(crack_accepted), [&](int index) {
                   return passes_electron_crack_acceptance(leptons[index]);
                 });
    if (crack_accepted.size() < 4) {
      fill_leading_quadruplet(kinematic);
    } else {
      result.cut_mask |= kElectronCrackAcceptance;
      std::vector<int> isolated;
      std::copy_if(crack_accepted.begin(), crack_accepted.end(),
                   std::back_inserter(isolated), [&](int index) {
                     return passes_isolation(leptons[index]);
                   });
      if (isolated.size() < 4) {
        fill_leading_quadruplet(crack_accepted);
      } else {
        result.cut_mask |= kIsolation;
        const PairingChoice pairing = best_pairing(leptons, isolated);
        if (!pairing.has_pairing) {
          fill_leading_quadruplet(isolated);
        } else {
          result.cut_mask |= kOssfPairing;
          result.output_indices = pairing.lepton_indices;
          result.channel = classify_channel(result.output_indices, leptons);
          result.four_lepton = pairing.four_lepton;
          result.z1 = pairing.z1;
          result.z2 = pairing.z2;
          result.z1_delta_r =
              delta_r(leptons[result.output_indices[0]].reconstructed,
                      leptons[result.output_indices[1]].reconstructed);
          result.z2_delta_r =
              delta_r(leptons[result.output_indices[2]].reconstructed,
                      leptons[result.output_indices[3]].reconstructed);

          if (all_ossf_masses_pass(result.output_indices, leptons)) {
            result.cut_mask |= kAllOssfMasses;
            if (result.z1.mass() > kZ1MassMinGeV &&
                result.z1.mass() < kZ1MassMaxGeV) {
              result.cut_mask |= kZ1Window;
              if (result.z2.mass() > kZ2MassMinGeV &&
                  result.z2.mass() < kZ2MassMaxGeV) {
                result.cut_mask |= kZ2Window;
                result.cut_mask |= kTrigger | kSelected;
              }
            }
          }
        }
      }
    }
  }

  result.min_delta_r = std::numeric_limits<double>::infinity();
  for (int first = 0; first < 4; ++first) {
    for (int second = first + 1; second < 4; ++second) {
      result.min_delta_r =
          std::min(result.min_delta_r,
                   delta_r(leptons[result.output_indices[first]].reconstructed,
                           leptons[result.output_indices[second]].reconstructed));
    }
  }
  return result;
}

double subset_probability(const std::vector<Lepton>& leptons,
                          ULong64_t mask) {
  double probability = 1.0;
  for (std::size_t index = 0; index < leptons.size(); ++index) {
    const bool reconstructed = (mask & (ULong64_t{1} << index)) != 0U;
    probability *= reconstructed ? leptons[index].efficiency
                                 : 1.0 - leptons[index].efficiency;
  }
  return probability;
}

int popcount(ULong64_t value) {
  int count = 0;
  while (value != 0U) {
    value &= value - 1U;
    ++count;
  }
  return count;
}

double trigger_efficiency(int channel, const Config& config) {
  if (config.response_profile == ResponseProfile::Perfect) {
    return 1.0;
  }
  return channel == 1 ? config.four_electron_trigger_efficiency
                      : config.other_trigger_efficiency;
}

RecoContext reconstructed_context(bool jets_available, int num_jets,
                                  const double jets[][kMaxRecoObjects],
                                  bool bjets_available, int num_bjets,
                                  const double bjets[][kMaxRecoObjects],
                                  bool met_available, const double etmiss[4]) {
  RecoContext context;
  if (jets_available || bjets_available) {
    // HwSim stores tagged b jets in thebJets and the remaining jets in
    // theJets.  Expose n_jets and leading_jet_pt as inclusive observables.
    context.n_jets = 0;
    context.leading_jet_pt = 0.0;
  }
  if (jets_available) {
    if (num_jets < 0 || num_jets > kMaxRecoObjects) {
      throw std::runtime_error("numJets is outside the supported range [0,100]");
    }
    for (int index = 0; index < num_jets; ++index) {
      const FourVector jet = from_reco_array(jets, index);
      if (finite_four_vector(jet) && jet.pt() > kJetPtMinGeV &&
          std::fabs(jet.eta()) < kJetEtaMax) {
        ++context.n_jets;
        context.leading_jet_pt = std::max(context.leading_jet_pt, jet.pt());
      }
    }
  }
  if (bjets_available) {
    if (num_bjets < 0 || num_bjets > kMaxRecoObjects) {
      throw std::runtime_error("numbJets is outside the supported range [0,100]");
    }
    context.n_bjets = 0;
    for (int index = 0; index < num_bjets; ++index) {
      const FourVector jet = from_reco_array(bjets, index);
      if (finite_four_vector(jet) && jet.pt() > kJetPtMinGeV &&
          std::fabs(jet.eta()) < kJetEtaMax) {
        ++context.n_jets;
        ++context.n_bjets;
        context.leading_jet_pt = std::max(context.leading_jet_pt, jet.pt());
      }
    }
  }
  if (met_available) {
    const FourVector missing{etmiss[1], etmiss[2], etmiss[3], etmiss[0]};
    context.met_pt = missing.pt();
    context.met_phi = missing.phi();
  }
  return context;
}

void book_source_tree(TTree& tree, SourceRecord& record) {
  tree.Branch("source_index", &record.source_index, "source_index/L");
  tree.Branch("generator_weight", &record.generator_weight, "generator_weight/D");
  tree.Branch("sample_weight", &record.sample_weight, "sample_weight/D");
  tree.Branch("base_event_weight", &record.base_event_weight, "base_event_weight/D");
  tree.Branch("response_failure_probability", &record.response_failure_probability,
              "response_failure_probability/D");
  tree.Branch("response_probability_in_tree", &record.response_probability_in_tree,
              "response_probability_in_tree/D");
  tree.Branch("selected_probability", &record.selected_probability,
              "selected_probability/D");
  tree.Branch("selected_event_weight", &record.selected_event_weight,
              "selected_event_weight/D");
  tree.Branch("response_closure_delta", &record.response_closure_delta,
              "response_closure_delta/D");
  tree.Branch("raw_lepton_count", &record.raw_lepton_count, "raw_lepton_count/I");
  tree.Branch("accepted_lepton_count", &record.accepted_lepton_count,
              "accepted_lepton_count/I");
  tree.Branch("subset_hypothesis_count", &record.subset_hypothesis_count,
              "subset_hypothesis_count/I");
  tree.Branch("cut_probabilities", record.cut_probabilities,
              "cut_probabilities[10]/D");
  tree.Branch("optional_weight_names", &record.optional_weight_names);
  tree.Branch("optional_weights", &record.optional_weights);
  tree.Branch("optional_selected_event_weights",
              &record.optional_selected_event_weights);
}

void book_four_lepton_tree(TTree& tree, FourLeptonRecord& record) {
  tree.Branch("source_index", &record.source_index, "source_index/L");
  tree.Branch("subset_mask", &record.subset_mask, "subset_mask/l");
  tree.Branch("cut_mask", &record.cut_mask, "cut_mask/i");
  tree.Branch("channel", &record.channel, "channel/I");
  tree.Branch("candidate_count", &record.candidate_count, "candidate_count/I");
  tree.Branch("hypothesis_probability", &record.hypothesis_probability,
              "hypothesis_probability/D");
  tree.Branch("generator_weight", &record.generator_weight, "generator_weight/D");
  tree.Branch("sample_weight", &record.sample_weight, "sample_weight/D");
  tree.Branch("response_weight", &record.response_weight, "response_weight/D");
  tree.Branch("trigger_weight", &record.trigger_weight, "trigger_weight/D");
  tree.Branch("pretrigger_event_weight", &record.pretrigger_event_weight,
              "pretrigger_event_weight/D");
  tree.Branch("event_weight", &record.event_weight, "event_weight/D");
  tree.Branch("optional_event_weights", &record.optional_event_weights);
  tree.Branch("m4l", &record.m4l, "m4l/D");
  tree.Branch("mZ1", &record.mZ1, "mZ1/D");
  tree.Branch("mZ2", &record.mZ2, "mZ2/D");
  tree.Branch("pt4l", &record.pt4l, "pt4l/D");
  tree.Branch("y4l", &record.y4l, "y4l/D");
  tree.Branch("min_delta_r", &record.min_delta_r, "min_delta_r/D");
  tree.Branch("z1_delta_r", &record.z1_delta_r, "z1_delta_r/D");
  tree.Branch("z2_delta_r", &record.z2_delta_r, "z2_delta_r/D");
  tree.Branch("n_jets", &record.n_jets, "n_jets/I");
  tree.Branch("n_bjets", &record.n_bjets, "n_bjets/I");
  tree.Branch("leading_jet_pt", &record.leading_jet_pt, "leading_jet_pt/D");
  tree.Branch("met_pt", &record.met_pt, "met_pt/D");
  tree.Branch("met_phi", &record.met_phi, "met_phi/D");
  tree.Branch("lepton_pt", record.lepton_pt, "lepton_pt[4]/D");
  tree.Branch("lepton_eta", record.lepton_eta, "lepton_eta[4]/D");
  tree.Branch("lepton_phi", record.lepton_phi, "lepton_phi[4]/D");
  tree.Branch("lepton_isolation", record.lepton_isolation,
              "lepton_isolation[4]/D");
  tree.Branch("lepton_flavour", record.lepton_flavour, "lepton_flavour[4]/I");
  tree.Branch("lepton_charge", record.lepton_charge, "lepton_charge[4]/I");
  tree.Branch("lepton_source_index", record.lepton_source_index,
              "lepton_source_index[4]/I");
  tree.Branch("lepton_dressed_photon_count",
              record.lepton_dressed_photon_count,
              "lepton_dressed_photon_count[4]/I");
}

void reset_four_lepton_record(FourLeptonRecord& record) {
  record = FourLeptonRecord{};
  record.m4l = -1.0;
  record.mZ1 = -1.0;
  record.mZ2 = -1.0;
  record.pt4l = -1.0;
  record.min_delta_r = -1.0;
  record.z1_delta_r = -1.0;
  record.z2_delta_r = -1.0;
  record.n_jets = -1;
  record.n_bjets = -1;
  record.leading_jet_pt = -1.0;
  record.met_pt = -1.0;
  std::fill(std::begin(record.lepton_source_index),
            std::end(record.lepton_source_index), -1);
}

std::string json_escape(const std::string& input) {
  std::ostringstream output;
  for (const unsigned char character : input) {
    switch (character) {
      case '"':
        output << "\\\"";
        break;
      case '\\':
        output << "\\\\";
        break;
      case '\b':
        output << "\\b";
        break;
      case '\f':
        output << "\\f";
        break;
      case '\n':
        output << "\\n";
        break;
      case '\r':
        output << "\\r";
        break;
      case '\t':
        output << "\\t";
        break;
      default:
        if (character < 0x20) {
          output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<int>(character) << std::dec;
        } else {
          output << static_cast<char>(character);
        }
    }
  }
  return output.str();
}

void write_weight_sums(std::ostream& output, const WeightSums& sums, int indent) {
  const std::string padding(static_cast<std::size_t>(indent), ' ');
  output << "{\n"
         << padding << "  \"entries\": " << sums.entries << ",\n"
         << padding << "  \"positive_entries\": " << sums.positive_entries << ",\n"
         << padding << "  \"negative_entries\": " << sums.negative_entries << ",\n"
         << padding << "  \"zero_entries\": " << sums.zero_entries << ",\n"
         << padding << "  \"negative_fraction\": " << sums.negative_fraction() << ",\n"
         << padding << "  \"sumw\": " << sums.sumw << ",\n"
         << padding << "  \"sumabsw\": " << sums.sumabsw << ",\n"
         << padding << "  \"sumw2\": " << sums.sumw2 << ",\n"
         << padding << "  \"sumw_positive\": " << sums.sumw_positive << ",\n"
         << padding << "  \"sumw_negative\": " << sums.sumw_negative << "\n"
         << padding << "}";
}

void write_summary(const std::filesystem::path& path, const Config& config,
                   const std::filesystem::path& root_path, long long events_read,
                   long long tree_rows, long long selected_rows,
                   const WeightSums& raw_generated,
                   const WeightSums& scaled_generated,
                   const WeightSums& scaled_selected,
                   const Diagnostics& diagnostics) {
  std::ofstream output(path);
  if (!output) {
    throw std::runtime_error("failed to create summary " + path.string());
  }
  const bool ssc_response =
      config.response_profile == ResponseProfile::Ssc;
  const double effective_electron_efficiency =
      ssc_response ? config.electron_efficiency : 1.0;
  const double effective_muon_efficiency =
      ssc_response ? config.muon_efficiency : 1.0;
  const double effective_four_electron_trigger =
      ssc_response ? config.four_electron_trigger_efficiency : 1.0;
  const double effective_other_trigger =
      ssc_response ? config.other_trigger_efficiency : 1.0;
  output << std::setprecision(17);
  output << "{\n"
         << "  \"schema_version\": 1,\n"
         << "  \"analysis\": \"HwSimPostAnalysis_fourlepton\",\n"
         << "  \"input\": \"" << json_escape(config.input_name) << "\",\n"
         << "  \"optional_weight_names_file\": \""
         << json_escape(config.optional_weight_names_file) << "\",\n"
         << "  \"output_root\": \"" << json_escape(root_path.string()) << "\",\n"
         << "  \"response_profile\": \"" << config.response_profile_name << "\",\n"
         << "  \"sample\": \"" << json_escape(config.sample) << "\",\n"
         << "  \"category\": \"" << json_escape(config.category) << "\",\n"
         << "  \"requested_channel\": \"" << json_escape(config.requested_channel)
         << "\",\n"
         << "  \"events_read\": " << events_read << ",\n"
         << "  \"four_lepton_rows\": " << tree_rows << ",\n"
         << "  \"selected_rows\": " << selected_rows << ",\n"
         << "  \"seed\": " << config.seed << ",\n"
         << "  \"weight_scale\": " << config.weight_scale << ",\n"
         << "  \"muon_resolution_scale\": " << config.muon_resolution_scale << ",\n"
         << "  \"weights_are_cross_sections\": false,\n"
         << "  \"event_weight_formula\": "
            "\"generator_weight * sample_weight * response_weight\",\n"
         << "  \"cut_mask_definition\": \"" << cut_mask_definition() << "\",\n"
         << "  \"cut_bit_zero_semantics\": \""
         << json_escape(cut_bit_zero_semantics()) << "\",\n"
         << "  \"channel_definition\": \"1:4e;2:4mu;3:2e2mu;0:unclassified\",\n"
         << "  \"generated_raw\": ";
  write_weight_sums(output, raw_generated, 2);
  output << ",\n  \"generated_scaled\": ";
  write_weight_sums(output, scaled_generated, 2);
  output << ",\n  \"selected_scaled\": ";
  write_weight_sums(output, scaled_selected, 2);
  output << ",\n"
         << "  \"response_closure\": {\n"
         << "    \"failure_count\": " << diagnostics.closure_failures << ",\n"
         << "    \"maximum_absolute_delta\": " << diagnostics.maximum_closure_delta
         << "\n"
         << "  },\n"
         << "  \"diagnostics\": {\n"
         << "    \"invalid_particle_records\": "
         << diagnostics.invalid_particle_records << ",\n"
         << "    \"accepted_lepton_overflow_events\": "
         << diagnostics.lepton_overflow_events << ",\n"
         << "    \"marginalized_nonbaseline_leptons\": "
         << diagnostics.marginalized_nonbaseline_leptons << ",\n"
         << "    \"muon_pt_below_table\": " << diagnostics.muon_pt_underflow << ",\n"
         << "    \"muon_pt_above_table\": " << diagnostics.muon_pt_overflow << ",\n"
         << "    \"optional_weight_name_syntheses\": "
         << diagnostics.optional_weight_name_syntheses << ",\n"
         << "    \"optional_weight_names_from_input\": "
         << diagnostics.optional_weight_names_from_input << ",\n"
         << "    \"optional_weight_names_from_file\": "
         << diagnostics.optional_weight_names_from_file << ",\n"
         << "    \"optional_weight_name_file_validations\": "
         << diagnostics.optional_weight_name_file_validations << ",\n"
         << "    \"optional_weight_name_mismatches\": "
         << diagnostics.optional_weight_name_mismatches << "\n"
         << "  },\n"
         << "  \"detector\": {\n"
         << "    \"parameter_source\": \""
         << (ssc_response
                 ? "GEM_EM_response_and_GEM_TDR_muon_digitization"
                 : "perfect_identity_response")
         << "\",\n"
         << "    \"electron_efficiency\": "
         << effective_electron_efficiency << ",\n"
         << "    \"muon_efficiency\": " << effective_muon_efficiency << ",\n"
         << "    \"trigger_4e_efficiency\": "
         << effective_four_electron_trigger << ",\n"
         << "    \"trigger_other_efficiency\": "
         << effective_other_trigger << ",\n"
         << "    \"em_pileup_noise_enabled\": "
         << (ssc_response && config.include_pileup_noise ? "true" : "false")
         << "\n"
         << "  }\n"
         << "}\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  try {
    if (argc == 1) {
      print_usage(std::cerr);
      return 1;
    }
    const Config config = parse_config(argc, argv);
    const std::vector<std::string> configured_optional_weight_names =
        load_optional_weight_names(config.optional_weight_names_file);

    TChain chain("Data");
    add_input_files(chain, config.input_name);
    const long long total_entries = chain.GetEntries();
    if (total_entries <= 0) {
      throw std::runtime_error("no events found in " + config.input_name);
    }

    const std::array<const char*, 3> required_branches =
        {{"numparticles", "objects", "evweight"}};
    for (const char* branch_name : required_branches) {
      if (chain.GetBranch(branch_name) == nullptr) {
        throw std::runtime_error(std::string("missing required HwSim branch ") +
                                 branch_name);
      }
    }

    const bool has_num_jets = chain.GetBranch("numJets") != nullptr;
    const bool has_jets = chain.GetBranch("theJets") != nullptr;
    const bool has_num_bjets = chain.GetBranch("numbJets") != nullptr;
    const bool has_bjets = chain.GetBranch("thebJets") != nullptr;
    if (has_num_jets != has_jets) {
      throw std::runtime_error(
          "HwSim jet input is incomplete: numJets and theJets must appear together");
    }
    if (has_num_bjets != has_bjets) {
      throw std::runtime_error(
          "HwSim b-jet input is incomplete: numbJets and thebJets must appear together");
    }
    const bool met_available = chain.GetBranch("theETmiss") != nullptr;
    const bool optional_weights_available =
        chain.GetBranch("theOptWeights") != nullptr;
    const bool optional_names_available =
        chain.GetBranch("theOptWeightsNames") != nullptr;
    if (optional_names_available && !optional_weights_available) {
      throw std::runtime_error(
          "theOptWeightsNames exists without theOptWeights");
    }

    int num_particles = 0;
    double objects[8][kMaxParticles] = {};
    double input_event_weight = 0.0;
    int num_jets = 0;
    double jets[5][kMaxRecoObjects] = {};
    int num_bjets = 0;
    double bjets[5][kMaxRecoObjects] = {};
    double etmiss[4] = {};
    std::vector<double>* input_optional_weights = nullptr;
    std::vector<std::string>* input_optional_weight_names = nullptr;

    chain.SetBranchAddress("numparticles", &num_particles);
    chain.SetBranchAddress("objects", objects);
    chain.SetBranchAddress("evweight", &input_event_weight);
    if (has_jets) {
      chain.SetBranchAddress("numJets", &num_jets);
      chain.SetBranchAddress("theJets", jets);
    }
    if (has_bjets) {
      chain.SetBranchAddress("numbJets", &num_bjets);
      chain.SetBranchAddress("thebJets", bjets);
    }
    if (met_available) {
      chain.SetBranchAddress("theETmiss", etmiss);
    }
    if (optional_weights_available) {
      chain.SetBranchAddress("theOptWeights", &input_optional_weights);
    }
    if (optional_names_available) {
      chain.SetBranchAddress("theOptWeightsNames",
                             &input_optional_weight_names);
    }

    long long first_event = std::max(0LL, config.first_event);
    long long last_event =
        config.last_event < 0 ? total_entries : std::min(config.last_event, total_entries);
    if (last_event <= first_event) {
      throw std::runtime_error("requested event range is empty");
    }

    const auto [root_path, summary_path] = output_paths(config);
    TFile output_file(root_path.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
      throw std::runtime_error("failed to create " + root_path.string());
    }

    SourceRecord source_record;
    FourLeptonRecord four_lepton_record;
    TTree source_tree("SourceEvents", "Four-lepton source-event response summary");
    TTree four_lepton_tree("FourLepton",
                           "Weighted four-lepton reconstruction hypotheses");
    book_source_tree(source_tree, source_record);
    book_four_lepton_tree(four_lepton_tree, four_lepton_record);

    int metadata_schema_version = 1;
    std::string metadata_response_profile = config.response_profile_name;
    std::string metadata_sample = config.sample;
    std::string metadata_category = config.category;
    std::string metadata_requested_channel = config.requested_channel;
    std::string metadata_cut_mask = cut_mask_definition();
    std::string metadata_cut_bit_zero_semantics = cut_bit_zero_semantics();
    std::string metadata_channel_definition =
        "1:4e;2:4mu;3:2e2mu;0:unclassified";
    std::string metadata_detector_parameter_source =
        config.response_profile == ResponseProfile::Perfect
            ? "perfect_identity_response"
            : "GEM_EM_response_and_GEM_TDR_muon_digitization";
    std::string metadata_optional_weight_names_file =
        config.optional_weight_names_file;
    std::string metadata_optional_weight_names_source;
    ULong64_t metadata_seed = config.seed;
    double metadata_weight_scale = config.weight_scale;
    double metadata_muon_resolution_scale = config.muon_resolution_scale;
    const bool ssc_response =
        config.response_profile == ResponseProfile::Ssc;
    double metadata_electron_efficiency =
        ssc_response ? config.electron_efficiency : 1.0;
    double metadata_muon_efficiency =
        ssc_response ? config.muon_efficiency : 1.0;
    double metadata_four_electron_trigger_efficiency =
        ssc_response ? config.four_electron_trigger_efficiency : 1.0;
    double metadata_other_trigger_efficiency =
        ssc_response ? config.other_trigger_efficiency : 1.0;
    bool metadata_em_pileup_noise_enabled =
        ssc_response && config.include_pileup_noise;
    TTree metadata_tree("AnalysisMetadata", "Four-lepton analysis metadata");
    metadata_tree.Branch("schema_version", &metadata_schema_version,
                         "schema_version/I");
    metadata_tree.Branch("response_profile", &metadata_response_profile);
    metadata_tree.Branch("sample", &metadata_sample);
    metadata_tree.Branch("category", &metadata_category);
    metadata_tree.Branch("requested_channel", &metadata_requested_channel);
    metadata_tree.Branch("cut_mask_definition", &metadata_cut_mask);
    metadata_tree.Branch("cut_bit_zero_semantics",
                         &metadata_cut_bit_zero_semantics);
    metadata_tree.Branch("channel_definition", &metadata_channel_definition);
    metadata_tree.Branch("detector_parameter_source",
                         &metadata_detector_parameter_source);
    metadata_tree.Branch("optional_weight_names_file",
                         &metadata_optional_weight_names_file);
    metadata_tree.Branch("optional_weight_names_source",
                         &metadata_optional_weight_names_source);
    metadata_tree.Branch("seed", &metadata_seed, "seed/l");
    metadata_tree.Branch("weight_scale", &metadata_weight_scale,
                         "weight_scale/D");
    metadata_tree.Branch("muon_resolution_scale",
                         &metadata_muon_resolution_scale,
                         "muon_resolution_scale/D");
    metadata_tree.Branch("electron_efficiency",
                         &metadata_electron_efficiency,
                         "electron_efficiency/D");
    metadata_tree.Branch("muon_efficiency", &metadata_muon_efficiency,
                         "muon_efficiency/D");
    metadata_tree.Branch("trigger_4e_efficiency",
                         &metadata_four_electron_trigger_efficiency,
                         "trigger_4e_efficiency/D");
    metadata_tree.Branch("trigger_other_efficiency",
                         &metadata_other_trigger_efficiency,
                         "trigger_other_efficiency/D");
    metadata_tree.Branch("em_pileup_noise_enabled",
                         &metadata_em_pileup_noise_enabled,
                         "em_pileup_noise_enabled/O");

    Diagnostics diagnostics;
    WeightSums raw_generated_sums;
    WeightSums scaled_generated_sums;
    WeightSums selected_scaled_sums;
    long long selected_rows = 0;

    std::cout << "Input entries: " << total_entries << '\n'
              << "Analyzing [" << first_event << ", " << last_event << ")\n"
              << "Response profile: " << config.response_profile_name << '\n'
              << "Output ROOT: " << root_path << '\n'
              << "Weight invariant: event_weight = generator_weight * sample_weight "
                 "* response_weight\n";

    for (long long event_index = first_event; event_index < last_event; ++event_index) {
      if (chain.GetEntry(event_index) <= 0) {
        throw std::runtime_error("failed to read source event " +
                                 std::to_string(event_index));
      }
      if (num_particles < 0 || num_particles > kMaxParticles) {
        throw std::runtime_error(
            "numparticles is outside the supported range [0,10000] at event " +
            std::to_string(event_index));
      }
      if (!std::isfinite(input_event_weight)) {
        throw std::runtime_error("non-finite evweight at event " +
                                 std::to_string(event_index));
      }

      std::vector<Particle> particles;
      particles.reserve(static_cast<std::size_t>(num_particles));
      int raw_lepton_count = 0;
      for (int index = 0; index < num_particles; ++index) {
        const FourVector vector = from_array(objects, index);
        const double id_value = objects[4][index];
        if (!finite_four_vector(vector) || !std::isfinite(id_value)) {
          ++diagnostics.invalid_particle_records;
          continue;
        }
        const int pdg_id = static_cast<int>(std::lround(id_value));
        if (std::abs(pdg_id) == 11 || std::abs(pdg_id) == 13) {
          ++raw_lepton_count;
        }
        particles.push_back({vector, pdg_id, index});
      }

      const RecoContext reco_context =
          reconstructed_context(has_jets, num_jets, jets, has_bjets, num_bjets,
                                bjets, met_available, etmiss);
      std::vector<Lepton> leptons =
          response_leptons(particles, config, event_index, diagnostics);
      const int accepted_lepton_count = static_cast<int>(
          std::count_if(leptons.begin(), leptons.end(),
                        [](const Lepton& lepton) {
                          return passes_baseline(lepton);
                        }));
      if (accepted_lepton_count > kMaxEnumeratedLeptons) {
        ++diagnostics.lepton_overflow_events;
        throw std::runtime_error(
            "source event " + std::to_string(event_index) + " contains " +
            std::to_string(accepted_lepton_count) +
            " baseline-accepted leptons; exact response enumeration supports "
            "at most " + std::to_string(kMaxEnumeratedLeptons));
      }
      if (leptons.size() > kMaxEnumeratedLeptons &&
          accepted_lepton_count <= kMaxEnumeratedLeptons) {
        // Candidates that deterministically fail the baseline cannot affect
        // the selected probability. Marginalize all but the most useful
        // cutflow representatives so high-multiplicity heavy-flavour events
        // remain analyzable without exceeding the exact 2^10 enumeration.
        std::stable_sort(
            leptons.begin(), leptons.end(),
            [](const Lepton& left, const Lepton& right) {
              const bool left_baseline = passes_baseline(left);
              const bool right_baseline = passes_baseline(right);
              if (left_baseline != right_baseline) {
                return left_baseline;
              }
              if (left.reconstructed.pt() != right.reconstructed.pt()) {
                return left.reconstructed.pt() > right.reconstructed.pt();
              }
              return left.source_index < right.source_index;
            });
        diagnostics.marginalized_nonbaseline_leptons +=
            static_cast<long long>(leptons.size() - kMaxEnumeratedLeptons);
        leptons.resize(kMaxEnumeratedLeptons);
        std::sort(leptons.begin(), leptons.end(),
                  [](const Lepton& left, const Lepton& right) {
                    return left.source_index < right.source_index;
                  });
      }

      source_record = SourceRecord{};
      source_record.source_index = event_index;
      source_record.generator_weight = input_event_weight;
      source_record.sample_weight = config.weight_scale;
      source_record.base_event_weight =
          input_event_weight * config.weight_scale;
      source_record.raw_lepton_count = raw_lepton_count;
      source_record.accepted_lepton_count = accepted_lepton_count;

      if (optional_weights_available && input_optional_weights != nullptr) {
        source_record.optional_weights = *input_optional_weights;
        if (optional_names_available &&
            input_optional_weight_names != nullptr &&
            !input_optional_weight_names->empty()) {
          source_record.optional_weight_names =
              *input_optional_weight_names;
          ++diagnostics.optional_weight_names_from_input;
          if (!configured_optional_weight_names.empty()) {
            ++diagnostics.optional_weight_name_file_validations;
            if (source_record.optional_weight_names !=
                configured_optional_weight_names) {
              ++diagnostics.optional_weight_name_mismatches;
              throw std::runtime_error(
                  "HwSim optional-weight names disagree with " +
                  config.optional_weight_names_file + " at event " +
                  std::to_string(event_index));
            }
          }
        } else if (!source_record.optional_weights.empty() &&
                   !configured_optional_weight_names.empty()) {
          if (configured_optional_weight_names.size() !=
              source_record.optional_weights.size()) {
            ++diagnostics.optional_weight_name_mismatches;
            throw std::runtime_error(
                "optional-weight names file contains " +
                std::to_string(configured_optional_weight_names.size()) +
                " names but HwSim event " + std::to_string(event_index) +
                " contains " +
                std::to_string(source_record.optional_weights.size()) +
                " weights");
          }
          source_record.optional_weight_names =
              configured_optional_weight_names;
          ++diagnostics.optional_weight_names_from_file;
        } else {
          if (!source_record.optional_weights.empty()) {
            ++diagnostics.optional_weight_name_syntheses;
          }
          source_record.optional_weight_names.reserve(
              source_record.optional_weights.size());
          for (std::size_t index = 0;
               index < source_record.optional_weights.size(); ++index) {
            source_record.optional_weight_names.push_back(
                "optional_" + std::to_string(index));
          }
        }
        if (source_record.optional_weight_names.size() !=
            source_record.optional_weights.size()) {
          ++diagnostics.optional_weight_name_mismatches;
          throw std::runtime_error(
              "optional-weight name/value size mismatch at event " +
              std::to_string(event_index));
        }
      }

      raw_generated_sums.add(input_event_weight);
      scaled_generated_sums.add(source_record.base_event_weight);

      bool source_has_selected_hypothesis = false;
      if (accepted_lepton_count > kMaxEnumeratedLeptons) {
        // Do not silently truncate: unsupported high-multiplicity events are
        // rejected above. Keep this defensive branch if the ordering changes.
        ++diagnostics.lepton_overflow_events;
        throw std::logic_error(
            "accepted-lepton overflow reached response enumeration");
      } else {
        source_record.response_failure_probability = 0.0;
        source_record.response_probability_in_tree = 0.0;
        const ULong64_t subset_count =
            ULong64_t{1} << static_cast<unsigned int>(leptons.size());
        for (ULong64_t mask = 0; mask < subset_count; ++mask) {
          const double probability = subset_probability(leptons, mask);
          if (!(probability > 0.0)) {
            continue;
          }
          const int candidate_count = popcount(mask);
          if (candidate_count < 4) {
            source_record.response_failure_probability += probability;
            continue;
          }
          source_record.response_probability_in_tree += probability;
          ++source_record.subset_hypothesis_count;

          std::vector<int> subset;
          subset.reserve(static_cast<std::size_t>(candidate_count));
          for (std::size_t index = 0; index < leptons.size(); ++index) {
            if ((mask & (ULong64_t{1} << index)) != 0U) {
              subset.push_back(static_cast<int>(index));
            }
          }
          const SubsetResult result = analyze_subset(leptons, subset);
          const bool selected = (result.cut_mask & kSelectedMask) == kSelectedMask;
          const double trigger =
              selected ? trigger_efficiency(result.channel, config) : 1.0;
          const double response_weight = probability * trigger;

          for (int stage = 0; stage < 8; ++stage) {
            const std::uint32_t prefix =
                (1U << static_cast<unsigned int>(stage + 1)) - 1U;
            if ((result.cut_mask & prefix) == prefix) {
              source_record.cut_probabilities[stage] += probability;
            }
          }
          if (selected) {
            source_has_selected_hypothesis = true;
            source_record.cut_probabilities[8] += response_weight;
            source_record.cut_probabilities[9] += response_weight;
            source_record.selected_probability += response_weight;
            ++selected_rows;
          }

          reset_four_lepton_record(four_lepton_record);
          four_lepton_record.source_index = event_index;
          four_lepton_record.subset_mask = mask;
          four_lepton_record.cut_mask = result.cut_mask;
          four_lepton_record.channel = result.channel;
          four_lepton_record.candidate_count = candidate_count;
          four_lepton_record.hypothesis_probability = probability;
          four_lepton_record.generator_weight = input_event_weight;
          four_lepton_record.sample_weight = config.weight_scale;
          four_lepton_record.response_weight = response_weight;
          four_lepton_record.trigger_weight = trigger;
          four_lepton_record.pretrigger_event_weight =
              source_record.base_event_weight * probability;
          four_lepton_record.event_weight =
              input_event_weight * config.weight_scale * response_weight;
          four_lepton_record.m4l = result.four_lepton.mass();
          four_lepton_record.mZ1 =
              (result.cut_mask & kOssfPairing) != 0U ? result.z1.mass() : -1.0;
          four_lepton_record.mZ2 =
              (result.cut_mask & kOssfPairing) != 0U ? result.z2.mass() : -1.0;
          four_lepton_record.pt4l = result.four_lepton.pt();
          four_lepton_record.y4l = result.four_lepton.rapidity();
          four_lepton_record.min_delta_r = result.min_delta_r;
          four_lepton_record.z1_delta_r = result.z1_delta_r;
          four_lepton_record.z2_delta_r = result.z2_delta_r;
          four_lepton_record.n_jets = reco_context.n_jets;
          four_lepton_record.n_bjets = reco_context.n_bjets;
          four_lepton_record.leading_jet_pt = reco_context.leading_jet_pt;
          four_lepton_record.met_pt = reco_context.met_pt;
          four_lepton_record.met_phi = reco_context.met_phi;

          for (int index = 0; index < kLeptonsPerQuadruplet; ++index) {
            const Lepton& lepton = leptons[result.output_indices[index]];
            four_lepton_record.lepton_pt[index] =
                lepton.reconstructed.pt();
            four_lepton_record.lepton_eta[index] =
                lepton.reconstructed.eta();
            four_lepton_record.lepton_phi[index] =
                lepton.reconstructed.phi();
            four_lepton_record.lepton_isolation[index] =
                lepton.isolation_et;
            four_lepton_record.lepton_flavour[index] = lepton.flavour;
            four_lepton_record.lepton_charge[index] = lepton.charge;
            four_lepton_record.lepton_source_index[index] =
                lepton.source_index;
            four_lepton_record.lepton_dressed_photon_count[index] =
                static_cast<int>(lepton.dressed_photon_indices.size());
          }
          four_lepton_record.optional_event_weights.reserve(
              source_record.optional_weights.size());
          for (const double optional_weight : source_record.optional_weights) {
            four_lepton_record.optional_event_weights.push_back(
                optional_weight * config.weight_scale * response_weight);
          }
          four_lepton_tree.Fill();
        }
      }

      source_record.response_closure_delta =
          source_record.response_failure_probability +
          source_record.response_probability_in_tree - 1.0;
      diagnostics.maximum_closure_delta =
          std::max(diagnostics.maximum_closure_delta,
                   std::fabs(source_record.response_closure_delta));
      if (std::fabs(source_record.response_closure_delta) > 1.0e-10) {
        ++diagnostics.closure_failures;
      }
      source_record.selected_event_weight =
          source_record.base_event_weight * source_record.selected_probability;
      if (source_has_selected_hypothesis) {
        selected_scaled_sums.add(source_record.selected_event_weight);
      }
      source_record.optional_selected_event_weights.reserve(
          source_record.optional_weights.size());
      for (const double optional_weight : source_record.optional_weights) {
        source_record.optional_selected_event_weights.push_back(
            optional_weight * config.weight_scale *
            source_record.selected_probability);
      }
      source_tree.Fill();

      if ((event_index - first_event) % 1000 == 0) {
        std::cout << "Processed source event " << event_index << "\r"
                  << std::flush;
      }
    }

    metadata_optional_weight_names_source =
        optional_weight_names_source(
            diagnostics, !configured_optional_weight_names.empty(),
            optional_weights_available);
    metadata_tree.Fill();

    output_file.cd();
    metadata_tree.Write();
    source_tree.Write();
    four_lepton_tree.Write();
    output_file.Close();

    write_summary(summary_path, config, root_path, last_event - first_event,
                  four_lepton_tree.GetEntries(), selected_rows,
                  raw_generated_sums, scaled_generated_sums,
                  selected_scaled_sums, diagnostics);

    std::cout << "\nWrote " << root_path << '\n'
              << "Wrote " << summary_path << '\n'
              << "Source events: " << source_tree.GetEntries() << '\n'
              << "FourLepton rows: " << four_lepton_tree.GetEntries() << '\n'
              << "Selected rows: " << selected_rows << '\n'
              << "Maximum reconstruction closure delta: "
              << diagnostics.maximum_closure_delta << std::endl;
    return diagnostics.closure_failures == 0 ? 0 : 2;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return 1;
  }
}
