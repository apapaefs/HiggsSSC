#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <TBranch.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>

#include "TopHist.h"

namespace {

constexpr int kMaxRecoObjects = 100;
constexpr int kVariableCount = 10;
constexpr double kPi = 3.14159265358979323846;
constexpr double kZMassGeV = 91.1876;

// GEM electromagnetic calorimeter response.  E and pT are in GeV.
constexpr double kEmBarrelSampling = 0.060;
constexpr double kEmEndcapSampling = 0.085;
constexpr double kEmConstant = 0.004;
constexpr double kEmBarrelThermalNoiseEt = 0.100;
constexpr double kEmEndcapThermalNoiseEt = 0.175;
constexpr double kEmCentralPileupNoiseEt = 0.120;

// GEM central and forward hadronic calorimeter response.
constexpr double kJetCentralSampling = 0.60;
constexpr double kJetCentralConstant = 0.04;
constexpr double kJetForwardSampling = 2.00;
constexpr double kJetForwardConstant = 0.06;

// GEM detector and H -> gamma gamma fiducial regions.
constexpr double kEmEtaMin = 0.10;
constexpr double kEmEtaMax = 3.00;
constexpr double kPhotonAnalysisEtaMax = 2.50;
constexpr double kCrackEtaMin = 1.01;
constexpr double kCrackEtaMax = 1.16;
constexpr double kHcalEtaMin = 0.10;
constexpr double kHcalEtaMax = 5.50;
constexpr double kTrackerMuonEtaMax = 2.50;

constexpr double kPhotonPtMinGeV = 20.0;
constexpr double kJetPtMinGeV = 20.0;
constexpr double kElectronPtMinGeV = 10.0;
constexpr double kMuonPtMinGeV = 10.0;
constexpr double kIsolationJetPtMinGeV = 10.0;
constexpr double kIsolationDeltaR = 0.40;
constexpr double kPartonMatchDeltaR = 0.40;

// Conditional efficiencies after explicit pT/eta/crack acceptance.
constexpr double kDefaultPhotonShowerIdEfficiency = 0.90;
constexpr double kDefaultPhotonElectronVetoEfficiency = 0.96;
constexpr double kDefaultPhotonElectronVetoEfficiencyNearZ = 0.86;
constexpr double kDefaultElectronEfficiency = 0.90;
constexpr double kDefaultMuonEfficiency = 0.85 * 0.95;
// The GEM sources do not give a universal offline jet efficiency.
constexpr double kDefaultJetEfficiency = 1.00;

constexpr double kDefaultElectronFakeRateNearZ = 0.0015;
constexpr double kDefaultElectronFakeRateContinuum = 0.020;
constexpr double kDefaultNearZHalfWidthGeV = 10.0;

enum class ResponseMode {
  Genuine,
  GammaJet,
  Dielectron,
};

enum class PhotonOrigin {
  None = 0,
  Genuine = 1,
  Jet = 2,
  Electron = 3,
};

struct FourVector {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double e = 0.0;

  double pt() const { return std::sqrt(px * px + py * py); }

  double pabs() const { return std::sqrt(px * px + py * py + pz * pz); }

  double eta() const {
    const double p = pabs();
    const double plus = p + pz;
    const double minus = p - pz;
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

struct Candidate {
  FourVector original;
  FourVector reconstructed;
  int input_index = -1;
  bool is_b_jet = false;
};

struct JetChoice {
  Candidate candidate;
  int matched_parton_id = 0;
  bool matched = false;
  bool flavor_inferred_from_hard_process = false;
};

struct SelectedHypothesis {
  bool valid = false;
  FourVector first;
  FourVector second;
  PhotonOrigin first_origin = PhotonOrigin::None;
  PhotonOrigin second_origin = PhotonOrigin::None;
  double first_probability = 0.0;
  double second_probability = 0.0;
  double probability = 0.0;
  double source_mass = -1.0;
  int matched_parton_id = 0;
};

struct Config {
  ResponseMode response_mode = ResponseMode::Genuine;
  std::string response_mode_name = "genuine";
  std::uint64_t seed = 14101983ULL;
  double weight_scale = 1.0;
  double photon_shower_id_efficiency = kDefaultPhotonShowerIdEfficiency;
  double photon_electron_veto_efficiency = kDefaultPhotonElectronVetoEfficiency;
  double photon_electron_veto_efficiency_near_z = kDefaultPhotonElectronVetoEfficiencyNearZ;
  double electron_efficiency = kDefaultElectronEfficiency;
  double muon_efficiency = kDefaultMuonEfficiency;
  double jet_efficiency = kDefaultJetEfficiency;
  double electron_fake_rate_near_z = kDefaultElectronFakeRateNearZ;
  double electron_fake_rate_continuum = kDefaultElectronFakeRateContinuum;
  double near_z_half_width_gev = kDefaultNearZHalfWidthGeV;
  double jet_fake_scale = 1.0;
  double electron_fake_scale = 1.0;
  double unmatched_jet_quark_fraction = 1.0;
  bool include_pileup_noise = true;
};

struct ObjectDiagnostics {
  double fiducial_photons = 0.0;
  double reconstructed_photons = 0.0;
  double fiducial_jets = 0.0;
  double reconstructed_jets = 0.0;
  double fiducial_electrons = 0.0;
  double reconstructed_electrons = 0.0;
  double fiducial_muons = 0.0;
  double reconstructed_muons = 0.0;
  int truncated_events = 0;
};

char* get_cmd_option(char** begin, char** end, const std::string& option) {
  char** iterator = std::find(begin, end, option);
  if (iterator != end && ++iterator != end) {
    return *iterator;
  }
  return nullptr;
}

bool cmd_option_exists(char** begin, char** end, const std::string& option) {
  return std::find(begin, end, option) != end;
}

std::string option_string(char** begin, char** end, const std::string& option,
                          const std::string& default_value) {
  if (!cmd_option_exists(begin, end, option)) {
    return default_value;
  }
  char* value = get_cmd_option(begin, end, option);
  if (value == nullptr) {
    throw std::runtime_error("missing value for " + option);
  }
  return value;
}

double option_double(char** begin, char** end, const std::string& option, double default_value) {
  const std::string text = option_string(begin, end, option, "");
  if (text.empty()) {
    return default_value;
  }
  std::size_t parsed = 0;
  const double value = std::stod(text, &parsed);
  if (parsed != text.size() || !std::isfinite(value)) {
    throw std::runtime_error("invalid numeric value for " + option + ": " + text);
  }
  return value;
}

std::uint64_t option_uint64(char** begin, char** end, const std::string& option,
                            std::uint64_t default_value) {
  const std::string text = option_string(begin, end, option, "");
  if (text.empty()) {
    return default_value;
  }
  std::size_t parsed = 0;
  const unsigned long long value = std::stoull(text, &parsed);
  if (parsed != text.size()) {
    throw std::runtime_error("invalid integer value for " + option + ": " + text);
  }
  return static_cast<std::uint64_t>(value);
}

long long option_int64(char** begin, char** end, const std::string& option,
                       long long default_value) {
  const std::string text = option_string(begin, end, option, "");
  if (text.empty()) {
    return default_value;
  }
  std::size_t parsed = 0;
  const long long value = std::stoll(text, &parsed);
  if (parsed != text.size()) {
    throw std::runtime_error("invalid integer value for " + option + ": " + text);
  }
  return value;
}

void require_probability(double value, const std::string& name) {
  if (!(value >= 0.0 && value <= 1.0)) {
    throw std::runtime_error(name + " must lie in [0,1]");
  }
}

double clamp_probability(double value) {
  return std::max(0.0, std::min(1.0, value));
}

ResponseMode parse_response_mode(const std::string& value) {
  if (value == "genuine") {
    return ResponseMode::Genuine;
  }
  if (value == "gammajet" || value == "gamma-jet") {
    return ResponseMode::GammaJet;
  }
  if (value == "dielectron" || value == "ee") {
    return ResponseMode::Dielectron;
  }
  throw std::runtime_error("unknown --response-mode '" + value +
                           "' (expected genuine, gammajet, or dielectron)");
}

Config parse_config(int argc, char* argv[]) {
  Config config;
  config.response_mode_name = option_string(argv, argv + argc, "--response-mode", "genuine");
  config.response_mode = parse_response_mode(config.response_mode_name);
  config.seed = option_uint64(argv, argv + argc, "--seed", config.seed);
  config.weight_scale = option_double(argv, argv + argc, "-w", config.weight_scale);
  config.weight_scale = option_double(argv, argv + argc, "--weight-scale", config.weight_scale);
  config.photon_shower_id_efficiency = option_double(
      argv, argv + argc, "--photon-id-efficiency", config.photon_shower_id_efficiency);
  config.photon_electron_veto_efficiency = option_double(
      argv, argv + argc, "--photon-electron-veto-efficiency",
      config.photon_electron_veto_efficiency);
  config.photon_electron_veto_efficiency_near_z = option_double(
      argv, argv + argc, "--photon-electron-veto-efficiency-near-z",
      config.photon_electron_veto_efficiency_near_z);
  config.electron_efficiency = option_double(
      argv, argv + argc, "--electron-efficiency", config.electron_efficiency);
  config.muon_efficiency = option_double(argv, argv + argc, "--muon-efficiency",
                                         config.muon_efficiency);
  config.jet_efficiency = option_double(argv, argv + argc, "--jet-efficiency",
                                        config.jet_efficiency);
  config.electron_fake_rate_near_z = option_double(
      argv, argv + argc, "--electron-fake-rate-near-z", config.electron_fake_rate_near_z);
  config.electron_fake_rate_continuum = option_double(
      argv, argv + argc, "--electron-fake-rate-continuum",
      config.electron_fake_rate_continuum);
  config.near_z_half_width_gev = option_double(
      argv, argv + argc, "--near-z-half-width", config.near_z_half_width_gev);
  config.jet_fake_scale = option_double(argv, argv + argc, "--jet-fake-scale",
                                        config.jet_fake_scale);
  config.electron_fake_scale = option_double(argv, argv + argc, "--electron-fake-scale",
                                             config.electron_fake_scale);
  config.unmatched_jet_quark_fraction = option_double(
      argv, argv + argc, "--unmatched-jet-quark-fraction",
      config.unmatched_jet_quark_fraction);
  if (cmd_option_exists(argv, argv + argc, "--no-pileup-noise")) {
    config.include_pileup_noise = false;
  }

  require_probability(config.photon_shower_id_efficiency, "photon ID efficiency");
  require_probability(config.photon_electron_veto_efficiency, "photon electron-veto efficiency");
  require_probability(config.photon_electron_veto_efficiency_near_z,
                      "near-Z photon electron-veto efficiency");
  require_probability(config.electron_efficiency, "electron efficiency");
  require_probability(config.muon_efficiency, "muon efficiency");
  require_probability(config.jet_efficiency, "jet efficiency");
  require_probability(config.electron_fake_rate_near_z, "near-Z electron fake rate");
  require_probability(config.electron_fake_rate_continuum, "continuum electron fake rate");
  require_probability(config.unmatched_jet_quark_fraction, "unmatched jet quark fraction");
  if (config.near_z_half_width_gev < 0.0) {
    throw std::runtime_error("near-Z half width must be non-negative");
  }
  if (config.jet_fake_scale < 0.0 || config.electron_fake_scale < 0.0) {
    throw std::runtime_error("fake-rate scale factors must be non-negative");
  }
  return config;
}

std::string replace_extension(std::string path, const std::string& replacement) {
  const std::array<std::string, 2> extensions = {{".root", ".input"}};
  for (const std::string& extension : extensions) {
    const std::size_t position = path.rfind(extension);
    if (position != std::string::npos && position + extension.size() == path.size()) {
      path.replace(position, extension.size(), replacement);
      return path;
    }
  }
  return path + replacement;
}

double delta_phi(const FourVector& first, const FourVector& second) {
  double value = second.phi() - first.phi();
  while (value > kPi) {
    value -= 2.0 * kPi;
  }
  while (value <= -kPi) {
    value += 2.0 * kPi;
  }
  return value;
}

double delta_r(const FourVector& first, const FourVector& second) {
  const double deta = first.eta() - second.eta();
  const double dphi = delta_phi(first, second);
  return std::sqrt(deta * deta + dphi * dphi);
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

FourVector massless_with_energy_direction(double energy, double eta, double phi) {
  const double safe_energy = std::max(1.0e-9, energy);
  const double transverse_momentum = safe_energy / std::cosh(eta);
  return {transverse_momentum * std::cos(phi), transverse_momentum * std::sin(phi),
          transverse_momentum * std::sinh(eta), safe_energy};
}

double gaussian_energy(double energy, double sigma, unsigned int seed) {
  if (!(energy > 0.0) || !std::isfinite(energy) || !(sigma >= 0.0) || !std::isfinite(sigma)) {
    return 1.0e-9;
  }
  TRandom3 random(seed);
  return std::max(1.0e-9, random.Gaus(energy, sigma));
}

double em_noise_et(double abs_eta, bool include_pileup_noise) {
  const double thermal = abs_eta < kCrackEtaMin ? kEmBarrelThermalNoiseEt
                                                : kEmEndcapThermalNoiseEt;
  if (!include_pileup_noise) {
    return thermal;
  }
  const double pileup = abs_eta < 1.4
                            ? kEmCentralPileupNoiseEt
                            : kEmCentralPileupNoiseEt + 0.366 * (abs_eta - 1.4);
  return std::sqrt(thermal * thermal + pileup * pileup);
}

FourVector smear_em(const FourVector& input, bool include_pileup_noise, unsigned int seed) {
  const double abs_eta = std::fabs(input.eta());
  const double sampling = abs_eta < kCrackEtaMin ? kEmBarrelSampling : kEmEndcapSampling;
  const double energy = input.e;
  const double transverse_momentum = input.pt();
  double relative_variance = sampling * sampling / std::max(energy, 1.0e-9) +
                             kEmConstant * kEmConstant;
  if (transverse_momentum > 0.0) {
    const double noise_fraction = em_noise_et(abs_eta, include_pileup_noise) /
                                  transverse_momentum;
    relative_variance += noise_fraction * noise_fraction;
  }
  const double smeared_energy = gaussian_energy(
      energy, energy * std::sqrt(std::max(0.0, relative_variance)), seed);
  return massless_with_energy_direction(smeared_energy, input.eta(), input.phi());
}

FourVector smear_jet(const FourVector& input, unsigned int seed) {
  const double abs_eta = std::fabs(input.eta());
  const bool central = abs_eta <= kEmEtaMax;
  const double sampling = central ? kJetCentralSampling : kJetForwardSampling;
  const double constant = central ? kJetCentralConstant : kJetForwardConstant;
  const double energy = input.e;
  const double sigma = energy * std::sqrt(sampling * sampling / std::max(energy, 1.0e-9) +
                                          constant * constant);
  const double smeared_energy = gaussian_energy(energy, sigma, seed);
  const double scale = energy > 0.0 ? smeared_energy / energy : 0.0;
  return {input.px * scale, input.py * scale, input.pz * scale, smeared_energy};
}

bool in_em_fiducial(const FourVector& object, double eta_max) {
  const double abs_eta = std::fabs(object.eta());
  return object.pt() >= kPhotonPtMinGeV && abs_eta > kEmEtaMin && abs_eta < eta_max &&
         !(abs_eta > kCrackEtaMin && abs_eta < kCrackEtaMax);
}

bool in_electron_fiducial(const FourVector& object) {
  const double abs_eta = std::fabs(object.eta());
  return object.pt() >= kElectronPtMinGeV && abs_eta > kEmEtaMin &&
         abs_eta < kTrackerMuonEtaMax &&
         !(abs_eta > kCrackEtaMin && abs_eta < kCrackEtaMax);
}

bool in_jet_fiducial(const FourVector& object) {
  const double abs_eta = std::fabs(object.eta());
  return object.pt() >= kJetPtMinGeV && abs_eta > kHcalEtaMin && abs_eta < kHcalEtaMax;
}

bool in_isolation_jet_fiducial(const FourVector& object) {
  const double abs_eta = std::fabs(object.eta());
  return object.pt() >= kIsolationJetPtMinGeV && abs_eta > kHcalEtaMin &&
         abs_eta < kHcalEtaMax;
}

bool in_muon_fiducial(const FourVector& object) {
  const double abs_eta = std::fabs(object.eta());
  return object.pt() >= kMuonPtMinGeV && abs_eta > kEmEtaMin &&
         abs_eta < kTrackerMuonEtaMax;
}

bool photon_isolated_from_jets(const FourVector& photon, const std::vector<Candidate>& jets) {
  for (const Candidate& jet : jets) {
    if (jet.reconstructed.pt() >= kIsolationJetPtMinGeV &&
        delta_r(photon, jet.reconstructed) < kIsolationDeltaR) {
      return false;
    }
  }
  return true;
}

double interpolate_clamped(double x, const std::array<double, 5>& xs,
                           const std::array<double, 5>& ys) {
  if (x <= xs.front()) {
    return ys.front();
  }
  if (x >= xs.back()) {
    return ys.back();
  }
  for (std::size_t index = 0; index + 1 < xs.size(); ++index) {
    if (x <= xs[index + 1]) {
      const double fraction = (x - xs[index]) / (xs[index + 1] - xs[index]);
      return ys[index] + fraction * (ys[index + 1] - ys[index]);
    }
  }
  return ys.back();
}

double jet_to_photon_fake_rate(double diphoton_mass, int parton_id,
                               double unmatched_quark_fraction, double scale) {
  const std::array<double, 5> masses = {{80.0, 100.0, 120.0, 140.0, 160.0}};
  const std::array<double, 5> quark_rates = {{4.1e-4, 3.9e-4, 3.7e-4, 3.6e-4, 3.4e-4}};
  const std::array<double, 5> gluon_rates = {{1.2e-4, 1.1e-4, 1.0e-4, 0.91e-4, 0.83e-4}};
  const double quark_rate = interpolate_clamped(diphoton_mass, masses, quark_rates);
  const double gluon_rate = interpolate_clamped(diphoton_mass, masses, gluon_rates);
  double rate = 0.0;
  if (std::abs(parton_id) >= 1 && std::abs(parton_id) <= 6) {
    rate = quark_rate;
  } else if (parton_id == 21) {
    rate = gluon_rate;
  } else {
    rate = unmatched_quark_fraction * quark_rate +
           (1.0 - unmatched_quark_fraction) * gluon_rate;
  }
  return clamp_probability(rate * scale);
}

bool is_near_z(double mass, const Config& config) {
  return std::fabs(mass - kZMassGeV) < config.near_z_half_width_gev;
}

double photon_efficiency_for_mass(double mass, const Config& config) {
  const double electron_veto = is_near_z(mass, config)
                                   ? config.photon_electron_veto_efficiency_near_z
                                   : config.photon_electron_veto_efficiency;
  return clamp_probability(config.photon_shower_id_efficiency * electron_veto);
}

double electron_to_photon_fake_rate(double mass, const Config& config) {
  // GEM's tracker gamma/e rejection is conditional on an isolated EM cluster.
  // Apply the separate common EM shower-ID efficiency here; the explicit
  // isolation requirement is imposed when the electron candidates are built.
  // Electron reconstruction efficiency and genuine-photon track-veto
  // acceptance do not enter this truth-electron -> photon transfer factor.
  const double rate = is_near_z(mass, config) ? config.electron_fake_rate_near_z
                                               : config.electron_fake_rate_continuum;
  return clamp_probability(config.photon_shower_id_efficiency * rate *
                           config.electron_fake_scale);
}

int bounded_count(int count, ObjectDiagnostics& diagnostics) {
  if (count < 0) {
    ++diagnostics.truncated_events;
    return 0;
  }
  if (count > kMaxRecoObjects) {
    ++diagnostics.truncated_events;
    return kMaxRecoObjects;
  }
  return count;
}

FourVector array_four_vector(const double objects[][kMaxRecoObjects], int index) {
  return {objects[1][index], objects[2][index], objects[3][index], objects[0][index]};
}

void sort_by_pt(std::vector<Candidate>& candidates) {
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate& left, const Candidate& right) {
              return left.reconstructed.pt() > right.reconstructed.pt();
            });
}

int match_parton_id(const Candidate& jet, int num_partons,
                    const double partons[5][kMaxRecoObjects]) {
  double best_delta_r = kPartonMatchDeltaR;
  int best_id = 0;
  for (int index = 0; index < std::min(num_partons, kMaxRecoObjects); ++index) {
    const int id = static_cast<int>(std::lround(partons[4][index]));
    if (!(id == 21 || (std::abs(id) >= 1 && std::abs(id) <= 6))) {
      continue;
    }
    const FourVector parton{partons[1][index], partons[2][index], partons[3][index],
                            partons[0][index]};
    const double distance = delta_r(jet.original, parton);
    if (distance < best_delta_r) {
      best_delta_r = distance;
      best_id = id;
    }
  }
  return best_id;
}

int unique_hard_colored_parton_id(int num_partons,
                                  const double partons[5][kMaxRecoObjects]) {
  int colored_id = 0;
  for (int index = 0; index < std::min(num_partons, kMaxRecoObjects); ++index) {
    const int id = static_cast<int>(std::lround(partons[4][index]));
    if (!(id == 21 || (std::abs(id) >= 1 && std::abs(id) <= 6))) {
      continue;
    }
    if (colored_id != 0) {
      return 0;
    }
    colored_id = id;
  }
  return colored_id;
}

bool matches_hard_photon(const Candidate& photon, int num_partons,
                         const double partons[5][kMaxRecoObjects]) {
  for (int index = 0; index < std::min(num_partons, kMaxRecoObjects); ++index) {
    const int id = static_cast<int>(std::lround(partons[4][index]));
    if (id != 22) {
      continue;
    }
    const FourVector parton{partons[1][index], partons[2][index], partons[3][index],
                            partons[0][index]};
    if (delta_r(photon.original, parton) < kPartonMatchDeltaR) {
      return true;
    }
  }
  return false;
}

JetChoice choose_fake_jet(const std::vector<Candidate>& fake_jet_candidates,
                          int num_partons, const double partons[5][kMaxRecoObjects],
                          bool partons_available) {
  JetChoice fallback;
  if (fake_jet_candidates.empty()) {
    return fallback;
  }
  fallback.candidate = fake_jet_candidates.front();
  if (!partons_available) {
    return fallback;
  }
  fallback.matched_parton_id = unique_hard_colored_parton_id(num_partons, partons);
  fallback.flavor_inferred_from_hard_process = fallback.matched_parton_id != 0;
  for (const Candidate& candidate : fake_jet_candidates) {
    const int id = match_parton_id(candidate, num_partons, partons);
    if (id != 0) {
      JetChoice choice;
      choice.candidate = candidate;
      choice.matched_parton_id = id;
      choice.matched = true;
      return choice;
    }
  }
  return fallback;
}

void reset_variables(double variables[kVariableCount]) {
  std::fill(variables, variables + kVariableCount, -999.0);
}

void order_selected_pair(SelectedHypothesis& hypothesis) {
  if (hypothesis.second.pt() > hypothesis.first.pt()) {
    std::swap(hypothesis.first, hypothesis.second);
    std::swap(hypothesis.first_origin, hypothesis.second_origin);
    std::swap(hypothesis.first_probability, hypothesis.second_probability);
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Use: ./HwSimPostAnalysis_gammagamma_SSC [input.root|input.input] [options]\n"
              << "Options: -t tag -n maxevents -nmin first -nmax last -w weight_scale\n"
              << "         --response-mode genuine|gammajet|dielectron --seed integer\n"
              << "         --photon-id-efficiency x --electron-efficiency x\n"
              << "         --muon-efficiency x --jet-efficiency x\n"
              << "         --jet-fake-scale x --electron-fake-scale x --no-pileup-noise"
              << std::endl;
    return 1;
  }

  try {
    const std::string input_name = argv[1];
    const Config config = parse_config(argc, argv);

    std::string tag;
    if (cmd_option_exists(argv, argv + argc, "-t")) {
      tag = "-" + option_string(argv, argv + argc, "-t", "");
    }

    TChain chain("Data");
    if (input_name.size() >= 6 && input_name.substr(input_name.size() - 6) == ".input") {
      std::ifstream input_list(input_name.c_str());
      if (!input_list) {
        throw std::runtime_error("failed to open input list " + input_name);
      }
      std::string path;
      while (input_list >> path) {
        if (path.empty() || path[0] == '#') {
          continue;
        }
        if (chain.Add(path.c_str()) == 0) {
          throw std::runtime_error("failed to add ROOT input " + path);
        }
        std::cout << "Adding " << path << std::endl;
      }
    } else if (input_name.size() >= 5 &&
               input_name.substr(input_name.size() - 5) == ".root") {
      if (chain.Add(input_name.c_str()) == 0) {
        throw std::runtime_error("failed to add ROOT input " + input_name);
      }
      std::cout << "Adding " << input_name << std::endl;
    } else {
      throw std::runtime_error("input must end in .root or .input");
    }

    const std::array<const char*, 15> required_branches = {{
        "thePhotons", "numPhotons", "theJets", "numJets", "thebJets", "numbJets",
        "theElectrons", "numElectrons", "thePositrons", "numPositrons", "theMuons",
        "numMuons", "theantiMuons", "numantiMuons", "evweight"}};
    for (const char* branch_name : required_branches) {
      if (chain.GetBranch(branch_name) == nullptr) {
        throw std::runtime_error(std::string("missing required HwSim branch ") + branch_name);
      }
    }

    int num_photons = 0;
    int num_jets = 0;
    int num_b_jets = 0;
    int num_electrons = 0;
    int num_positrons = 0;
    int num_muons = 0;
    int num_antimuons = 0;
    int num_outgoing = 0;
    double photons[4][kMaxRecoObjects] = {};
    double jets[5][kMaxRecoObjects] = {};
    double b_jets[5][kMaxRecoObjects] = {};
    double electrons[4][kMaxRecoObjects] = {};
    double positrons[4][kMaxRecoObjects] = {};
    double muons[4][kMaxRecoObjects] = {};
    double antimuons[4][kMaxRecoObjects] = {};
    double partons[5][kMaxRecoObjects] = {};
    double event_weight_input = 0.0;

    chain.SetBranchAddress("thePhotons", &photons);
    chain.SetBranchAddress("numPhotons", &num_photons);
    chain.SetBranchAddress("theJets", &jets);
    chain.SetBranchAddress("numJets", &num_jets);
    chain.SetBranchAddress("thebJets", &b_jets);
    chain.SetBranchAddress("numbJets", &num_b_jets);
    chain.SetBranchAddress("theElectrons", &electrons);
    chain.SetBranchAddress("numElectrons", &num_electrons);
    chain.SetBranchAddress("thePositrons", &positrons);
    chain.SetBranchAddress("numPositrons", &num_positrons);
    chain.SetBranchAddress("theMuons", &muons);
    chain.SetBranchAddress("numMuons", &num_muons);
    chain.SetBranchAddress("theantiMuons", &antimuons);
    chain.SetBranchAddress("numantiMuons", &num_antimuons);
    chain.SetBranchAddress("evweight", &event_weight_input);

    const bool partons_available = chain.GetBranch("numoutgoing") != nullptr &&
                                   chain.GetBranch("partons") != nullptr;
    if (partons_available) {
      chain.SetBranchAddress("numoutgoing", &num_outgoing);
      chain.SetBranchAddress("partons", &partons);
    }

    const long long event_count = chain.GetEntries();
    if (event_count <= 0) {
      throw std::runtime_error("no events found in " + input_name);
    }
    long long first_event = 0;
    long long last_event = event_count;
    if (cmd_option_exists(argv, argv + argc, "-n")) {
      last_event = option_int64(argv, argv + argc, "-n", event_count);
    } else if (cmd_option_exists(argv, argv + argc, "-nmax")) {
      last_event = option_int64(argv, argv + argc, "-nmax", event_count);
    }
    if (cmd_option_exists(argv, argv + argc, "-nmin")) {
      first_event = option_int64(argv, argv + argc, "-nmin", 0LL);
    }
    first_event = std::max(0LL, first_event);
    last_event = std::min(event_count, last_event);
    if (last_event <= first_event) {
      throw std::runtime_error("requested event range is empty");
    }

    std::cout << "Total number of events in " << input_name << " : " << event_count << '\n'
              << "Analyzing [" << first_event << ", " << last_event << ")\n"
              << "SSC response mode: " << config.response_mode_name << '\n'
              << "Weighted exclusive response hypotheses: enabled\n"
              << "Hard-process parton matching: " << (partons_available ? "enabled" : "unavailable")
              << '\n'
              << "Analysis random seed: " << config.seed << '\n'
              << "Event weight scale: " << config.weight_scale << std::endl;

    const std::string output_top = replace_extension(input_name, tag + ".top");
    const std::string output_dat = replace_extension(input_name, tag + ".dat");
    const std::string output_evp = replace_extension(input_name, tag + ".evp");
    const std::string output_root = replace_extension(input_name, tag + "_var.root");

    std::ofstream dat_output(output_dat.c_str());
    std::ofstream event_output(output_evp.c_str());
    if (!dat_output || !event_output) {
      throw std::runtime_error("failed to create .dat or .evp output");
    }

    TFile root_output(output_root.c_str(), "RECREATE");
    if (root_output.IsZombie()) {
      throw std::runtime_error("failed to create " + output_root);
    }
    TTree output_tree("Data2", "SSC-smeared gamma gamma data tree");
    double variables[kVariableCount] = {};
    double eventweight[1] = {};
    double generatorweight[1] = {};
    double responseweight[1] = {};
    int photonorigin[2] = {};
    long long sourceevent[1] = {};
    output_tree.Branch("variables", variables, "variables[10]/D");
    output_tree.Branch("eventweight", eventweight, "eventweight[1]/D");
    output_tree.Branch("generatorweight", generatorweight, "generatorweight[1]/D");
    output_tree.Branch("responseweight", responseweight, "responseweight[1]/D");
    output_tree.Branch("photonorigin", photonorigin, "photonorigin[2]/I");
    output_tree.Branch("sourceevent", sourceevent, "sourceevent[1]/L");

    TopHist h_dummy(10, output_top, "dummy histo", 0, 1);
    TopHist h_nphotons(11, output_top, "number of selected photons", -0.5, 10.5);
    TopHist h_pt_photons(80, output_top, "pT of selected photons", 0, 2000);
    TopHist h_eta_photons(60, output_top, "eta of selected photons", -3, 3);
    TopHist h_pt_gamma1(80, output_top, "pT of leading photon", 0, 2000);
    TopHist h_pt_gamma2(80, output_top, "pT of subleading photon", 0, 2000);
    TopHist h_eta_gamma1(60, output_top, "eta of leading photon", -3, 3);
    TopHist h_eta_gamma2(60, output_top, "eta of subleading photon", -3, 3);
    TopHist h_mgg(100, output_top, "diphoton invariant mass", 0, 1000);
    TopHist h_delta_r_gg(80, output_top, "DeltaR of two leading photons", 0, 12);
    TopHist h_delta_phi_gg(64, output_top, "DeltaPhi of two leading photons", -3.2, 3.2);
    TopHist h_pt_gg(80, output_top, "pT of diphoton system", 0, 2000);
    TopHist h_y_gg(60, output_top, "rapidity of diphoton system", -6, 6);

    double sum_weight = 0.0;
    double sum_tree_weight = 0.0;
    double sum_diphoton_weight = 0.0;
    double sum_unweighted_diphoton_probability = 0.0;
    long long selected_hypotheses = 0;
    long long quark_matched_hypotheses = 0;
    long long gluon_matched_hypotheses = 0;
    long long quark_flavor_hypotheses = 0;
    long long gluon_flavor_hypotheses = 0;
    long long hard_parton_flavor_fallback_hypotheses = 0;
    long long unknown_jet_flavor_hypotheses = 0;
    long long unmatched_jet_hypotheses = 0;
    long long matched_hard_photon_hypotheses = 0;
    long long unmatched_hard_photon_hypotheses = 0;
    ObjectDiagnostics diagnostics;

    for (long long event_index = first_event; event_index < last_event; ++event_index) {
      if (chain.GetEntry(event_index) <= 0) {
        throw std::runtime_error("failed to read event " + std::to_string(event_index));
      }
      if ((event_index - first_event) % 1000 == 0) {
        std::cout << "Event number: " << event_index << "\r" << std::flush;
      }
      const double base_weight = event_weight_input * config.weight_scale;
      sum_weight += base_weight;

      std::vector<Candidate> reconstructed_jets;
      std::vector<Candidate> isolation_jets;
      std::vector<Candidate> all_fake_jet_candidates;
      std::vector<Candidate> fake_jet_candidates;
      const int safe_num_jets = bounded_count(num_jets, diagnostics);
      const int safe_num_b_jets = bounded_count(num_b_jets, diagnostics);
      for (int source = 0; source < 2; ++source) {
        const int count = source == 0 ? safe_num_jets : safe_num_b_jets;
        const double (*array)[kMaxRecoObjects] = source == 0 ? jets : b_jets;
        for (int index = 0; index < count; ++index) {
          const FourVector original = array_four_vector(array, index);
          if (!(original.e > 0.0) || !std::isfinite(original.e)) {
            continue;
          }
          Candidate jet;
          jet.original = original;
          jet.reconstructed = smear_jet(
              original, object_seed(config.seed, event_index, 20U + source, index));
          jet.input_index = index;
          jet.is_b_jet = source == 1;
          if (in_isolation_jet_fiducial(jet.reconstructed)) {
            isolation_jets.push_back(jet);
          }
          if (in_jet_fiducial(jet.reconstructed)) {
            reconstructed_jets.push_back(jet);
          }

          Candidate fake = jet;
          fake.reconstructed = smear_em(
              original, config.include_pileup_noise,
              object_seed(config.seed, event_index, 30U + source, index));
          all_fake_jet_candidates.push_back(fake);
          if (in_em_fiducial(fake.reconstructed, kPhotonAnalysisEtaMax)) {
            fake_jet_candidates.push_back(fake);
          }
        }
      }
      sort_by_pt(reconstructed_jets);
      sort_by_pt(all_fake_jet_candidates);
      sort_by_pt(fake_jet_candidates);

      std::vector<Candidate> all_photon_candidates;
      std::vector<Candidate> photon_candidates;
      const int safe_num_photons = bounded_count(num_photons, diagnostics);
      for (int index = 0; index < safe_num_photons; ++index) {
        Candidate candidate;
        candidate.original = array_four_vector(photons, index);
        if (!(candidate.original.e > 0.0) || !std::isfinite(candidate.original.e)) {
          continue;
        }
        candidate.reconstructed = smear_em(
            candidate.original, config.include_pileup_noise,
            object_seed(config.seed, event_index, 1U, index));
        candidate.input_index = index;
        all_photon_candidates.push_back(candidate);
        if (in_em_fiducial(candidate.reconstructed, kPhotonAnalysisEtaMax) &&
            photon_isolated_from_jets(candidate.reconstructed, isolation_jets)) {
          photon_candidates.push_back(candidate);
        }
      }
      sort_by_pt(all_photon_candidates);
      sort_by_pt(photon_candidates);

      std::vector<Candidate> all_electron_candidates;
      std::vector<Candidate> all_positron_candidates;
      std::vector<Candidate> electron_candidates;
      std::vector<Candidate> positron_candidates;
      std::size_t fiducial_electron_count = 0;
      const int safe_num_electrons = bounded_count(num_electrons, diagnostics);
      const int safe_num_positrons = bounded_count(num_positrons, diagnostics);
      for (int charge = 0; charge < 2; ++charge) {
        const int count = charge == 0 ? safe_num_electrons : safe_num_positrons;
        const double (*array)[kMaxRecoObjects] = charge == 0 ? electrons : positrons;
        std::vector<Candidate>& all_destination =
            charge == 0 ? all_electron_candidates : all_positron_candidates;
        std::vector<Candidate>& destination =
            charge == 0 ? electron_candidates : positron_candidates;
        for (int index = 0; index < count; ++index) {
          Candidate candidate;
          candidate.original = array_four_vector(array, index);
          if (!(candidate.original.e > 0.0) || !std::isfinite(candidate.original.e)) {
            continue;
          }
          candidate.reconstructed = smear_em(
              candidate.original, config.include_pileup_noise,
              object_seed(config.seed, event_index, 40U + charge, index));
          candidate.input_index = index;
          all_destination.push_back(candidate);
          if (in_electron_fiducial(candidate.reconstructed)) {
            ++fiducial_electron_count;
          }
          if (in_em_fiducial(candidate.reconstructed, kPhotonAnalysisEtaMax) &&
              photon_isolated_from_jets(candidate.reconstructed, isolation_jets)) {
            destination.push_back(candidate);
          }
        }
      }
      sort_by_pt(all_electron_candidates);
      sort_by_pt(all_positron_candidates);
      sort_by_pt(electron_candidates);
      sort_by_pt(positron_candidates);

      // Weighted expected reconstructed-object counts provide diagnostics for
      // efficiencies that do not enter a particular diphoton response mode.
      diagnostics.fiducial_photons += base_weight * photon_candidates.size();
      diagnostics.reconstructed_photons +=
          base_weight * photon_candidates.size() * config.photon_shower_id_efficiency *
          config.photon_electron_veto_efficiency;
      diagnostics.fiducial_jets += base_weight * reconstructed_jets.size();
      diagnostics.reconstructed_jets +=
          base_weight * reconstructed_jets.size() * config.jet_efficiency;
      diagnostics.fiducial_electrons += base_weight * fiducial_electron_count;
      diagnostics.reconstructed_electrons +=
          base_weight * fiducial_electron_count * config.electron_efficiency;

      std::size_t fiducial_muons = 0;
      const int safe_num_muons = bounded_count(num_muons, diagnostics);
      const int safe_num_antimuons = bounded_count(num_antimuons, diagnostics);
      for (int charge = 0; charge < 2; ++charge) {
        const int count = charge == 0 ? safe_num_muons : safe_num_antimuons;
        const double (*array)[kMaxRecoObjects] = charge == 0 ? muons : antimuons;
        for (int index = 0; index < count; ++index) {
          if (in_muon_fiducial(array_four_vector(array, index))) {
            ++fiducial_muons;
          }
        }
      }
      diagnostics.fiducial_muons += base_weight * fiducial_muons;
      diagnostics.reconstructed_muons +=
          base_weight * fiducial_muons * config.muon_efficiency;

      SelectedHypothesis hypothesis;
      if (config.response_mode == ResponseMode::Genuine && !photon_candidates.empty()) {
        hypothesis.valid = true;
        hypothesis.first = photon_candidates[0].reconstructed;
        hypothesis.first_origin = PhotonOrigin::Genuine;
        double pair_mass = -1.0;
        if (photon_candidates.size() >= 2) {
          hypothesis.second = photon_candidates[1].reconstructed;
          hypothesis.second_origin = PhotonOrigin::Genuine;
          hypothesis.source_mass = (photon_candidates[0].original +
                                    photon_candidates[1].original).mass();
          pair_mass = (hypothesis.first + hypothesis.second).mass();
        } else {
          for (const Candidate& photon : all_photon_candidates) {
            if (photon.input_index != photon_candidates[0].input_index) {
              hypothesis.source_mass =
                  (photon_candidates[0].original + photon.original).mass();
              pair_mass = (hypothesis.first + photon.reconstructed).mass();
              break;
            }
          }
        }
        const double efficiency = photon_efficiency_for_mass(pair_mass, config);
        hypothesis.first_probability = efficiency;
        if (photon_candidates.size() >= 2) {
          hypothesis.second_probability = efficiency;
        }
      } else if (config.response_mode == ResponseMode::GammaJet &&
                 !photon_candidates.empty() && !fake_jet_candidates.empty()) {
        const JetChoice jet_choice = choose_fake_jet(
            fake_jet_candidates, num_outgoing, partons, partons_available);
        const Candidate* hard_photon = &photon_candidates.front();
        bool hard_photon_matched = false;
        if (partons_available) {
          for (const Candidate& photon : photon_candidates) {
            if (matches_hard_photon(photon, num_outgoing, partons)) {
              hard_photon = &photon;
              hard_photon_matched = true;
              break;
            }
          }
        }
        hypothesis.valid = true;
        hypothesis.first = hard_photon->reconstructed;
        hypothesis.second = jet_choice.candidate.reconstructed;
        hypothesis.first_origin = PhotonOrigin::Genuine;
        hypothesis.second_origin = PhotonOrigin::Jet;
        hypothesis.matched_parton_id = jet_choice.matched_parton_id;
        hypothesis.source_mass = (hard_photon->original +
                                  jet_choice.candidate.original).mass();
        const double pair_mass = (hypothesis.first + hypothesis.second).mass();
        const double real_photon_efficiency = photon_efficiency_for_mass(pair_mass, config);
        const double fake_rate = jet_to_photon_fake_rate(
            pair_mass, hypothesis.matched_parton_id, config.unmatched_jet_quark_fraction,
            config.jet_fake_scale);
        hypothesis.first_probability = real_photon_efficiency;
        hypothesis.second_probability = fake_rate;
        if (jet_choice.matched_parton_id == 21) {
          ++gluon_flavor_hypotheses;
        } else if (std::abs(jet_choice.matched_parton_id) >= 1 &&
                   std::abs(jet_choice.matched_parton_id) <= 6) {
          ++quark_flavor_hypotheses;
        }
        if (jet_choice.matched && jet_choice.matched_parton_id == 21) {
          ++gluon_matched_hypotheses;
        } else if (jet_choice.matched &&
                   std::abs(jet_choice.matched_parton_id) >= 1 &&
                   std::abs(jet_choice.matched_parton_id) <= 6) {
          ++quark_matched_hypotheses;
        } else {
          ++unmatched_jet_hypotheses;
          if (jet_choice.flavor_inferred_from_hard_process) {
            ++hard_parton_flavor_fallback_hypotheses;
          } else {
            ++unknown_jet_flavor_hypotheses;
          }
        }
        if (hard_photon_matched) {
          ++matched_hard_photon_hypotheses;
        } else {
          ++unmatched_hard_photon_hypotheses;
        }
      } else if (config.response_mode == ResponseMode::GammaJet &&
                 !photon_candidates.empty()) {
        const Candidate* hard_photon = &photon_candidates.front();
        if (partons_available) {
          for (const Candidate& photon : photon_candidates) {
            if (matches_hard_photon(photon, num_outgoing, partons)) {
              hard_photon = &photon;
              break;
            }
          }
        }
        double pair_mass = -1.0;
        if (!all_fake_jet_candidates.empty()) {
          pair_mass = (hard_photon->reconstructed +
                       all_fake_jet_candidates.front().reconstructed).mass();
          hypothesis.source_mass = (hard_photon->original +
                                    all_fake_jet_candidates.front().original).mass();
        }
        hypothesis.valid = true;
        hypothesis.first = hard_photon->reconstructed;
        hypothesis.first_origin = PhotonOrigin::Genuine;
        hypothesis.first_probability = photon_efficiency_for_mass(pair_mass, config);
      } else if (config.response_mode == ResponseMode::GammaJet &&
                 !fake_jet_candidates.empty()) {
        const JetChoice jet_choice = choose_fake_jet(
            fake_jet_candidates, num_outgoing, partons, partons_available);
        double pair_mass = 125.0;
        if (!all_photon_candidates.empty()) {
          pair_mass = (all_photon_candidates.front().reconstructed +
                       jet_choice.candidate.reconstructed).mass();
          hypothesis.source_mass = (all_photon_candidates.front().original +
                                    jet_choice.candidate.original).mass();
        }
        hypothesis.valid = true;
        hypothesis.first = jet_choice.candidate.reconstructed;
        hypothesis.first_origin = PhotonOrigin::Jet;
        hypothesis.matched_parton_id = jet_choice.matched_parton_id;
        hypothesis.first_probability = jet_to_photon_fake_rate(
            pair_mass, hypothesis.matched_parton_id,
            config.unmatched_jet_quark_fraction, config.jet_fake_scale);
      } else if (config.response_mode == ResponseMode::Dielectron &&
                 !electron_candidates.empty() && !positron_candidates.empty()) {
        hypothesis.valid = true;
        hypothesis.first = electron_candidates.front().reconstructed;
        hypothesis.second = positron_candidates.front().reconstructed;
        hypothesis.first_origin = PhotonOrigin::Electron;
        hypothesis.second_origin = PhotonOrigin::Electron;
        hypothesis.source_mass = (electron_candidates.front().original +
                                  positron_candidates.front().original).mass();
        const double fake_rate = electron_to_photon_fake_rate(hypothesis.source_mass, config);
        hypothesis.first_probability = fake_rate;
        hypothesis.second_probability = fake_rate;
      } else if (config.response_mode == ResponseMode::Dielectron &&
                 (!electron_candidates.empty() || !positron_candidates.empty())) {
        const bool use_electron = !electron_candidates.empty();
        const Candidate& selected = use_electron ? electron_candidates.front()
                                                 : positron_candidates.front();
        const std::vector<Candidate>& opposite =
            use_electron ? all_positron_candidates : all_electron_candidates;
        hypothesis.valid = true;
        hypothesis.first = selected.reconstructed;
        hypothesis.first_origin = PhotonOrigin::Electron;
        if (!opposite.empty()) {
          hypothesis.source_mass = (selected.original + opposite.front().original).mass();
        }
        hypothesis.first_probability =
            electron_to_photon_fake_rate(hypothesis.source_mass, config);
      }

      if (hypothesis.valid) {
        hypothesis.first_probability = clamp_probability(hypothesis.first_probability);
        hypothesis.second_probability = clamp_probability(hypothesis.second_probability);
        hypothesis.probability = hypothesis.first_probability * hypothesis.second_probability;
        order_selected_pair(hypothesis);
      }

      const double selected_probability = hypothesis.valid ? hypothesis.probability : 0.0;
      sum_unweighted_diphoton_probability += selected_probability;
      sourceevent[0] = event_index;
      generatorweight[0] = base_weight;

      // Enumerate the exclusive 0-, 1-, and (when both source objects pass
      // kinematic/isolation acceptance) 2-candidate detector outcomes.  Their
      // probabilities sum to one, so the output-tree weights close to the
      // generator input without dividing away rare fake probabilities.
      const double first_probability = hypothesis.valid ? hypothesis.first_probability : 0.0;
      const double second_probability = hypothesis.valid ? hypothesis.second_probability : 0.0;
      const double no_candidate_probability =
          (1.0 - first_probability) * (1.0 - second_probability);
      const double first_only_probability = first_probability * (1.0 - second_probability);
      const double second_only_probability = (1.0 - first_probability) * second_probability;

      reset_variables(variables);
      variables[9] = 0.0;
      eventweight[0] = base_weight * no_candidate_probability;
      responseweight[0] = no_candidate_probability;
      photonorigin[0] = static_cast<int>(PhotonOrigin::None);
      photonorigin[1] = static_cast<int>(PhotonOrigin::None);
      output_tree.Fill();
      sum_tree_weight += eventweight[0];
      h_nphotons.thfill(0.0, eventweight[0]);

      const auto fill_single_candidate = [&](const FourVector& candidate, PhotonOrigin origin,
                                             double probability) {
        if (!(probability > 0.0)) {
          return;
        }
        const double single_weight = base_weight * probability;
        reset_variables(variables);
        variables[1] = candidate.pt();
        variables[2] = candidate.eta();
        variables[9] = 1.0;
        eventweight[0] = single_weight;
        responseweight[0] = probability;
        photonorigin[0] = static_cast<int>(origin);
        photonorigin[1] = static_cast<int>(PhotonOrigin::None);
        output_tree.Fill();
        sum_tree_weight += single_weight;
        h_nphotons.thfill(1.0, single_weight);
        h_pt_photons.thfill(candidate.pt(), single_weight);
        h_eta_photons.thfill(candidate.eta(), single_weight);
        h_pt_gamma1.thfill(candidate.pt(), single_weight);
        h_eta_gamma1.thfill(candidate.eta(), single_weight);
      };

      if (hypothesis.valid) {
        fill_single_candidate(hypothesis.first, hypothesis.first_origin,
                              first_only_probability);
        fill_single_candidate(hypothesis.second, hypothesis.second_origin,
                              second_only_probability);
      }

      if (hypothesis.valid && selected_probability > 0.0) {
        ++selected_hypotheses;
        const FourVector diphoton = hypothesis.first + hypothesis.second;
        const double diphoton_weight = base_weight * selected_probability;
        const double pair_delta_r = delta_r(hypothesis.first, hypothesis.second);
        const double pair_delta_phi = delta_phi(hypothesis.first, hypothesis.second);

        reset_variables(variables);
        variables[0] = diphoton.mass();
        variables[1] = hypothesis.first.pt();
        variables[2] = hypothesis.first.eta();
        variables[3] = hypothesis.second.pt();
        variables[4] = hypothesis.second.eta();
        variables[5] = pair_delta_r;
        variables[6] = pair_delta_phi;
        variables[7] = diphoton.pt();
        variables[8] = diphoton.rapidity();
        variables[9] = 2.0;
        eventweight[0] = diphoton_weight;
        responseweight[0] = selected_probability;
        photonorigin[0] = static_cast<int>(hypothesis.first_origin);
        photonorigin[1] = static_cast<int>(hypothesis.second_origin);
        output_tree.Fill();
        sum_tree_weight += diphoton_weight;
        sum_diphoton_weight += diphoton_weight;
        event_output << event_index << '\n';

        h_nphotons.thfill(2.0, diphoton_weight);
        h_pt_photons.thfill(hypothesis.first.pt(), diphoton_weight);
        h_pt_photons.thfill(hypothesis.second.pt(), diphoton_weight);
        h_eta_photons.thfill(hypothesis.first.eta(), diphoton_weight);
        h_eta_photons.thfill(hypothesis.second.eta(), diphoton_weight);
        h_pt_gamma1.thfill(hypothesis.first.pt(), diphoton_weight);
        h_pt_gamma2.thfill(hypothesis.second.pt(), diphoton_weight);
        h_eta_gamma1.thfill(hypothesis.first.eta(), diphoton_weight);
        h_eta_gamma2.thfill(hypothesis.second.eta(), diphoton_weight);
        h_mgg.thfill(diphoton.mass(), diphoton_weight);
        h_delta_r_gg.thfill(pair_delta_r, diphoton_weight);
        h_delta_phi_gg.thfill(pair_delta_phi, diphoton_weight);
        h_pt_gg.thfill(diphoton.pt(), diphoton_weight);
        h_y_gg.thfill(diphoton.rapidity(), diphoton_weight);
      }
    }

    const long long output_entries = output_tree.GetEntries();
    root_output.cd();
    output_tree.Write();
    root_output.Close();
    std::cout << "\nA ROOT tree has been written to: " << output_root << std::endl;

    h_dummy.thfill(0.5);
    h_dummy.plot(1, 0);
    h_nphotons.add(output_top, 0, 0);
    h_pt_photons.add(output_top, 0, 0);
    h_eta_photons.add(output_top, 0, 0);
    h_pt_gamma1.add(output_top, 0, 0);
    h_eta_gamma1.add(output_top, 0, 0);
    h_pt_gamma2.add(output_top, 0, 0);
    h_eta_gamma2.add(output_top, 0, 0);
    h_mgg.add(output_top, 0, 0);
    h_delta_r_gg.add(output_top, 0, 0);
    h_delta_phi_gg.add(output_top, 0, 0);
    h_pt_gg.add(output_top, 0, 0);
    h_y_gg.add(output_top, 0, 0);

    const double closure_difference = sum_tree_weight - sum_weight;
    dat_output << std::setprecision(15);
    dat_output << "# HwSimPostAnalysis_gammagamma_SSC summary\n";
    dat_output << "input " << input_name << '\n';
    dat_output << "analysis SSC_GEM_weighted_response\n";
    dat_output << "detector_response ssc\n";
    dat_output << "response_mode " << config.response_mode_name << '\n';
    dat_output << "events_read " << (last_event - first_event) << '\n';
    dat_output << "valid_diphoton_hypotheses " << selected_hypotheses << '\n';
    // Retain the legacy key consumed by make_gammagamma_report.py.
    dat_output << "events_with_two_selected_photons " << selected_hypotheses << '\n';
    dat_output << "weighted_hypotheses 1\n";
    dat_output << "tree_entries " << output_entries << '\n';
    dat_output << "weight_scale " << config.weight_scale << '\n';
    dat_output << "sum_weight " << sum_weight << '\n';
    dat_output << "sum_tree_weight " << sum_tree_weight << '\n';
    dat_output << "tree_weight_closure_difference " << closure_difference << '\n';
    dat_output << "sum_diphoton_weight " << sum_diphoton_weight << '\n';
    dat_output << "sum_unweighted_diphoton_probability "
               << sum_unweighted_diphoton_probability << '\n';
    dat_output << "seed " << config.seed << '\n';
    dat_output << "hard_partons_available " << (partons_available ? 1 : 0) << '\n';
    dat_output << "quark_matched_hypotheses " << quark_matched_hypotheses << '\n';
    dat_output << "gluon_matched_hypotheses " << gluon_matched_hypotheses << '\n';
    dat_output << "quark_flavor_hypotheses " << quark_flavor_hypotheses << '\n';
    dat_output << "gluon_flavor_hypotheses " << gluon_flavor_hypotheses << '\n';
    dat_output << "hard_parton_flavor_fallback_hypotheses "
               << hard_parton_flavor_fallback_hypotheses << '\n';
    dat_output << "unknown_jet_flavor_hypotheses "
               << unknown_jet_flavor_hypotheses << '\n';
    dat_output << "unmatched_jet_hypotheses " << unmatched_jet_hypotheses << '\n';
    dat_output << "matched_hard_photon_hypotheses "
               << matched_hard_photon_hypotheses << '\n';
    dat_output << "unmatched_hard_photon_hypotheses "
               << unmatched_hard_photon_hypotheses << '\n';
    dat_output << "photon_pt_min_gev " << kPhotonPtMinGeV << '\n';
    dat_output << "photon_abs_eta_min " << kEmEtaMin << '\n';
    dat_output << "photon_abs_eta_max " << kPhotonAnalysisEtaMax << '\n';
    dat_output << "photon_crack_eta_min " << kCrackEtaMin << '\n';
    dat_output << "photon_crack_eta_max " << kCrackEtaMax << '\n';
    dat_output << "photon_isolation_delta_r " << kIsolationDeltaR << '\n';
    dat_output << "photon_isolation_jet_pt_min_gev " << kIsolationJetPtMinGeV << '\n';
    dat_output << "photon_id_efficiency " << config.photon_shower_id_efficiency << '\n';
    dat_output << "photon_electron_veto_efficiency "
               << config.photon_electron_veto_efficiency << '\n';
    dat_output << "photon_electron_veto_efficiency_near_z "
               << config.photon_electron_veto_efficiency_near_z << '\n';
    dat_output << "electron_efficiency " << config.electron_efficiency << '\n';
    dat_output << "muon_efficiency " << config.muon_efficiency << '\n';
    dat_output << "jet_efficiency " << config.jet_efficiency << '\n';
    dat_output << "jet_efficiency_is_assumption 1\n";
    dat_output << "electron_fake_rate_near_z " << config.electron_fake_rate_near_z << '\n';
    dat_output << "electron_fake_rate_continuum "
               << config.electron_fake_rate_continuum << '\n';
    dat_output << "electron_fake_includes_photon_id_efficiency 1\n";
    dat_output << "electron_fake_requires_jet_isolation 1\n";
    dat_output << "near_z_half_width_gev " << config.near_z_half_width_gev << '\n';
    dat_output << "jet_fake_scale " << config.jet_fake_scale << '\n';
    dat_output << "electron_fake_scale " << config.electron_fake_scale << '\n';
    dat_output << "unmatched_jet_quark_fraction "
               << config.unmatched_jet_quark_fraction << '\n';
    dat_output << "em_pileup_noise_enabled " << (config.include_pileup_noise ? 1 : 0) << '\n';
    dat_output << "fiducial_photon_object_weight " << diagnostics.fiducial_photons << '\n';
    dat_output << "reconstructed_photon_object_weight "
               << diagnostics.reconstructed_photons << '\n';
    dat_output << "fiducial_jet_object_weight " << diagnostics.fiducial_jets << '\n';
    dat_output << "reconstructed_jet_object_weight " << diagnostics.reconstructed_jets << '\n';
    dat_output << "fiducial_electron_object_weight " << diagnostics.fiducial_electrons << '\n';
    dat_output << "reconstructed_electron_object_weight "
               << diagnostics.reconstructed_electrons << '\n';
    dat_output << "fiducial_muon_object_weight " << diagnostics.fiducial_muons << '\n';
    dat_output << "reconstructed_muon_object_weight " << diagnostics.reconstructed_muons << '\n';
    dat_output << "truncated_object_count_collections " << diagnostics.truncated_events << '\n';

    std::cout << "------------------\n"
              << "valid diphoton hypotheses    = " << selected_hypotheses << '\n'
              << "sum input weight             = " << sum_weight << '\n'
              << "sum output-tree weight       = " << sum_tree_weight << '\n'
              << "sum diphoton weight          = " << sum_diphoton_weight << '\n'
              << "tree-weight closure delta    = " << closure_difference << '\n'
              << "------------------" << std::endl;
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << std::endl;
    return 1;
  }
}
