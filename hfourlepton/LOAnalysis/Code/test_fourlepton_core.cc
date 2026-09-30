#define main fourlepton_analyzer_program_main
#include "HwSimPostAnalysis_fourlepton.cc"
#undef main

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int checks = 0;

void require(bool condition, const std::string& message) {
  ++checks;
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void require_close(double actual, double expected, double tolerance,
                   const std::string& message) {
  ++checks;
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::ostringstream detail;
    detail << message << ": expected " << std::setprecision(17) << expected
           << ", got " << actual;
    throw std::runtime_error(detail.str());
  }
}

Lepton test_lepton(int pdg_id, double pt, double eta, double phi,
                    int source_index, double efficiency = 1.0) {
  Lepton lepton;
  const int flavour = std::abs(pdg_id);
  const double mass = flavour == 13 ? kMuonMassGeV : 0.0;
  lepton.bare = with_pt_eta_phi_mass(pt, eta, phi, mass);
  lepton.dressed = lepton.bare;
  lepton.reconstructed = lepton.bare;
  lepton.pdg_id = pdg_id;
  lepton.flavour = flavour;
  lepton.charge = pdg_id > 0 ? -1 : 1;
  lepton.source_index = source_index;
  lepton.isolation_et = 0.0;
  lepton.efficiency = efficiency;
  return lepton;
}

Particle test_particle(int pdg_id, double pt, double eta, double phi,
                       int source_index, double mass = 0.0) {
  return {with_pt_eta_phi_mass(pt, eta, phi, mass), pdg_id, source_index};
}

double back_to_back_pt_for_mass(double pair_mass, double particle_mass = 0.0) {
  return std::sqrt(0.25 * pair_mass * pair_mass -
                   particle_mass * particle_mass);
}

std::vector<Lepton> mixed_flavour_quadruplet(double electron_pair_mass,
                                             double muon_pair_mass) {
  const double muon_pt = std::max(20.0, 0.60 * muon_pair_mass);
  const double muon_sine_half_angle =
      std::sqrt(muon_pair_mass * muon_pair_mass -
                4.0 * kMuonMassGeV * kMuonMassGeV) /
      (2.0 * muon_pt);
  const double muon_opening_angle =
      2.0 * std::asin(muon_sine_half_angle);
  return {
      test_lepton(11, 0.5 * electron_pair_mass, 0.0, 0.0, 0),
      test_lepton(-11, 0.5 * electron_pair_mass, 0.0, kPi, 1),
      test_lepton(13, muon_pt, 0.0, 1.0, 2),
      test_lepton(-13, muon_pt, 0.0, 1.0 + muon_opening_angle, 3),
  };
}

void test_muon_resolution_table() {
  Diagnostics diagnostics;
  require_close(muon_fractional_resolution(10.0, 0.1, diagnostics), 0.028,
                1.0e-14, "10 GeV central muon-resolution anchor");
  require_close(muon_fractional_resolution(500.0, 0.0, diagnostics), 0.05,
                1.0e-14, "500 GeV eta=0 muon-resolution anchor");
  require_close(muon_fractional_resolution(500.0, 2.5, diagnostics), 0.12,
                1.0e-14, "500 GeV eta=2.5 muon-resolution anchor");
  require_close(muon_fractional_resolution(500.0, 1.25, diagnostics), 0.085,
                1.0e-14, "500 GeV eta interpolation");

  const double geometric_midpoint_pt = std::sqrt(10.0 * 25.0);
  require_close(
      muon_fractional_resolution(geometric_midpoint_pt, 0.1, diagnostics),
      0.5 * (0.028 + 0.0167), 1.0e-14,
      "logarithmic pT interpolation");

  Diagnostics clamp_diagnostics;
  const double below =
      muon_fractional_resolution(1.0, 0.1, clamp_diagnostics);
  const double at_low =
      muon_fractional_resolution(10.0, 0.1, clamp_diagnostics);
  const double above =
      muon_fractional_resolution(5000.0, 2.5, clamp_diagnostics);
  const double at_high =
      muon_fractional_resolution(500.0, 2.5, clamp_diagnostics);
  require_close(below, at_low, 1.0e-14, "low-pT resolution clamp");
  require_close(above, at_high, 1.0e-14, "high-pT resolution clamp");
  require(clamp_diagnostics.muon_pt_underflow == 1,
          "low-pT clamp diagnostic count");
  require(clamp_diagnostics.muon_pt_overflow == 1,
          "high-pT clamp diagnostic count");
}

void test_deterministic_smearing() {
  const FourVector electron =
      with_pt_eta_phi_mass(50.0, 0.45, 1.1, 0.0);
  const FourVector electron_a =
      smear_electron(electron, true, 123456U);
  const FourVector electron_b =
      smear_electron(electron, true, 123456U);
  const FourVector electron_c =
      smear_electron(electron, true, 654321U);
  require_close(electron_a.e, electron_b.e, 0.0,
                "electron smearing is deterministic");
  require_close(electron_a.eta(), electron.eta(), 1.0e-12,
                "electron smearing preserves eta");
  require_close(delta_phi(electron_a, electron), 0.0, 1.0e-12,
                "electron smearing preserves phi");
  require(std::fabs(electron_a.e - electron_c.e) > 1.0e-10,
          "electron seed changes the response draw");

  const FourVector muon =
      with_pt_eta_phi_mass(100.0, 0.7, -0.4, kMuonMassGeV);
  Diagnostics diagnostics_a;
  Diagnostics diagnostics_b;
  Diagnostics diagnostics_c;
  const FourVector muon_a =
      smear_muon(muon, -1, 1.0, 10101U, diagnostics_a);
  const FourVector muon_b =
      smear_muon(muon, -1, 1.0, 10101U, diagnostics_b);
  const FourVector muon_c =
      smear_muon(muon, -1, 1.0, 20202U, diagnostics_c);
  require_close(muon_a.pt(), muon_b.pt(), 0.0,
                "muon smearing is deterministic");
  require_close(muon_a.eta(), muon.eta(), 1.0e-12,
                "muon smearing preserves eta");
  require_close(delta_phi(muon_a, muon), 0.0, 1.0e-12,
                "muon smearing preserves phi");
  require_close(muon_a.mass(), kMuonMassGeV, 1.0e-9,
                "muon smearing preserves mass");
  const FourVector composite_input =
      with_pt_eta_phi_mass(100.0, 0.7, -0.4, 0.8);
  Diagnostics composite_diagnostics;
  const FourVector track_only =
      smear_muon(composite_input, -1, 1.0, 10101U,
                 composite_diagnostics);
  require_close(track_only.mass(), kMuonMassGeV, 1.0e-9,
                "curvature response always uses the physical muon mass");
  require(std::fabs(muon_a.pt() - muon_c.pt()) > 1.0e-10,
          "muon seed changes the curvature draw");
}

void test_dressing_and_isolation() {
  std::vector<Particle> particles = {
      test_particle(11, 20.0, 0.0, 0.0, 0),
      test_particle(-11, 20.0, 0.0, 1.0, 1),
      test_particle(22, 2.0, 0.01, 0.04, 2),
      test_particle(211, 1.0, 0.02, 0.20, 3, 0.13957),
      test_particle(12, 8.0, 0.02, 0.22, 4),
      test_particle(211, 3.0, 0.0, 2.0, 5, 0.13957),
  };
  std::vector<Lepton> leptons = build_dressed_leptons(particles);
  require(leptons.size() == 2, "dressing finds both stable leptons");
  require(leptons[0].dressed_photon_indices.size() == 1,
          "nearest electron receives the dressing photon");
  require(leptons[0].dressed_photon_indices[0] == 2,
          "dressing records the photon source index");
  require(leptons[1].dressed_photon_indices.empty(),
          "non-nearest electron does not receive the photon");
  require(leptons[0].dressed.e > leptons[0].bare.e,
          "dressing adds photon four-momentum");

  leptons[0].reconstructed = leptons[0].dressed;
  require_close(isolation_et(leptons[0], particles),
                particles[3].p4.transverse_energy(), 1.0e-12,
                "isolation excludes candidate, dressing photon, neutrino, "
                "and out-of-cone particle");

  std::vector<Particle> muon_particles = {
      test_particle(13, 40.0, 0.2, 0.3, 10, kMuonMassGeV),
      test_particle(22, 3.0, 0.21, 0.31, 11),
  };
  Config config;
  config.response_profile = ResponseProfile::Ssc;
  config.seed = 424242U;
  Diagnostics response_diagnostics;
  const std::vector<Lepton> response =
      response_leptons(muon_particles, config, 7, response_diagnostics);
  Diagnostics expected_diagnostics;
  const FourVector expected =
      smear_muon(
          muon_particles[0].p4, -1, 1.0,
          object_seed(config.seed, 7, 13U, 10U), expected_diagnostics) +
      muon_particles[1].p4;
  require_close(response[0].reconstructed.px, expected.px, 1.0e-12,
                "muon response smears the bare track before photon dressing");
  require_close(response[0].reconstructed.e, expected.e, 1.0e-12,
                "muon response adds the unsmeared dressing photon");
}

void test_open_boundaries() {
  std::vector<Particle> dressing_boundary = {
      test_particle(11, 20.0, 0.0, 0.0, 0),
      test_particle(
          22, 1.0, 0.0,
          std::nextafter(kDressingDeltaR,
                         std::numeric_limits<double>::infinity()),
          1),
  };
  const std::vector<Lepton> boundary_dressed =
      build_dressed_leptons(dressing_boundary);
  require_close(delta_r(dressing_boundary[0].p4, dressing_boundary[1].p4),
                kDressingDeltaR, 1.0e-15,
                "dressing boundary fixture is at DeltaR=0.1");
  require(delta_r(dressing_boundary[0].p4, dressing_boundary[1].p4) >=
              kDressingDeltaR,
          "dressing boundary fixture lies on the excluded side");
  require(boundary_dressed[0].dressed_photon_indices.empty(),
          "photon at DeltaR=0.1 is not dressed");

  dressing_boundary[1] =
      test_particle(22, 1.0, 0.0,
                    std::nextafter(kDressingDeltaR, 0.0), 1);
  const std::vector<Lepton> inside_dressed =
      build_dressed_leptons(dressing_boundary);
  require(inside_dressed[0].dressed_photon_indices.size() == 1,
          "photon immediately inside DeltaR=0.1 is dressed");

  Lepton isolation_candidate = test_lepton(11, 20.0, 0.0, 0.0, 0);
  std::vector<Particle> isolation_particles = {
      test_particle(11, 20.0, 0.0, 0.0, 0),
      test_particle(211, kIsolationConstituentEtMinGeV, 0.0, 0.2, 1),
  };
  require_close(isolation_et(isolation_candidate, isolation_particles), 0.0,
                1.0e-15,
                "constituent with ET=0.5 GeV is excluded from isolation");
  isolation_particles[1] =
      test_particle(211,
                    std::nextafter(kIsolationConstituentEtMinGeV,
                                   std::numeric_limits<double>::infinity()),
                    0.0, 0.2, 1);
  require(isolation_et(isolation_candidate, isolation_particles) >
              kIsolationConstituentEtMinGeV,
          "constituent immediately above ET=0.5 GeV enters isolation");
  isolation_particles[1] =
      test_particle(
          211, 1.0, 0.0,
          std::nextafter(kIsolationDeltaR,
                         std::numeric_limits<double>::infinity()),
          1);
  require_close(delta_r(isolation_candidate.reconstructed,
                        isolation_particles[1].p4),
                kIsolationDeltaR, 1.0e-15,
                "isolation boundary fixture is at DeltaR=0.35");
  require(delta_r(isolation_candidate.reconstructed,
                  isolation_particles[1].p4) >= kIsolationDeltaR,
          "isolation boundary fixture lies on the excluded side");
  require_close(isolation_et(isolation_candidate, isolation_particles), 0.0,
                1.0e-15,
                "constituent at DeltaR=0.35 is outside isolation cone");

  Lepton kinematic = test_lepton(11, kLeptonPtMinGeV, 0.0, 0.0, 0);
  require(!passes_kinematic_acceptance(kinematic),
          "lepton with pT=10 GeV fails the open threshold");
  kinematic =
      test_lepton(11, 20.0, kLeptonEtaMax + 1.0e-12, 0.0, 0);
  require_close(std::fabs(kinematic.reconstructed.eta()), kLeptonEtaMax,
                2.0e-12, "eta boundary fixture is at |eta|=2.5");
  require(std::fabs(kinematic.reconstructed.eta()) >= kLeptonEtaMax,
          "eta boundary fixture lies on the excluded side");
  require(!passes_kinematic_acceptance(kinematic),
          "lepton with |eta|=2.5 fails the open acceptance");

  Lepton crack_low =
      test_lepton(11, 20.0, kElectronCrackEtaMin - 1.0e-12, 0.0, 0);
  Lepton crack_high =
      test_lepton(11, 20.0, kElectronCrackEtaMax + 1.0e-12, 0.0, 0);
  Lepton crack_inside =
      test_lepton(11, 20.0,
                  0.5 * (kElectronCrackEtaMin + kElectronCrackEtaMax),
                  0.0, 0);
  require(passes_electron_crack_acceptance(crack_low),
          "lower electron-crack endpoint is accepted");
  require(passes_electron_crack_acceptance(crack_high),
          "upper electron-crack endpoint is accepted");
  require(!passes_electron_crack_acceptance(crack_inside),
          "electron strictly inside the crack is rejected");

  isolation_candidate.isolation_et = kIsolationEtMaxGeV;
  require(!passes_isolation(isolation_candidate),
          "isolation ET=5 GeV fails the open threshold");
  isolation_candidate.isolation_et =
      std::nextafter(kIsolationEtMaxGeV, 0.0);
  require(passes_isolation(isolation_candidate),
          "isolation immediately below 5 GeV passes");

  std::vector<Lepton> low_mass_boundary = {
      test_lepton(11, 0.5 * kAllOssfMassMinGeV, 0.0, 0.0, 0),
      test_lepton(-11, 0.5 * kAllOssfMassMinGeV, 0.0, kPi, 1),
      test_lepton(13,
                  back_to_back_pt_for_mass(20.0, kMuonMassGeV),
                  0.0, 1.0, 2),
      test_lepton(-13,
                  back_to_back_pt_for_mass(20.0, kMuonMassGeV),
                  0.0, 1.0 + kPi, 3),
  };
  const std::array<int, 4> indices = {{0, 1, 2, 3}};
  require_close(
      (low_mass_boundary[0].reconstructed +
       low_mass_boundary[1].reconstructed)
          .mass(),
      kAllOssfMassMinGeV, 1.0e-14,
      "low-mass boundary fixture has mll=4 GeV");
  require(!all_ossf_masses_pass(indices, low_mass_boundary),
          "OSSF pair with mll=4 GeV fails the open veto");
  low_mass_boundary[0] =
      test_lepton(11, 0.5 * kAllOssfMassMinGeV + 1.0e-6,
                  0.0, 0.0, 0);
  low_mass_boundary[1] =
      test_lepton(-11, 0.5 * kAllOssfMassMinGeV + 1.0e-6,
                  0.0, kPi, 1);
  require(all_ossf_masses_pass(indices, low_mass_boundary),
          "OSSF pair immediately above mll=4 GeV passes");

  const std::vector<int> subset = {0, 1, 2, 3};
  for (const double boundary : {kZ1MassMinGeV, kZ1MassMaxGeV}) {
    const std::vector<Lepton> leptons =
        mixed_flavour_quadruplet(boundary, 40.0);
    const SubsetResult result = analyze_subset(leptons, subset);
    require_close(result.z1.mass(), boundary, 1.0e-12,
                  "Z1 mass boundary fixture");
    require((result.cut_mask & kZ1Window) == 0U,
            "Z1 endpoint fails the open mass window");
  }
  for (const double inside :
       {kZ1MassMinGeV + 1.0e-6, kZ1MassMaxGeV - 1.0e-6}) {
    const SubsetResult result =
        analyze_subset(mixed_flavour_quadruplet(inside, 40.0), subset);
    require((result.cut_mask & kSelectedMask) == kSelectedMask,
            "Z1 immediately inside its mass window is selected");
  }

  for (const double boundary : {kZ2MassMinGeV, kZ2MassMaxGeV}) {
    const std::vector<Lepton> leptons =
        mixed_flavour_quadruplet(90.0, boundary);
    const SubsetResult result = analyze_subset(leptons, subset);
    require_close(result.z2.mass(), boundary, 1.0e-12,
                  "Z2 mass boundary fixture");
    require((result.cut_mask & kZ2Window) == 0U,
            "Z2 endpoint fails the open mass window");
  }
  for (const double inside :
       {kZ2MassMinGeV + 1.0e-6, kZ2MassMaxGeV - 1.0e-6}) {
    const SubsetResult result =
        analyze_subset(mixed_flavour_quadruplet(90.0, inside), subset);
    require((result.cut_mask & kSelectedMask) == kSelectedMask,
            "Z2 immediately inside its mass window is selected");
  }
}

void test_ossf_pairing_and_cutflow() {
  std::vector<Lepton> leptons = {
      test_lepton(11, 45.0, 0.0, 0.0, 0),
      test_lepton(-11, 45.0, 0.0, kPi, 1),
      test_lepton(13, 20.0, 0.0, 1.0, 2),
      test_lepton(-13, 20.0, 0.0, 1.0 + kPi, 3),
  };
  const std::vector<int> subset = {0, 1, 2, 3};
  const SubsetResult selected = analyze_subset(leptons, subset);
  require(selected.cut_mask == kSelectedMask,
          "2e2mu quadruplet passes the complete cut mask");
  require(selected.channel == 3, "2e2mu channel classification");
  require_close(selected.z1.mass(), 90.0, 1.0e-10,
                "Z1 is the pair nearest the Z mass");
  require_close(selected.z2.mass(),
                2.0 * std::sqrt(20.0 * 20.0 +
                                kMuonMassGeV * kMuonMassGeV),
                1.0e-10,
                "remaining OSSF pair is Z2");
  require(leptons[selected.output_indices[0]].charge == -1 &&
              leptons[selected.output_indices[1]].charge == 1,
          "Z1 pair has deterministic negative-positive ordering");

  std::vector<Lepton> same_sign = {
      test_lepton(11, 45.0, 0.0, 0.0, 0),
      test_lepton(11, 45.0, 0.0, kPi, 1),
      test_lepton(13, 20.0, 0.0, 1.0, 2),
      test_lepton(13, 20.0, 0.0, 1.0 + kPi, 3),
  };
  const SubsetResult no_pair = analyze_subset(same_sign, subset);
  require(no_pair.cut_mask ==
              (kAtLeastFour | kKinematicAcceptance |
               kElectronCrackAcceptance | kIsolation),
          "same-sign quadruplet stops at the OSSF stage");
}

void test_identical_flavour_pairing_and_extra_leptons() {
  std::vector<Lepton> ambiguous = {
      test_lepton(11, 45.0, 0.0, 0.0, 0),
      test_lepton(-11, 45.0, 0.0, kPi, 1),
      test_lepton(11, 20.0, 0.0, 1.0, 2),
      test_lepton(-11, 20.0, 0.0, 1.0 + kPi, 3),
  };
  const std::vector<int> four = {0, 1, 2, 3};
  const PairingChoice ambiguous_choice = best_pairing(ambiguous, four);
  require(ambiguous_choice.has_pairing,
          "ambiguous 4e topology has a disjoint OSSF pairing");
  require_close(ambiguous_choice.z1.mass(), 90.0, 1.0e-10,
                "ambiguous 4e topology chooses the pair nearest mZ");
  require(ambiguous_choice.pairing_source_indices ==
              std::array<int, 4>{{0, 1, 2, 3}},
          "ambiguous 4e pairing is deterministic");
  require(analyze_subset(ambiguous, four).channel == 1,
          "identical-electron topology is classified as 4e");

  std::vector<Lepton> higher_sum_pt = ambiguous;
  higher_sum_pt.push_back(test_lepton(11, 25.0, 0.0, 1.0, 4));
  const PairingChoice higher_sum_choice =
      best_pairing(higher_sum_pt, {0, 1, 2, 3, 4});
  require(higher_sum_choice.sorted_source_indices ==
              std::array<int, 4>{{0, 1, 3, 4}},
          "extra-lepton ambiguity prefers decreasing quadruplet sum pT");

  std::vector<Lepton> source_tie = ambiguous;
  source_tie.push_back(test_lepton(11, 20.0, 0.0, 1.0, 4));
  const PairingChoice source_tie_choice =
      best_pairing(source_tie, {0, 1, 2, 3, 4});
  require(source_tie_choice.sorted_source_indices ==
              std::array<int, 4>{{0, 1, 2, 3}},
          "exact extra-lepton tie prefers stable source indices");

  std::vector<Lepton> cross_pair_veto = {
      test_lepton(11, 45.0, 0.0, 0.0, 0),
      test_lepton(-11, 45.0, 0.0, kPi, 1),
      test_lepton(11, 20.0, 0.0, kPi + 0.01, 2),
      test_lepton(-11, 20.0, 0.0, 0.01, 3),
  };
  const SubsetResult vetoed = analyze_subset(cross_pair_veto, four);
  require((vetoed.cut_mask & kOssfPairing) != 0U,
          "cross-pair-veto fixture has two selected OSSF pairs");
  require((vetoed.cut_mask & kAllOssfMasses) == 0U,
          "all-OSSF veto checks cross-pairs in identical-flavour events");
}

void test_subset_probability_closure() {
  std::vector<Lepton> leptons = {
      test_lepton(11, 45.0, 0.0, 0.0, 0, 0.90),
      test_lepton(-11, 45.0, 0.0, kPi, 1, 0.90),
      test_lepton(13, 20.0, 0.0, 1.0, 2, 0.8075),
      test_lepton(-13, 20.0, 0.0, 1.0 + kPi, 3, 0.8075),
      test_lepton(11, 12.0, 2.0, 2.0, 4, 0.90),
  };
  const ULong64_t count = ULong64_t{1} << leptons.size();
  double total = 0.0;
  double lower_multiplicity = 0.0;
  double stored = 0.0;
  for (ULong64_t mask = 0; mask < count; ++mask) {
    const double probability = subset_probability(leptons, mask);
    total += probability;
    if (popcount(mask) < 4) {
      lower_multiplicity += probability;
    } else {
      stored += probability;
    }
  }
  require_close(total, 1.0, 1.0e-14,
                "exclusive reconstruction subsets sum to unity");
  require_close(lower_multiplicity + stored, 1.0, 1.0e-14,
                "stored plus response-failure probability closes");
  require_close(subset_probability(leptons, count - 1U),
                0.90 * 0.90 * 0.8075 * 0.8075 * 0.90, 1.0e-15,
                "all-reconstructed subset probability");
}

void test_inclusive_jet_context() {
  double jets[5][kMaxRecoObjects] = {};
  double bjets[5][kMaxRecoObjects] = {};
  double missing[4] = {};
  jets[0][0] = 25.0;
  jets[1][0] = 25.0;
  bjets[0][0] = 30.0;
  bjets[2][0] = 30.0;
  const RecoContext context =
      reconstructed_context(true, 1, jets, true, 1, bjets, false, missing);
  require(context.n_jets == 2,
          "inclusive jet count includes HwSim's separate b-jet collection");
  require(context.n_bjets == 1, "b-jet count retains the tagged subset");
  require_close(context.leading_jet_pt, 30.0, 1.0e-12,
                "leading inclusive jet can be a tagged b jet");
}

}  // namespace

int main() {
  try {
    test_muon_resolution_table();
    test_deterministic_smearing();
    test_dressing_and_isolation();
    test_open_boundaries();
    test_ossf_pairing_and_cutflow();
    test_identical_flavour_pairing_and_extra_leptons();
    test_subset_probability_closure();
    test_inclusive_jet_context();
    std::cout << "four-lepton core self-test passed (" << checks
              << " checks)" << std::endl;
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "four-lepton core self-test failed after " << checks
              << " checks: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
