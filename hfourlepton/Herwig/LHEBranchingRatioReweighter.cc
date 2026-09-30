// -*- C++ -*-
//
// Exact forced-decay importance weight for Les Houches event handlers.
//
// Herwig::BranchingRatioReweighter computes the desired selector-sum
// product, but applies it only when the active handler is a
// ThePEG::StandardEventHandler.  LHE shower runs use
// ThePEG::LesHouchesEventHandler, so these handlers apply the corresponding
// factors directly to the Event and to every optional event weight.
//
// Prompt W/Z selections and late heavy-flavour selections deliberately use
// separate handlers.  Prompt resonances are complete after the cascade but
// before heavy-hadron decays, while selected ground-state B mesons can be
// created recursively during the decay stage.  Keeping the scopes separate
// prevents virtual W/Z phase-space lines inside hadron decays from acquiring
// an unrelated prompt branching factor.

#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/EnumParticles.h"

#include <cmath>
#include <set>

namespace HiggsSSC {

using namespace ThePEG;

namespace {

enum class ReweightScope { Prompt, HeavyFlavor };
constexpr double kHeavyFlavorSemileptonicBias = 4.0;

bool in_scope(const tcPPtr& particle, ReweightScope scope) {
  const long absolute_id = std::abs(particle->id());
  if (scope == ReweightScope::Prompt) {
    return absolute_id == 23 || absolute_id == 24;
  }
  return absolute_id == 511 || absolute_id == 521 || absolute_id == 531;
}

bool is_terminal_copy(const tcPPtr& particle, ReweightScope scope) {
  const tcPDPtr charge_conjugate = particle->dataPtr()->CC();
  for (const auto& child : particle->children()) {
    if (child->id() == particle->id()) {
      return false;
    }
    // HwDecayHandler represents neutral B mixing by a single
    // charge-conjugate continuation.  It exists only in the late
    // heavy-flavour scope and must not contribute a second factor.
    if (scope == ReweightScope::HeavyFlavor &&
        particle->children().size() == 1 && charge_conjugate &&
        child->dataPtr() == charge_conjugate) {
      return false;
    }
  }
  return true;
}

bool is_biased_semileptonic_mode(const tcDMPtr& mode) {
  if (!mode) {
    return false;
  }
  int charged_leptons = 0;
  int electron_muon_neutrinos = 0;
  int other_products = 0;
  for (const tcPDPtr product : mode->orderedProducts()) {
    const long absolute_id = std::abs(product->id());
    if (absolute_id == 11 || absolute_id == 13) {
      ++charged_leptons;
    } else if (absolute_id == 12 || absolute_id == 14) {
      ++electron_muon_neutrinos;
    } else {
      ++other_products;
    }
  }
  return charged_leptons == 1 && electron_muon_neutrinos == 1 &&
         other_products == 1;
}

double heavy_flavour_importance_factor(const tcPPtr& particle,
                                       const char* handler_name) {
  // The campaign multiplies only the selected direct one-hadron e/mu
  // branching ratios by kHeavyFlavorSemileptonicBias.  All other active
  // modes retain b_m=1, so the proposal has full support.  Recover the
  // original active branching sum from the biased table and apply p_m/q_m:
  //
  //   p_m/q_m = [sum_n b_n p_n] / [sum_n p_n * b_m].
  const double biased_sum = particle->dataPtr()->decaySelector().sum();
  double original_sum = 0.0;
  for (const tcDMPtr mode : particle->dataPtr()->decayModes()) {
    if (!mode->on() || mode->brat() <= 0.0) {
      continue;
    }
    const double mode_bias =
        is_biased_semileptonic_mode(mode)
            ? kHeavyFlavorSemileptonicBias
            : 1.0;
    original_sum += mode->brat() / mode_bias;
  }
  const double selected_bias =
      is_biased_semileptonic_mode(particle->decayMode())
          ? kHeavyFlavorSemileptonicBias
          : 1.0;
  const double factor = biased_sum / (original_sum * selected_bias);
  if (!std::isfinite(factor) || factor <= 0.0) {
    throw Exception() << handler_name
                      << " obtained invalid heavy-flavour importance factor "
                      << factor << " (biased sum " << biased_sum
                      << ", original sum " << original_sum
                      << ", selected bias " << selected_bias << ")"
                      << Exception::runerror;
  }
  return factor;
}

void apply_reweighting(EventHandler& event_handler, ReweightScope scope,
                       const char* handler_name) {
  tEventPtr event = event_handler.currentEvent();
  if (!event) {
    return;
  }

  double factor = 1.0;
  std::set<tcPPtr> particles;
  event->select(std::inserter(particles, particles.end()), AllSelector());
  for (const tcPPtr& particle : particles) {
    if (!in_scope(particle, scope) || particle->dataPtr()->stable()) {
      continue;
    }
    if (particle->id() == ParticleID::Remnant ||
        particle->id() == ParticleID::Cluster) {
      continue;
    }
    if (particle->mass() < ZERO) {
      continue;
    }
    if (particle == event->incoming().first ||
        particle == event->incoming().second) {
      continue;
    }
    // Ground-state B factors are evaluated after all recursive decays.  A
    // non-null DecayMode distinguishes a selector-driven B decay from a
    // history or phase-space line.  Shower-reconstructed prompt W/Z copies
    // do not retain this marker, hence their separate earlier handler.
    if (scope == ReweightScope::HeavyFlavor && !particle->decayMode()) {
      continue;
    }
    if (!is_terminal_copy(particle, scope)) {
      continue;
    }
    factor *=
        scope == ReweightScope::HeavyFlavor
            ? heavy_flavour_importance_factor(particle, handler_name)
            : particle->dataPtr()->decaySelector().sum();
  }

  if (!std::isfinite(factor) || factor < 0.0) {
    throw Exception() << handler_name << " obtained invalid factor " << factor
                      << Exception::runerror;
  }

  event->weight(event->weight() * factor);
  for (auto& named_weight : event->optionalWeights()) {
    named_weight.second *= factor;
  }
}

}  // namespace

class LHEBranchingRatioReweighter final : public StepHandler {
 public:
  LHEBranchingRatioReweighter() = default;

  void handle(EventHandler& event_handler, const tPVector&, const Hint&) override {
    apply_reweighting(event_handler, ReweightScope::Prompt,
                      "LHEBranchingRatioReweighter");
  }

  static void Init() {
    static ClassDocumentation<LHEBranchingRatioReweighter> documentation(
        "Applies selected prompt W/Z branching-ratio factors to Les Houches "
        "event weights after the perturbative cascade.");
  }

 protected:
  IBPtr clone() const override { return new_ptr(*this); }
  IBPtr fullclone() const override { return new_ptr(*this); }
};

class LHEHeavyFlavorBranchingRatioReweighter final : public StepHandler {
 public:
  LHEHeavyFlavorBranchingRatioReweighter() = default;

  void handle(EventHandler& event_handler, const tPVector&, const Hint&) override {
    apply_reweighting(event_handler, ReweightScope::HeavyFlavor,
                      "LHEHeavyFlavorBranchingRatioReweighter");
  }

  static void Init() {
    static ClassDocumentation<LHEHeavyFlavorBranchingRatioReweighter>
        documentation(
            "Applies exact full-support ground-state B-meson decay-importance "
            "weights to Les Houches events after recursive hadron decays.");
  }

 protected:
  IBPtr clone() const override { return new_ptr(*this); }
  IBPtr fullclone() const override { return new_ptr(*this); }
};

DescribeNoPIOClass<LHEBranchingRatioReweighter, StepHandler>
    describeHiggsSSCLHEBranchingRatioReweighter(
        "HiggsSSC::LHEBranchingRatioReweighter",
        "LHEBranchingRatioReweighter.so");
DescribeNoPIOClass<LHEHeavyFlavorBranchingRatioReweighter, StepHandler>
    describeHiggsSSCLHEHeavyFlavorBranchingRatioReweighter(
        "HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter",
        "LHEBranchingRatioReweighter.so");

}  // namespace HiggsSSC
