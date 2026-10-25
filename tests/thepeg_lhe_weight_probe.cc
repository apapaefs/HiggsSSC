// Read LHE weight maps with the installed ThePEG parser. No event generation,
// showering, repository save, or physics initialization is performed.
#include "ThePEG/Interface/InterfaceBase.h"
#include "ThePEG/LesHouches/LesHouchesFileReader.h"
#include "ThePEG/Repository/Repository.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
  try {
    if (argc != 3) throw std::runtime_error("usage: thepeg_lhe_weight_probe LesHouches.so FILE.lhe[.gz]");
    // Load the plugin as Herwig does (also works with macOS loadable bundles).
    for (const auto& command : {std::string("library ") + argv[1],
                               std::string("create ThePEG::LesHouchesFileReader /LHEWeightProbe")}) {
      const auto error = ThePEG::Repository::exec(command, std::cerr);
      if (!error.empty()) throw std::runtime_error(error);
    }
    auto object = ThePEG::BaseRepository::GetPointer("/LHEWeightProbe");
    if (!object) throw std::runtime_error("ThePEG reader creation failed");
    // The create command above fixes the concrete type without linking a plugin.
    auto& reader = static_cast<ThePEG::LesHouchesFileReader&>(*object);
    const auto* filename = ThePEG::BaseRepository::FindInterface(object, "FileName");
    if (!filename) throw std::runtime_error("ThePEG FileName interface is unavailable");
    filename->exec(reader, "set", argv[2]);
    reader.open();
    std::cout << std::setprecision(17);
    long index = 0;
    while (reader.doReadEvent()) {
      std::cout << "event\t" << index++ << '\n';
      for (const auto& weight : reader.optionalEventWeights()) {
        std::cout << "weight\t" << weight.first << '\t' << weight.second << '\n';
      }
    }
    reader.close();
    if (!index) throw std::runtime_error("ThePEG read no LHE events");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
