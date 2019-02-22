#pragma once

#include "../backend/features.h"
#include "../backend/protein.h"
#include "freesasa.h"

namespace common {
  namespace sasa {
    class SasaFeature : public inspire::backend::Feature<float> {
      private:
      // Atom_name => [Residue_name => radius]
      std::map<std::string, std::map<std::string, float>> DISTANCES;
      // Model_name => [Chain_name => [Aminoacid_name => interface]]
      std::map<std::string, std::map<std::string, std::map<std::string, float>>> SASAS;

      public:
      SasaFeature(inspire::backend::ProteinIterator* iterator, std::string radiuses, std::string composition) : inspire::backend::Feature<float>(iterator) {
        std::map<std::string, float> elements;
        {
          std::ifstream input(radiuses);
          std::string line;
          while (!input.eof() && std::getline(input, line)) {
            size_t tab = line.find('\t');
            if (tab == line.npos) {
              std::cerr << "File with van der Waals radiuses '" << radiuses << "' contains invalid line '" << line << "'" << std::endl;
            } else {
              auto ins = elements.insert({line.substr(0, tab), std::stof(line.substr(tab+1))});
              if (!ins.second) {
                std::cerr << "File with van der Waals radiuses '" << radiuses << "' contains multiple definitions of element '" << line.substr(0, tab) << "'" << std::endl;
              }
            }
          }
          input.close();
        }
        {
          std::ifstream input(composition);
          std::string line;
          while (!input.eof() && std::getline(input, line)) {
            size_t first = line.find_first_of("\t ");
            if (first == line.npos) {
              std::cerr << "File with composition '" << composition << "' contains invalid line '" << line << "'" << std::endl;
            } else {
              size_t second = line.find_last_of("\t ");
              if (second == line.npos || second == first) {
                std::cerr << "File with composition radiuses '" << composition << "' contains invalid line '" << line << "'" << std::endl;
              } else {
                std::string key = line.substr(0, first);
                if (key == "ANY") {
                  key = "";
                }
                auto elements_it = elements.find(line.substr(second+1));
                if (elements_it == elements.end()) {
                  std::cerr << "File with van der Waals radiuses '" << radiuses << "does not contains a definition of '" << line.substr(second+1) << "' required in file with composition '" << composition << "'" << std::endl;
                } else {
                  DISTANCES[line.substr(first+1, second-first-1)][key] = elements_it->second;
                }
              }
            }
          }
          input.close();
        }
      }

      std::string title() override { return "sasa"; }

      void init(inspire::backend::Protein* protein) override {
        SASAS.clear();

        std::set<std::string> warned;
        ITERATOR->init(protein);
        if (ITERATOR->resetModel()) {
          do {
            auto &model = SASAS[ITERATOR->getModelName()];
            if (ITERATOR->resetChain()) {
              do {
                auto &chain = model[ITERATOR->getChainName()];
                std::map<std::string, std::vector<size_t>> indices;
                std::vector<double> coordinates;
                std::vector<double> radiuses;
                if (ITERATOR->resetAminoacid()) {
                  do {
                    auto &aminoacid = indices[ITERATOR->getAminoacidName()];
                    if (ITERATOR->resetAtom()) {
                      do {
                        auto distances_it = DISTANCES.find(ITERATOR->atom_name());
                        if (distances_it == DISTANCES.end()) {
                          if (warned.insert(ITERATOR->atom_name()).second) {
                            std::cerr << ITERATOR->getProteinName() << ": no radius is defined for atom '" << ITERATOR->atom_name() << "' (" << ITERATOR->element() << ")" << std::endl;
                          }
                        } else {
                          auto radius = distances_it->second.find(ITERATOR->residue_name());
                          if (radius == distances_it->second.end()) {
                            radius = distances_it->second.find("");
                          }
                          if (radius == distances_it->second.end()) {
                            if (warned.insert(ITERATOR->residue_name() + " " + ITERATOR->atom_name()).second) {
                              std::cerr << ITERATOR->getProteinName() << ": radius of atom '" << ITERATOR->atom_name() << "' (" << ITERATOR->element() << ") is not defined for residue '" << ITERATOR->residue_name() << "' nor default residue" << std::endl;
                            }
                          } else if (ITERATOR->computeCharacteristics()) {
                            aminoacid.push_back(radiuses.size());
                            coordinates.push_back(ITERATOR->x());
                            coordinates.push_back(ITERATOR->y());
                            coordinates.push_back(ITERATOR->z());
                            radiuses.push_back(radius->second);
                          } else {
                            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not contains a valid characteristic for an atom '" << ITERATOR->getAtomName() << "' (" << ITERATOR->element() << ") in aminoacid '"
                              << ITERATOR->getAminoacidName() << "' in chain '" << ITERATOR->getChainName() << "' in model '" << ITERATOR->getModelName() << "'" << std::endl;
                          }
                        }
                      } while (ITERATOR->nextAtom());
                    }
                  } while (ITERATOR->nextAminoacid());
                }
                if (radiuses.empty()) {
                  for (auto indices_it = indices.begin(); indices_it != indices.end(); ++indices_it) {
                    chain[indices_it->first] = UNDEFINED;
                  }
                } else {
                  freesasa_result *result = freesasa_calc_coord(&coordinates[0], &radiuses[0], radiuses.size(), nullptr);
                  for (auto indices_it = indices.begin(); indices_it != indices.end(); ++indices_it) {
                    if (indices_it->second.empty()) {
                      chain[indices_it->first] = UNDEFINED;
                    } else {
                      float &sasa = chain[indices_it->first];
                      sasa = 0;
                      for (auto atom = indices_it->second.begin(); atom != indices_it->second.end(); ++atom) {
                        sasa += result->sasa[*atom];
                      }
                    }
                  }
                  freesasa_result_free(result);
                }
              } while (ITERATOR->nextChain());
            }
          } while (ITERATOR->nextModel());
        }
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      float feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        auto sasas_it = SASAS.find(model);
        if (sasas_it == SASAS.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto chain_it = sasas_it->second.find(chain);
        if (chain_it == sasas_it->second.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto aminoacid_it = chain_it->second.find(aminoacid);
        if (aminoacid_it == chain_it->second.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" + aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        return aminoacid_it->second;
      }
    };
  }
}
