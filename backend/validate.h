#pragma once

#include "index.h"
#include "features.h"
#include "../common/filesystem.h"
#include <set>
#include <map>

namespace inspire {
  namespace backend {
    class Validate {
      private:

      std::set<std::string> PROTEINS;
      std::set<std::string> PROTEINS_VALID;

      static inline void relative_complement(std::set<std::string> &from, std::set<std::string> what) {
        for (auto what_it = what.begin(); what_it != what.end(); ++what_it) {
          from.erase(*what_it);
        }
      }

      static inline void print_serialized(const std::set<std::string> &list) {
        for (auto it = list.begin(); it != list.end(); ++it) {
          std::cout << '\t' << *it;
        }
        std::cout << '\n';
      }

      bool contains_heteroatoms(const std::string &file, std::string &id) {
        id = ProteinParser::parse_id(file);
        return PROTEINS.find(id) != PROTEINS.end() && ProteinParser::contains_hetatoms(file);
      }

      public:
      Validate(std::string index_file) {
        Index index(index_file);
        if (index.reset()) {
          std::string protein_name = index.protein();
          PROTEINS.emplace(protein_name);
          do {
            if (protein_name != index.protein()) {
              protein_name = index.protein();
              PROTEINS.emplace(protein_name);
            }
          } while (index.next());
        }
        PROTEINS_VALID.insert(PROTEINS.begin(), PROTEINS.end());
      }

      void validate_heteroatoms(std::set<std::string> paths) {
        std::set<std::string> proteins_heteroatoms;

        std::string id;
        for (auto paths_it = paths.begin(); paths_it != paths.end(); ++paths_it) {
          if (common::filesystem::is_directory(*paths_it)) {
            common::filesystem::RecursiveDirectoryFileIterator file_iterator(*paths_it);
            if (file_iterator.has_file()) {
              do {
                try {
                  if (contains_heteroatoms(file_iterator.filename(), id)) {
                    proteins_heteroatoms.emplace(id);
                  }
                } catch (...) { }
              } while (file_iterator.has_next());
            } else {
              std::cerr << "There is no file in the given directory." << std::endl;
            }
          } else {
            try {
              if (contains_heteroatoms(*paths_it, id)) {
                proteins_heteroatoms.emplace(id);
              }
            } catch (...) { }
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Heteroatoms\t" << proteins_heteroatoms.size();
        print_serialized(proteins_heteroatoms);
        relative_complement(PROTEINS_VALID, proteins_heteroatoms);
      }

      void validate_remark_350(std::string index_log) {
        std::set<std::string> proteins_remark_350;
        
        std::string line;
        std::ifstream index_log_file(index_log);
        while (std::getline(index_log_file, line)) {
          if (line.find("Protein does not contain a definition of biomolecules (REMARK 350)") != line.npos) {
            proteins_remark_350.emplace(line.substr(9, 4));
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Remark 350 \t" << proteins_remark_350.size();
        print_serialized(proteins_remark_350);
        relative_complement(PROTEINS_VALID, proteins_remark_350);
        
      }

      void validate_chains_count(size_t min_chains, std::string index_file) {
        std::set<std::string> proteins_multimer;

        std::string protein_name;
        Index index(index_file);

        if (index.reset()) {
          protein_name = index.protein();
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t chains = 1;
          do {
            if (protein_name != index.protein()) {
              if (chains < min_chains) {
                proteins_multimer.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
            } else if (model_name != index.model()) {
              if (chains < min_chains) {
                proteins_multimer.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
            } else if (chain_name != index.chain()) {
              chain_name = index.chain();
              ++chains;
            }
          } while (index.next());
          if (chains < min_chains) {
            proteins_multimer.emplace(protein_name);
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Monomer    \t" << proteins_multimer.size();
        print_serialized(proteins_multimer);
        relative_complement(PROTEINS_VALID, proteins_multimer);
      }

      void validate_ambiguity(std::string index_file) {
        std::set<std::string> proteins_ambiguity;

        std::string protein_name;
        Index index(index_file);
        
        if (index.reset()) {
          protein_name = index.protein();
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t chains = 1;
          std::set<size_t> sizes;
          do {
            if (protein_name != index.protein()) {
              if (sizes.size() > 1 || (sizes.size() == 1 && chains != *sizes.begin())) {
                proteins_ambiguity.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
              sizes.clear();
            } else if (model_name != index.model()) {
              sizes.emplace(chains);
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
            } else if (chain_name != index.chain()) {
              chain_name = index.chain();
              ++chains;
            }
          } while (index.next());
          if (sizes.size() > 1 || (sizes.size() == 1 && chains != *sizes.begin())) {
            proteins_ambiguity.emplace(protein_name);
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Ambiguous  \t" << proteins_ambiguity.size();
        print_serialized(proteins_ambiguity);
        relative_complement(PROTEINS_VALID, proteins_ambiguity);
      }

      void validate_residues(std::string index_file, std::string aminoacids_file, std::string aminoacids_valid) {
        std::set<std::string> proteins_unknown;

        std::string protein_name;
        Index index(index_file);
        FeatureReader aminoacid(aminoacids_file, "aminoacid"); // NOTE: It is also possible to check aminoacid's log
        std::set<std::string> valid;

        {
          std::ifstream input(aminoacids_valid);
          std::string line;
          while (std::getline(input, line)) {
            size_t tab = line.find('\t');
            if (tab == line.npos || tab != line.rfind('\t')) {
              std::cerr << "Unexpected format of line in '" << aminoacids_valid << "': '" << line << "'" << std::endl;
            } else {
              valid.emplace(line.substr(0, tab));
            }
          }
        }

        if (index.reset()) {
          protein_name = index.protein();
          bool unknown = true;
          do {
            if (protein_name != index.protein()) {
              protein_name = index.protein();
              unknown = true;
            }
            if (!aminoacid.next(index.index()) || aminoacid.value().empty() || valid.find(aminoacid.value()) == valid.end()) {
              if (unknown) {
                proteins_unknown.emplace(index.protein());
                unknown = false;
              }
            }
          } while (index.next());
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Unknown    \t" << proteins_unknown.size();
        print_serialized(proteins_unknown);
        relative_complement(PROTEINS_VALID, proteins_unknown);
      }

      void validate_length(size_t min_length, std::string index_file, ProteinIterator* iterator, std::string index_log) {
        std::set<std::string> proteins_length;

        std::string protein_name;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>> missing;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>>::const_iterator missing_protein_it = missing.end();
        std::map<std::string, std::map<std::string, size_t>>::const_iterator missing_model_it;

        {
          std::string line;
          std::ifstream index_log_file(index_log);
          while (std::getline(index_log_file, line)) {
            size_t first;
            size_t second;
            size_t third;
            size_t fourth;
            size_t fifth;
            size_t sixth;
            if ((first = line.find("WARNING [")) != line.npos && (second = line.find("]:", first+9)) != line.npos && (third = line.find("Chain '", second+2)) != line.npos &&
              (fourth = line.find("' in model '", third+7)) != line.npos && (fifth = line.find("' does not contain information about ", fourth+12)) != line.npos &&
                (sixth = line.find(" residues", fifth+37)) != line.npos) {
              missing[line.substr(first+9, second-first-9)][line.substr(fourth+12, fifth-fourth-12)][line.substr(third+7, fourth-third-7)] = std::stol(line.substr(fifth+37, sixth-fifth-37));
            }
          }
        }

        Index index(index_file);
        
        if (index.reset()) {
          protein_name = index.protein();
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t count = 0;
          do {
            if (protein_name != index.protein()) {
              if (count < min_length) {
                proteins_length.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              count = 0;
              missing_protein_it = missing.find(protein_name);
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (model_name != index.model()) {
              if (count < min_length) {
                proteins_length.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              count = 0;
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (chain_name != index.chain()) {
              if (count < min_length) {
                proteins_length.emplace(protein_name);
              }
              chain_name = index.chain();
              count = 0;
              if (missing_protein_it != missing.end() && missing_model_it != missing_protein_it->second.end()) {
                auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                if (missing_chain_it != missing_model_it->second.end()) {
                  count += missing_chain_it->second;
                }
              }
            }
            ++count;
          } while (index.next());
          if (count < min_length) {
            proteins_length.emplace(protein_name);
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Length     \t" << proteins_length.size();
        print_serialized(proteins_length);
        relative_complement(PROTEINS_VALID, proteins_length);
      }

      void validate_incomplete_abs(size_t max_incomplete, std::string index_file, ProteinIterator* iterator, std::string index_log, std::string aminoacids_file, std::string elements_file, std::string elements_full) {
        std::set<std::string> proteins_abs_carbon_alpha;

        std::string protein_name;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>> missing;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>>::const_iterator missing_protein_it = missing.end();
        std::map<std::string, std::map<std::string, size_t>>::const_iterator missing_model_it;

        {
          std::string line;
          std::ifstream index_log_file(index_log);
          while (std::getline(index_log_file, line)) {
            size_t first;
            size_t second;
            size_t third;
            size_t fourth;
            size_t fifth;
            size_t sixth;
            if ((first = line.find("WARNING [")) != line.npos && (second = line.find("]:", first+9)) != line.npos && (third = line.find("Chain '", second+2)) != line.npos &&
              (fourth = line.find("' in model '", third+7)) != line.npos && (fifth = line.find("' does not contain information about ", fourth+12)) != line.npos &&
                (sixth = line.find(" residues", fifth+37)) != line.npos) {
              missing[line.substr(first+9, second-first-9)][line.substr(fourth+12, fifth-fourth-12)][line.substr(third+7, fourth-third-7)] = std::stol(line.substr(fifth+37, sixth-fifth-37));
            }
          }
        }

        Index index(index_file);
        FeatureReader aminoacid(aminoacids_file, "aminoacid");
        FeatureReader elements(elements_file, "composition");
        std::map<std::string, std::set<std::string>> full;

        {
          std::ifstream input(elements_full);
          std::string line;
          while (std::getline(input, line)) {
            std::string residue;
            std::stringstream parts(line);
            if (std::getline(parts, residue, '\t')) {
              std::string elements;
              std::set<std::string> &sets = full[residue];
              while (std::getline(parts, elements, '\t')) {
                sets.emplace(elements);
              }
            } else {
              std::cerr << "Unexpected format of line in '" << elements_full << "': '" << line << "'" << std::endl;
            }
          }
        }

        if (index.reset()) {
          protein_name = index.protein();
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t abs = 0;
          do {
            if (protein_name != index.protein()) {
              if (abs > max_incomplete) {
                proteins_abs_carbon_alpha.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              abs = 0;
              missing_protein_it = missing.find(protein_name);
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    abs += missing_chain_it->second;
                  }
                }
              }
            } else if (model_name != index.model()) {
              if (abs > max_incomplete) {
                proteins_abs_carbon_alpha.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              abs = 0;
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    abs += missing_chain_it->second;
                  }
                }
              }
            } else if (chain_name != index.chain()) {
              chain_name = index.chain();
              if (missing_protein_it != missing.end() && missing_model_it != missing_protein_it->second.end()) {
                auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                if (missing_chain_it != missing_model_it->second.end()) {
                  abs += missing_chain_it->second;
                }
              }
            }
            if (!aminoacid.next(index.index()) || aminoacid.value().empty()) {
              ++abs;
            } else if (elements.next(index.index())) {
              auto full_it = full.find(aminoacid.value());
              if (full_it == full.end()) {
                std::cerr << "Residue '" << aminoacid.value() << "' is not specified in file '" << elements_full << "'" << std::endl;
                ++abs;
              } else {
                bool remove = true;
                for (auto variants_it = full_it->second.begin(); variants_it != full_it->second.end(); ++variants_it) {
                  if (elements.value() == *variants_it) {
                    remove = false;
                    break;
                  }
                }
                if (remove) {
                  ++abs;
                }
              }
            }
          } while (index.next());
          if (abs > max_incomplete) {
            proteins_abs_carbon_alpha.emplace(protein_name);
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Abs C_alpha\t" << proteins_abs_carbon_alpha.size();
        print_serialized(proteins_abs_carbon_alpha);
        relative_complement(PROTEINS_VALID, proteins_abs_carbon_alpha);
      }

      void validate_incomplete_rel(double max_incomplete, std::string index_file, ProteinIterator* iterator, std::string index_log, std::string aminoacids_file, std::string elements_file, std::string elements_full) {
        std::set<std::string> proteins_rel_carbon_alpha;

        std::string protein_name;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>> missing;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>>::const_iterator missing_protein_it = missing.end();
        std::map<std::string, std::map<std::string, size_t>>::const_iterator missing_model_it;

        {
          std::string line;
          std::ifstream index_log_file(index_log);
          while (std::getline(index_log_file, line)) {
            size_t first;
            size_t second;
            size_t third;
            size_t fourth;
            size_t fifth;
            size_t sixth;
            if ((first = line.find("WARNING [")) != line.npos && (second = line.find("]:", first+9)) != line.npos && (third = line.find("Chain '", second+2)) != line.npos &&
              (fourth = line.find("' in model '", third+7)) != line.npos && (fifth = line.find("' does not contain information about ", fourth+12)) != line.npos &&
                (sixth = line.find(" residues", fifth+37)) != line.npos) {
              missing[line.substr(first+9, second-first-9)][line.substr(fourth+12, fifth-fourth-12)][line.substr(third+7, fourth-third-7)] = std::stol(line.substr(fifth+37, sixth-fifth-37));
            }
          }
        }

        Index index(index_file);
        FeatureReader aminoacid(aminoacids_file, "aminoacid");
        FeatureReader elements(elements_file, "composition");
        std::map<std::string, std::set<std::string>> full;

        {
          std::ifstream input(elements_full);
          std::string line;
          while (std::getline(input, line)) {
            std::string residue;
            std::stringstream parts(line);
            if (std::getline(parts, residue, '\t')) {
              std::string elements;
              std::set<std::string> &sets = full[residue];
              while (std::getline(parts, elements, '\t')) {
                sets.emplace(elements);
              }
            } else {
              std::cerr << "Unexpected format of line in '" << elements_full << "': '" << line << "'" << std::endl;
            }
          }
        }

        if (index.reset()) {
          protein_name = index.protein();
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t count = 0;
          size_t rel = 0;
          do {
            if (protein_name != index.protein()) {
              if (rel > max_incomplete*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              count = 0;
              rel = 0;
              missing_protein_it = missing.find(protein_name);
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    rel += missing_chain_it->second;
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (model_name != index.model()) {
              if (rel > max_incomplete*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              count = 0;
              rel = 0;
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    rel += missing_chain_it->second;
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (chain_name != index.chain()) {
              if (rel > max_incomplete*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              chain_name = index.chain();
              count = 0;
              rel = 0;
              if (missing_protein_it != missing.end() && missing_model_it != missing_protein_it->second.end()) {
                auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                if (missing_chain_it != missing_model_it->second.end()) {
                  rel += missing_chain_it->second;
                  count += missing_chain_it->second;
                }
              }
            }
            ++count;
            if (!aminoacid.next(index.index()) || aminoacid.value().empty()) {
              ++rel;
            } else if (elements.next(index.index())) {
              auto full_it = full.find(aminoacid.value());
              if (full_it == full.end()) {
                std::cerr << "Residue '" << aminoacid.value() << "' is not specified in file '" << elements_full << "'" << std::endl;
                ++rel;
              } else {
                bool remove = true;
                for (auto variants_it = full_it->second.begin(); variants_it != full_it->second.end(); ++variants_it) {
                  if (elements.value() == *variants_it) {
                    remove = false;
                    break;
                  }
                }
                if (remove) {
                  ++rel;
                }
              }
            }
          } while (index.next());
          if (rel > max_incomplete*count) {
            proteins_rel_carbon_alpha.emplace(protein_name);
          }
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Rel C_alpha\t" << proteins_rel_carbon_alpha.size();
        print_serialized(proteins_rel_carbon_alpha);
        relative_complement(PROTEINS_VALID, proteins_rel_carbon_alpha);
      }

      void validate_interfaces(std::string index_file, std::string interfaces_file, size_t interacting, size_t noninteracting) {
        std::set<std::string> proteins_monotonous;

        std::string protein_name;
        std::string model_name;
        std::string chain_name;
        Index index(index_file);
        FeatureReader interfaces(interfaces_file, "interface");

        if (index.reset()) {
          protein_name = index.protein();
          model_name = index.model();
          chain_name = index.chain();
          size_t length = 0;
          size_t count = 0;
          do {
            if (protein_name != index.protein()) {
              if (count < interacting) {
                proteins_monotonous.emplace(protein_name);
              }
              if (length-count < noninteracting) {
                proteins_monotonous.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              length = 0;
              count = 0;
            }
            if (model_name != index.model()) {
              if (count < interacting) {
                proteins_monotonous.emplace(protein_name);
              }
              if (length-count < noninteracting) {
                proteins_monotonous.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              length = 0;
              count = 0;
            }
            if (chain_name != index.chain()) {
              if (count < interacting) {
                proteins_monotonous.emplace(protein_name);
              }
              if (length-count < noninteracting) {
                proteins_monotonous.emplace(protein_name);
              }
              chain_name = index.chain();
              length = 0;
              count = 0;
            }
            ++length;
            if (interfaces.next(index.index()) && interfaces.value() == "I") {
              ++count;
            }
          } while (index.next());
        }

        std::cout << "Total      \t" << PROTEINS_VALID.size() << std::endl;
        std::cout << "Interfaces \t" << proteins_monotonous.size();
        print_serialized(proteins_monotonous);
        relative_complement(PROTEINS_VALID, proteins_monotonous);
      }



      void print_proteins() {
        std::cout << "Total      \t" << PROTEINS_VALID.size();
        for (auto proteins_it = PROTEINS_VALID.begin(); proteins_it != PROTEINS_VALID.end(); ++proteins_it) {
          std::cout << '\t' << *proteins_it;
        }
        std::cout << std::endl;
      }
    };
  }
}