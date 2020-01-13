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
      static const size_t MIN_LENGTH = 20;
      static const size_t MAX_INCOMPLETE_ABSOLUTE = 50;
      static const size_t MAX_INCOMPLETE_RELATIVE_NUMERATOR = 5;
      static const size_t MAX_INCOMPLETE_RELATIVE_DENOMINATOR = 100;
      static const size_t MIN_CHAINS = 2;

      std::set<std::string> files;

      static inline void relative_complement(std::map<std::string, std::string> &from, const std::set<std::string> what) {
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

      public:
      void add_files(std::string path) {
        if (common::filesystem::is_directory(path)) {
          common::filesystem::RecursiveDirectoryFileIterator file_iterator(path);
          if (file_iterator.has_file()) {
            do {
              files.emplace(file_iterator.filename());
            } while (file_iterator.has_next());
          } else {
            std::cerr << "There is no file in the given directory." << std::endl;
          }
        } else {
          files.emplace(path);
        }
      }

      void validate(std::string index_file, ProteinIterator* iterator, std::string index_log, std::string aminoacids_file, std::string aminoacids_valid, std::string elements_file, std::string elements_full) {
        std::set<std::string> proteins;
        std::set<std::string> proteins_remark_350;
        std::set<std::string> proteins_multimer; /*Some biomolecules could be valid*/
        std::set<std::string> proteins_unknown; /*Some biomolecules could be valid and some invalid?*/
        std::set<std::string> proteins_heteroatoms; /**/
        std::set<std::string> proteins_length; /*Some biomolecules/models could be valid?*/
        std::set<std::string> proteins_abs_carbon_alpha; /*Some biomolecules/models could be valid*/
        std::set<std::string> proteins_rel_carbon_alpha; /*Some biomolecules/models could be valid*/
        std::set<std::string> proteins_ambiguity; /*Some biomolecules/models could be valid*/

        std::string protein_name;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>> missing;
        std::map<std::string, std::map<std::string, std::map<std::string, size_t>>>::const_iterator missing_protein_it = missing.end();
        std::map<std::string, std::map<std::string, size_t>>::const_iterator missing_model_it;

        {
          std::string line;
          std::ifstream index_log_file(index_log);
          while (std::getline(index_log_file, line)) {
            if (line.find("Protein does not contain a definition of biomolecules (REMARK 350)") != line.npos) {
              proteins_remark_350.emplace(line.substr(9, 4));
            } else {
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
        }

        Index index(index_file);
        FeatureReader aminoacid(aminoacids_file, "aminoacid"); // NOTE: It is also possible to check aminoacid's log
        std::set<std::string> valid;
        FeatureReader elements(elements_file, "composition");
        std::map<std::string, std::set<std::string>> full;

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
          proteins.emplace(protein_name);
          std::string model_name = index.model();
          std::string chain_name = index.chain();
          size_t chains = 1;
          bool unknown = true;
          size_t count = 0;
          size_t abs = 0;
          size_t rel = 0;
          std::set<size_t> sizes;
          size_t size = 1;
          do {
            if (protein_name != index.protein()) {
              if (chains < MIN_CHAINS) {
                proteins_multimer.emplace(protein_name);
              }
              if (count < MIN_LENGTH) {
                proteins_length.emplace(protein_name);
              }
              if (abs >= MAX_INCOMPLETE_ABSOLUTE) {
                proteins_abs_carbon_alpha.emplace(protein_name);
              }
              if (MAX_INCOMPLETE_RELATIVE_DENOMINATOR*rel >= MAX_INCOMPLETE_RELATIVE_NUMERATOR*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              if (sizes.size() > 1 || (sizes.size() == 1 && size != *sizes.begin())) {
                proteins_ambiguity.emplace(protein_name);
              }
              protein_name = index.protein();
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
              unknown = true;
              count = 0;
              abs = 0;
              rel = 0;
              sizes.clear();
              size = 1;
              proteins.emplace(protein_name);
              missing_protein_it = missing.find(protein_name);
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    abs += missing_chain_it->second;
                    rel += missing_chain_it->second;
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (model_name != index.model()) {
              sizes.emplace(size);
              if (chains < MIN_CHAINS) {
                proteins_multimer.emplace(protein_name);
              }
              if (count < MIN_LENGTH) {
                proteins_length.emplace(protein_name);
              }
              if (MAX_INCOMPLETE_RELATIVE_DENOMINATOR*rel >= MAX_INCOMPLETE_RELATIVE_NUMERATOR*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              model_name = index.model();
              chain_name = index.chain();
              chains = 1;
              count = 0;
              rel = 0;
              size = 1;
              if (missing_protein_it != missing.end()) {
                missing_model_it = missing_protein_it->second.find(iterator->getBasicModelName(model_name));
                if (missing_model_it != missing_protein_it->second.end()) {
                  auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                  if (missing_chain_it != missing_model_it->second.end()) {
                    abs += missing_chain_it->second;
                    rel += missing_chain_it->second;
                    count += missing_chain_it->second;
                  }
                }
              }
            } else if (chain_name != index.chain()) {
              if (count < MIN_LENGTH) {
                proteins_length.emplace(protein_name);
              }
              if (MAX_INCOMPLETE_RELATIVE_DENOMINATOR*rel >= MAX_INCOMPLETE_RELATIVE_NUMERATOR*count) {
                proteins_rel_carbon_alpha.emplace(protein_name);
              }
              chain_name = index.chain();
              ++chains;
              count = 0;
              rel = 0;
              ++size;
              if (missing_protein_it != missing.end() && missing_model_it != missing_protein_it->second.end()) {
                auto missing_chain_it = missing_model_it->second.find(iterator->getBasicChainName(chain_name));
                if (missing_chain_it != missing_model_it->second.end()) {
                  abs += missing_chain_it->second;
                  rel += missing_chain_it->second;
                  count += missing_chain_it->second;
                }
              }
            }
            ++count;
            if (!aminoacid.next(index.index()) || aminoacid.value().empty() || valid.find(aminoacid.value()) == valid.end()) {
              if (unknown) {
                proteins_unknown.emplace(index.protein());
                unknown = false;
              }
            } else if (elements.next(index.index())) {
              auto full_it = full.find(aminoacid.value());
              if (full_it == full.end()) {
                std::cerr << "Residue '" << aminoacid.value() << "' is not specified in file '" << elements_full << "'" << std::endl;
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
                  ++rel;
                }
              }
            }
          } while (index.next());
          if (chains < MIN_CHAINS) {
            proteins_multimer.emplace(protein_name);
          }
          if (count < MIN_LENGTH) {
            proteins_length.emplace(protein_name);
          }
          if (abs >= MAX_INCOMPLETE_ABSOLUTE) {
            proteins_abs_carbon_alpha.emplace(protein_name);
          }
          if (MAX_INCOMPLETE_RELATIVE_DENOMINATOR*rel >= MAX_INCOMPLETE_RELATIVE_NUMERATOR*count) {
            proteins_rel_carbon_alpha.emplace(protein_name);
          }
          if (sizes.size() > 1 || (sizes.size() == 1 && size != *sizes.begin())) {
            proteins_ambiguity.emplace(protein_name);
          }
        }

        std::map<std::string, std::string> paths;
        for (auto files_it = files.begin(); files_it != files.end(); ++files_it) {
          try {
            std::string id = ProteinParser::parse_id(*files_it);
            if (proteins.find(id) != proteins.end()) {
              paths[id] = *files_it;
              if (ProteinParser::contains_hetatoms(*files_it)) {
                proteins_heteroatoms.emplace(id);
              }
            }
          } catch (...) { }
        }

        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Rel C_alpha\t" << proteins_rel_carbon_alpha.size();
        print_serialized(proteins_rel_carbon_alpha);
        relative_complement(paths, proteins_rel_carbon_alpha);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Abs C_alpha\t" << proteins_abs_carbon_alpha.size();
        print_serialized(proteins_abs_carbon_alpha);
        relative_complement(paths, proteins_abs_carbon_alpha);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Heteroatoms\t" << proteins_heteroatoms.size();
        print_serialized(proteins_heteroatoms);
        //relative_complement(paths, proteins_heteroatoms);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Unknown    \t" << proteins_unknown.size();
        print_serialized(proteins_unknown);
        relative_complement(paths, proteins_unknown);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Length     \t" << proteins_length.size();
        print_serialized(proteins_length);
        relative_complement(paths, proteins_length);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Ambiguous  \t" << proteins_ambiguity.size();
        print_serialized(proteins_ambiguity);
        relative_complement(paths, proteins_ambiguity);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Monomer    \t" << proteins_multimer.size();
        print_serialized(proteins_multimer);
        relative_complement(paths, proteins_multimer);
        std::cout << "Total      \t" << paths.size() << std::endl;
        std::cout << "Remark 350 \t" << proteins_remark_350.size();
        print_serialized(proteins_remark_350);
        relative_complement(paths, proteins_remark_350);
        std::cout << "Total      \t" << paths.size();
        for (auto paths_it = paths.begin(); paths_it != paths.end(); ++paths_it) {
          std::cout << '\t' << paths_it->first;
        }
        std::cout << std::endl;
      }
    };
  }
}