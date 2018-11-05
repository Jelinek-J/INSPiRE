#pragma once

#include "index.h"
#include <map>
#include <unordered_map>

namespace inspire {
  namespace backend {
    class Bimap {
      private:
      std::vector<std::string> names;
      std::unordered_map<std::string, size_t> indices;

      public:
      size_t add_or_get(const std::string &name) {
        auto ret = indices.insert({name, names.size()});
        if (ret.second) {
          names.push_back(name);
        }
        return ret.first->second;
      }

      std::string get_name(const size_t index) {
        return names[index];
      }

      size_t get_index(const std::string name) {
        indices[name];
      }
    };

    class Combinator {
      public:
      virtual double combine(double query, double temp) = 0;
    };
    class MinimalCombinator : public Combinator {
      public:
      double combine(double query, double temp) override {
        return std::min(query, temp);
      }
    };
    class MaximalCombinator : public Combinator {
      public:
      double combine(double query, double temp) override {
        return std::max(query, temp);
      }
    };

    class Similariter {
      protected:
      virtual std::pair<std::string, std::string> label(std::string protein, std::string model, std::string chain, std::string residue) = 0;

      private:
      // Dictionary of proteins: name <=> index
      Bimap KB_PROTEIN_NAMES;
      // Dictionary of chains: name <=> index
      Bimap KB_RESIDUE_NAMES;
      // What residue and protein corresponds to the given index
      std::vector<std::pair<size_t, size_t>> KB_PROTEINS;
      // How many residues is in the given protein
      std::map<size_t, size_t> KB_SIZES;
      // Combines scores for query and template
      Combinator* COMBINATOR;

      inline void initialize(const std::string &filename, Bimap &protein_names, Bimap &residue_names, std::vector<std::pair<size_t, size_t>> &proteins, std::map<size_t, size_t> &sizes) {
        Index index(filename);
        if (!index.reset()) {
          throw elemental::exception::TitledException("The index file is empty or it is not possible to read it.");
        }

        do {
          std::pair<std::string, std::string> key = label(index.protein(), index.model(), index.chain(), index.aminoacid());
          size_t protein_index = protein_names.add_or_get(key.first);
          size_t residue_index = residue_names.add_or_get(key.second);
          proteins.push_back(std::make_pair(protein_index, residue_index));
          if (proteins.size() != index.index()) {
            throw elemental::exception::TitledException("Unexpected format of index file: expected continuous aritmetic sequence starting at 1 with step 1");
          }
          auto sizes_it = sizes.insert({protein_index, 1});
          if (!sizes_it.second) {
            ++(sizes_it.first->second);
          }
        } while (index.next());
      }

      protected:
      inline void initialize(const std::string &filename) {
        initialize(filename, KB_PROTEIN_NAMES, KB_RESIDUE_NAMES, KB_PROTEINS, KB_SIZES);
      }

      public:
      Similariter(Combinator *combinator) : COMBINATOR(combinator) { }

      void analyze(std::string query_index, std::string mined, std::string output) {
        // Dictionary of proteins: name <=> index
        Bimap query_protein_names;
        // Dictionary of chains: name <=> index
        Bimap query_residue_names;
        // What residue and protein corresponds to the given index
        std::vector<std::pair<size_t, size_t>> query_proteins;
        // How many residues is in the given protein
        std::map<size_t, size_t> query_sizes;

        initialize(query_index, query_protein_names, query_residue_names, query_proteins, query_sizes);

        // Mined protein_id => residue_id => distance => residues
        std::map<size_t, std::map<size_t, std::map<size_t, std::vector<size_t>>>> similar;
        {
          std::ifstream input(mined);
          std::string line;
          while (!input.eof() && std::getline(input, line)) {
            size_t id = std::stoll(line)-1;
            size_t protein_id = query_proteins[id].first;
            std::map<size_t, std::vector<size_t>> &residue = similar[protein_id][id];
            while (!input.eof() && std::getline(input, line) && !line.empty()) {
              size_t index = line.find('\t');
              if (index == line.npos) {
                std::cerr << "Unexpected format of a line: '" << line << "'";
              } else {
                residue[std::stoll(line.substr(index + 1))].push_back(std::stoll(line.substr(0, index))-1);
              }
            }
          }
        }

        // Output stream
        if (output.empty() || output.back() == elemental::filesystem::directory_separator) {
          size_t i = mined.rfind(elemental::filesystem::directory_separator);
          std::string tmp = (i == mined.npos ? mined : mined.substr(i+1));
          if (elemental::string::ends_with(tmp, ".med")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".rty";
        } else if (!elemental::string::ends_with(output, ".rty")) {
          output += ".rty";
        }
        std::ofstream stream(output);

        // Iterate over query proteins
        for (auto query_sizes_it = query_sizes.begin(); query_sizes_it != query_sizes.end(); ++query_sizes_it) {
          stream << query_protein_names.get_name(query_sizes_it->first) << '\n';
          // Maximal score from the query protein view
          double max = 0;
          // residue_id => distance => residues
          auto &subsimilar = similar[query_sizes_it->first];
          // How many hits have individual proteins: protein_id => ratio;
          std::map<size_t, std::pair<double, std::map<size_t, double>>> counts;
          // Iterate over residues
          for (auto protein_similar_it = subsimilar.begin(); protein_similar_it != subsimilar.end(); ++protein_similar_it) {
            std::map<size_t, std::vector<size_t>> different;
            for (auto residue_similar_it = protein_similar_it->second.begin(); residue_similar_it != protein_similar_it->second.end(); ++residue_similar_it) {
              for (auto level_similar_it = residue_similar_it->second.begin(); level_similar_it != residue_similar_it->second.end(); ++level_similar_it) {
                different[KB_PROTEINS[*level_similar_it].first].push_back(*level_similar_it);
              }
            }
            if (different.size() > 0) {
              double weight = 1.0/different.size();
              max += weight;
              for (auto different_it = different.begin(); different_it != different.end(); ++different_it) {
                auto ins = counts.insert({different_it->first, {weight, std::map<size_t, double>()}});
                if (!ins.second) {
                  ins.first->second.first += weight;
                }
                for (auto same_it = different_it->second.begin(); same_it != different_it->second.end(); ++same_it) {
                  auto subins = ins.first->second.second.insert({*same_it, weight});
                  if (!subins.second && subins.first->second < weight) {
                    subins.first->second = weight;
                  }
                }
              }
            }
          }
          std::set<std::pair<double, size_t>> sorted;
          for (auto counts_it = counts.begin(); counts_it != counts.end(); ++counts_it) {
            double sum = 0;
            for (auto subcounts_it = counts_it->second.second.begin(); subcounts_it != counts_it->second.second.end(); ++subcounts_it) {
              sum += subcounts_it->second;
            }
            sorted.emplace(std::make_pair(COMBINATOR->combine(counts_it->second.first/max, sum/(sum+max/subsimilar.size()*(KB_SIZES[counts_it->first]-counts_it->second.second.size()))), counts_it->first));
          }
          for (auto sorted_rit = sorted.rbegin(); sorted_rit != sorted.rend(); ++sorted_rit) {
            stream << sorted_rit->first << '\t' << KB_PROTEIN_NAMES.get_name(sorted_rit->second) <<  '\n';
          }
          stream << std::endl;
        }
      }
    };

    class ProteinSimilariter : public Similariter {
      protected:
      std::pair<std::string, std::string> label(std::string protein, std::string model, std::string chain, std::string residue) override {
        if (model.size() == 0) {
          return std::make_pair(protein, chain + "\t" + residue);
        } else {
          return std::make_pair(protein + "\t" + model, chain + "\t" + residue);
        }
      }

      public:
      ProteinSimilariter(Combinator *combinator, std::string index) : Similariter(combinator) {
        initialize(index);
      }

    };

    class ChainSimilariter : public Similariter {
      protected:
      std::pair<std::string, std::string> label(std::string protein, std::string model, std::string chain, std::string residue) override {
        if (model.size() == 0) {
          return std::make_pair(protein + "\t" + chain, residue);
        } else {
          return std::make_pair(protein + "\t" + chain + "\t" + model, residue);
        }
      }

      public:
      ChainSimilariter(Combinator *combinator, std::string index) : Similariter(combinator) {
        initialize(index);
      }
    };
  }
}