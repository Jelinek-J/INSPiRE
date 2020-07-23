#pragma once

#include "index.h"
#include "iterators.h"
#include "features.h"
#include "../common/exception.h"
#include "../common/random.h"
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>

namespace inspire {
  namespace backend {
    class Random {
      private:
      std::ofstream OUTPUT;

      void gselect(size_t count, double threshold, std::string header, ProteinIterator* iterator, std::string index_path, std::string neighbours_path, std::string similarity_path) {
        std::map<std::string, std::set<std::string>> collisions;
        std::ifstream similarities(similarity_path);
        std::string line;
        while (std::getline(similarities, line)) {
          std::set<std::string> &id = collisions[line.substr(0, 4)];
          while (std::getline(similarities, line) && !line.empty()) {
            size_t first = line.find('\t');
            if (first == line.npos) {
              throw common::exception::TitledException("Unexpected format of similarity file: '" + line + "'");
            }
            size_t second = line.rfind('\t');
            if (first == second) {
              throw common::exception::TitledException("Unexpected format of similarity file: '" + line + "'");
            }
            if (std::stod(line.substr(0, first)) >= threshold) {
              id.emplace(line.substr(first+1, second-first-1));
            }
          }
        }

        std::vector<std::string> order;
        std::map<std::string, std::string> proteins;
        inspire::backend::Index index(index_path);
        if (index.reset()) {
          std::string prev_protein;
          std::string prev_chain;
          std::vector<std::string> chains;
          bool add = false;
          do {
            std::string new_protein = index.protein();
            if (prev_protein != new_protein) {
              if (!chains.empty()) {
                proteins[prev_protein] = chains[common::random::random_size_t(chains.size())];
                order.push_back(prev_protein);
                chains.clear();
              }
              prev_protein = new_protein;
              prev_chain = "";
              add = collisions.find(new_protein) != collisions.end();
            }
            if (add && index.model().empty()) {
              std::string new_chain = index.chain();
              if (prev_chain != new_chain && iterator->getBasicChainName(new_chain) == new_chain) {
                chains.push_back(new_chain);
                prev_chain = new_chain;
              }
            }
          } while (index.next());
          if (!chains.empty()) {
            proteins[prev_protein] = chains[common::random::random_size_t(chains.size())];
            order.push_back(prev_protein);
          }
        }
        std::random_shuffle(order.begin(), order.end());

        std::set<std::string> selected;
        std::set<std::string> forbidden;
        while (selected.size() < count && selected.size() + forbidden.size() < order.size()) {
          std::string next = order[common::random::random_size_t(order.size())];
          if (selected.find(next) == selected.end() && forbidden.find(next) == forbidden.end()) {
            std::set<std::string> &collision = collisions[next];
            bool add = true;
            for (auto collision_it = collision.begin(); collision_it != collision.end(); ++collision_it) {
              if (selected.find(*collision_it) != selected.end()) {
                add = false;
                break;
              }
            }
            if (add) {
              selected.emplace(next);
              forbidden.insert(collision.begin(), collision.end());
            } else {
              forbidden.emplace(next);
            }
          }
        }
        if (selected.size() < count) {
          std::cerr << "Not enough unrelated structures, only " << selected.size() << " structures taken" << std::endl;
        }

        inspire::backend::FeatureReader neighbours(neighbours_path, "neighbours");
        // ProteinId => { chainId => { neighbourChainId => count }}
        std::map<std::string, std::map<std::string, std::map<std::string, double>>> pairs;
        if (index.reset()) {
          std::string prev_protein;
          std::string prev_chain;
          bool parse = true;
          do {
            if (prev_protein != index.protein() || prev_chain != index.chain()) {
              prev_protein = index.protein();
              prev_chain = index.chain();
              parse = selected.find(index.protein()) != selected.end() && proteins[index.protein()] == iterator->getBasicChainName(index.chain());
            }
            if (parse && neighbours.next(index.index())) {
              std::set<std::string> chains;
              std::stringstream stream(neighbours.value());
              std::string value;
              while (std::getline(stream, value, ' ')) {
                chains.emplace(value);
              }
              std::map<std::string, double> &counts = pairs[index.protein()][index.chain()];
              for (auto chains_it = chains.begin(); chains_it != chains.end(); ++chains_it) {
                auto counts_it = counts.find(*chains_it);
                if (counts_it == counts.end()) {
                  counts[*chains_it] = 1.0/chains.size();
                } else {
                  counts[*chains_it] += 1.0/chains.size();
                }
              }
            }
          } while (index.next());
        }

        OUTPUT << header;
        for (auto selected_id = selected.begin(); selected_id != selected.end(); ++selected_id) {
          if (selected_id != selected.begin()) {
            OUTPUT << ", ";
          }
          double max_c = 0;
          std::string max_v;
          auto &chain = pairs[*selected_id];
          for (auto chain_it = chain.begin(); chain_it != chain.end(); ++chain_it) {
            for (auto chain_it_it = chain_it->second.begin(); chain_it_it != chain_it->second.end(); ++chain_it_it) {
              if (chain_it_it->second > max_c) {
                max_c = chain_it_it->second;
                max_v = iterator->getBasicChainName(chain_it_it->first);
              }
            }
          }
          OUTPUT << *selected_id << '.' << proteins[*selected_id] << ':' << max_v;
        }
        OUTPUT << std::endl;
      }

      public:
      Random(std::string output) : OUTPUT(output) { }

      void select(size_t count, double threshold, std::string header, ProteinIterator* iterator, std::string index, std::string neighbours, std::string similarity) {
        gselect(count, threshold, header + ":\t", iterator, index, neighbours, similarity);
      }

      void select(size_t count, double threshold, ProteinIterator* iterator, std::string index, std::string neighbours, std::string similarity){
        gselect(count, threshold, "", iterator, index, neighbours, similarity);
      }
    };
  }
}