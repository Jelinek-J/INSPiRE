#pragma once

#include "index.h"
#include "../common/exception.h"
#include "../common/random.h"
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <fstream>

namespace inspire {
  namespace backend {
    class Random {
      private:
      std::ofstream OUTPUT;

      void gselect(size_t count, double limit, std::string header, std::string index_path, std::string similarity_path) {
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
            if (std::stod(line.substr(0, first)) >= limit) {
              id.emplace(line.substr(first+1, second-first-1));
            }
          }
        }

        std::vector<std::string> order;
        std::map<std::string, char> proteins;
        inspire::backend::Index index(index_path);
        if (index.reset()) {
          std::string prev_protein;
          std::string prev_chain;
          std::vector<char> chains;
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
              if (new_chain.size() == 1 && prev_chain != new_chain) {
                chains.push_back(new_chain[0]);
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
        std::cout << header;
        for (auto selected_id = selected.begin(); selected_id != selected.end(); ++selected_id) {
          if (selected_id != selected.begin()) {
            std::cout << ", ";
          }
          std::cout << *selected_id << '.' << proteins[*selected_id];
        }
        std::cout << std::endl;
      }

      public:
      Random(std::string output) : OUTPUT(output) { }

      void select(size_t count, double limit, std::string header, std::string index, std::string similarity) {
        gselect(count, limit, header + ":\t", index, similarity);
      }

      void select(size_t count, double limit, std::string index, std::string similarity){
        gselect(count, limit, "", index, similarity);
      }
    };
  }
}