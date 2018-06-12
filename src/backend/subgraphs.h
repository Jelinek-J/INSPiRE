#pragma once

#include "features.h"
#include "index.h"

namespace inspire {
  namespace backend {
    // NOTE: Due to time and/or space demands, it is required that lines in index file are grouped by proteinId, modelId and chainId with proteinId has the highest priority and chainId the lowest priority
    //       and coordinates file is sorted by indices. Both is naturally met in the current implementation and it also results in the smallest sizes of both sizes.
    class Subgraphs{
      private:
      Index INDEX;
      FeatureReader FEATURES;

      std::string PATH;
      // k-nearest subgraphs that should be identified.
      std::map<int, std::ofstream> COUNT;
      // Radial distance subgraphs that should be identified.
      // NOTE: Distance is squared to avoid doing it during extraction (and also to avoid square root)
      std::map<double, std::ofstream> DISTANCE;
      // Path length subgraphs that should be identified; sorted by the length limit for a edge.
      // NOTE: Distance is squared to avoid doing it during extraction (and also to avoid square root)
      std::map<double, std::map<int, std::ofstream>> EDGES;

      void compute(const std::string &protein_id, const std::string &model_id, const std::vector<std::pair<std::string, std::map<size_t, Coordinate>>> &model) {
        int max_count;
        if (COUNT.rbegin() == COUNT.rend()) {
          max_count = 0;
        } else {
          max_count = COUNT.rbegin()->first;
        }
        double max_distance;
        if (DISTANCE.rbegin() == DISTANCE.rend()) {
          max_distance = 0;
        } else {
          max_distance = DISTANCE.rbegin()->first;
        }
        double max_length;
        if (EDGES.rbegin() == EDGES.rend()) {
          max_length = 0;
        } else {
          max_length = EDGES.rbegin()->first;
        }
        std::map<double, int> max_edges;
        for (auto &pair : EDGES) {
          if (pair.second.rbegin() == pair.second.rend()) {
            max_edges[pair.first] = 0;
          } else {
            max_edges[pair.first] = pair.second.rbegin()->first;
          }
        }

        for (auto chain_it = model.begin(); chain_it != model.end(); ++chain_it) {
          // TODO: Check time demands in comparision with list/vector/set of pairs and multimap as one final sort can be quicker than sort after every insert
          std::map<size_t, std::vector<std::pair<double, size_t>>> distances;
          for (auto residue_it = chain_it->second.begin(); residue_it != chain_it->second.end(); ++residue_it) {
            auto& working = distances[residue_it->first];
            // NOTE: This is not necessary, however, now will be the k-nearest element at the index k, thus it is much easier to think about it and check conditions
            working.push_back({0, residue_it->first});
            for (auto residue2_it = residue_it; ++residue2_it != chain_it->second.end(); ) {
              // NOTE: Squared value does not change the order of residues
              double distance = std::pow(std::get<0>(residue_it->second) - std::get<0>(residue2_it->second), 2)
                              + std::pow(std::get<1>(residue_it->second) - std::get<1>(residue2_it->second), 2)
                              + std::pow(std::get<2>(residue_it->second) - std::get<2>(residue2_it->second), 2);
              working.push_back({distance, residue2_it->first});
              distances[residue2_it->first].push_back({distance, residue_it->first});
            }
            std::sort(working.begin(), working.end());

            if (working[0].second == residue_it->first && (working.size() == 1 || working[1].first > 0)) {
              { // COUNT subgraphs
                for (auto count_it = COUNT.begin(); count_it != COUNT.end(); ++count_it) {
                  count_it->second << working[0].second;
                }
                // NOTE: If there is more elements in the same distance, they all will be taken.
                for (int i = 1; i < working.size() && (i <= max_count || working[i].first == working[max_count].first); ++i) {
                  // Streams that have a higher or equal count limit
                  for (auto count_it = COUNT.lower_bound(i); count_it != COUNT.end(); ++count_it) {
                    count_it->second << '\t' << working[i].second;
                  }
                  // Streams that have a lower limit, but the last residue within the limit is in the same distance as the current limit
                  for (auto it = COUNT.lower_bound(i); it != COUNT.begin() && (working[i].first == working[(--it)->first].first); ) {
                    it->second << '\t' << working[i].second;
                  }
                }
                for (auto count_it = COUNT.begin(); count_it != COUNT.end(); ++count_it) {
                  count_it->second << "\n";
                }
              }
              { // DISTANCE subgraphs
                for (auto distance_it = DISTANCE.begin(); distance_it != DISTANCE.end(); ++distance_it) {
                  distance_it->second << working[0].second;
                }
                for (size_t i = 1; i < working.size() && working[i].first < max_distance; ++i) {
                  // Streams that have a higher (or equal) distance limit
                  for (auto distance_it = DISTANCE.lower_bound(working[i].first); distance_it != DISTANCE.end(); ++distance_it) {
                    distance_it->second << '\t' << working[i].second;
                  }
                }
                for (auto distance_it = DISTANCE.begin(); distance_it != DISTANCE.end(); ++distance_it) {
                  distance_it->second << "\n";
                }
              }
            } else {
              std::cerr << protein_id << '-' << model_id << '-' << chain_it->first << ": Elements " << working[0].second << " and " << working[1].second
                << " have the same position (" << std::get<0>(residue_it->second) << "; " << std::get<1>(residue_it->second) << "; " << std::get<2>(residue_it->second) << ").\n";
            }
            // NOTE: Do not delete distant neighbours is not space-saving, but expected maximal size is less than 1GB thus time-saving is preferred.
          }
          for (auto edges_it = EDGES.rbegin(); edges_it != EDGES.rend(); ++edges_it) {
            for (auto residue_it = distances.begin(); residue_it != distances.end(); ++residue_it) {
              // NOTE: Consideration that vectors are sorted thanks previous stage
              // NOTE: In the reverse order, we can filter distances maps according to the length of edge
              // TODO: Is it really effective do delete them?
              while (!residue_it->second.empty() && residue_it->second.rbegin()->first > edges_it->first) {
                residue_it->second.pop_back();
              }
            }
            for (auto residue_it = distances.begin(); residue_it != distances.end(); ++residue_it) {
              for (auto edge_it = edges_it->second.begin(); edge_it != edges_it->second.end(); ++edge_it) {
                edge_it->second << residue_it->first;
              }
              // TODO: Effectivity of unordered_set compared to the (ordered) set  for (probably) small sets?
              std::set<size_t> neighbours;
              neighbours.insert(residue_it->first);
              for (int i = 1; i <= max_edges[edges_it->first]; i++) {
                std::set<size_t> next;
                for (auto neighbours_it = neighbours.begin(); neighbours_it != neighbours.end(); ++neighbours_it) {
                  auto add = distances.find(*neighbours_it);
                  if (add == distances.end()) {
                    std::cerr << "Unexpected exception during processing model '" << model_id << "' from protein '" << protein_id << "': residue no. " << *neighbours_it << " not found" << std::endl;
                  } else {
                    // TODO: Effectivity compared to the set<pair<double, int>> i.e. batch insert possible, but bigger objects are used (pair<double, int> instead of int)
                    for (auto add_it = add->second.begin(); add_it != add->second.end(); ++add_it) {
                      next.insert(add_it->second);
                    }
                  }
                }
                neighbours.insert(next.begin(), next.end());
                // NOTE: Expectatition that find and erase is cheaper that check every written residue; it must be removed to satisfy the condition the central residue first and only first
                neighbours.erase(residue_it->first);
                auto stream = edges_it->second.find(i);
                if (stream != edges_it->second.end()) {
                  for (auto neighbours_it = neighbours.begin(); neighbours_it != neighbours.end(); ++neighbours_it) {
                    stream->second << '\t' << *neighbours_it;
                  }
                }
              }
              for (auto edge_it = edges_it->second.begin(); edge_it != edges_it->second.end(); ++edge_it) {
                edge_it->second << "\n";
              }
            }
          }
        }

        for (auto count_it = COUNT.begin(); count_it != COUNT.end(); ++count_it) {
          count_it->second.flush();
        }
        for (auto distance_it = DISTANCE.begin(); distance_it != DISTANCE.end(); ++distance_it) {
          distance_it->second.flush();
        }
        for (auto edges_it = EDGES.begin(); edges_it != EDGES.end(); ++edges_it) {
          for (auto edge_it = edges_it->second.begin(); edge_it != edges_it->second.end(); ++edge_it) {
            edge_it->second.flush();
          }
        }
      }


      protected:
      // Classify groups within them complete subgraphs should lie.
      // I.E. <chain> in the case of subgraphs where all residues are in the same chain; and arbitrary const in the case that residues can be from whole protein
      virtual std::string chain_class(std::string chain) = 0;


      public:
      Subgraphs(const std::string index, const std::string coordinates, const std::string path) : INDEX(index), FEATURES(coordinates, CoordinateFeature(nullptr).title()), PATH(path) {
        if (!INDEX.reset()) {
          throw elemental::exception::TitledException("The index file is empty or it is not possible to read it.");
        }
      }
      ~Subgraphs() {
        for (auto count_it = COUNT.begin(); count_it != COUNT.end(); ++count_it) {
          count_it->second.close();
        }
        for (auto distance_it = DISTANCE.begin(); distance_it != DISTANCE.end(); ++distance_it) {
          distance_it->second.close();
        }
        for (auto edges_it = EDGES.begin(); edges_it != EDGES.end(); ++edges_it) {
          for (auto edge_it = edges_it->second.begin(); edge_it != edges_it->second.end(); ++edge_it) {
            edge_it->second.close();
          }
        }
      }

      void k_nearest(int k) {
        if (!COUNT.insert({k, std::ofstream(PATH + "c" + std::to_string(k) + ".sup")}).second) {
          std::cerr << "Duplicit request to generate subgraphs with " << k << "-nearest neighbours" << std::endl;
        }
      }

      void distance_limit(double d) {
        if (!DISTANCE.insert({std::pow(d, 2), std::ofstream(PATH + "d" + std::to_string(d) + ".sup")}).second) {
          std::cerr << "Duplicit request to generate subgraphs from neighbours within " << d << " distance" << std::endl;
        }
      }

      void edge_limit(double d, int e) {
        std::map<int, std::ofstream> &edges = EDGES[std::pow(d, 2)];
        if (!edges.insert({e, std::ofstream(PATH + "e" + std::to_string(d) + "_" + std::to_string(e) + ".sup")}).second) {
          std::cerr << "Duplicit request to generate subgraphs from neighbours within " << e << " edges distance where edge is between residues with distance at most " << d << std::endl;
        }
      }

      bool empty() {
        if (COUNT.size() > 0 || DISTANCE.size() > 0) {
          return false;
        }
        for (auto edges_it = EDGES.begin(); edges_it != EDGES.end(); ++edges_it) {
          if (edges_it->second.size() > 0) {
            return false;
          }
        }
        return true;
      }

      // TODO: Maybe it could be possible to reuse results from one chain to its transformations and thus save time if <fragmented> switcher is set
      void extract_subgraphs(ProteinIterator* iterator) {
        std::string protein_id = INDEX.protein();
        std::string model_id = INDEX.model();

        // Chain => { id => coordinates }
        std::vector<std::pair<std::string, std::map<size_t, Coordinate>>> model;
        model.push_back({chain_class(INDEX.chain()), std::map<size_t, Coordinate>()});
        do {
          if (protein_id != INDEX.protein()) {
            if (protein_id.size() > 0) {
              compute(protein_id, model_id, model);
              model.clear();
            }
            std::cout << protein_id << '\r';
            protein_id = INDEX.protein();
            model_id = INDEX.model();
            model.push_back({chain_class(INDEX.chain()), std::map<size_t, Coordinate>()});
          } else if (model_id != INDEX.model()) {
            compute(protein_id, model_id, model);
            model.clear();
            model_id = INDEX.model();
            model.push_back({chain_class(INDEX.chain()), std::map<size_t, Coordinate>()});
          } else {
            std::string ch_c = chain_class(INDEX.chain());
            if (model.back().first != ch_c) {
              model.push_back({ch_c, std::map<size_t, Coordinate>()});
            }
          }
          if (FEATURES.next(INDEX.index())) {
            model.back().second[INDEX.index()] = CoordinateFeature::parse(FEATURES.value()); 
          } else {
            std::cerr << "Residue no." << INDEX.index() << " is not presented in the coordinates file" << std::endl;
          }
        } while (INDEX.next());
        if (protein_id.size() > 0) {
          compute(protein_id, model_id, model);
        }
      }
    };

    class ProteinSubgraphs : public Subgraphs {
      protected:
      std::string chain_class(std::string chain) override { return ""; }
      public:
      ProteinSubgraphs(const std::string index, const std::string coordinates, const std::string path) : Subgraphs(index, coordinates, path) { }
    };

    class ChainSubgraphs : public Subgraphs {
      protected:
      std::string chain_class(std::string chain) override { return chain; }
      public:
      ChainSubgraphs(const std::string index, const std::string coordinates, const std::string path) : Subgraphs(index, coordinates, path) { }
    };

    class SubgraphReader {
      private:
      std::string NAME;
      std::ifstream INPUT;

      size_t INDEX;
      bool USED;
      // NOTE: NODES are withou central residue
      std::vector<int> NODES;

      protected:
      virtual void fill() {
        NODES.push_back(INDEX);
      }

      public:
      SubgraphReader(std::string file) : NAME(file), INPUT(file), INDEX(0), USED(true) { }

      bool next(size_t index) {
        if (index >= INDEX && (USED || index > INDEX)) {
          USED = false;
          do {
            std::string line;
            if (INPUT.eof() || !std::getline(INPUT, line)) {
              INDEX = std::numeric_limits<int>::max();
              return false;
            }
            std::stringstream parts(line);
            std::string part;
            if (!std::getline(parts, part, '\t')) {
              throw elemental::exception::TitledException("Unexpected empty row in subgraphs file '" + NAME + "'");
            }
            INDEX = std::stoi(part);
            NODES.clear();
            fill();
            while (std::getline(parts, part, '\t')) {
              NODES.push_back(std::stoi(part));
            }
          } while (index > INDEX);
        }
        USED = index == INDEX;
        return USED;
      }

      // NOTE: NODES are without central residue
      std::vector<int> value() {
        return NODES;
      }
    };
    class EdgesReader : public SubgraphReader {
      protected:
      void fill() override { }

      public:
      EdgesReader(std::string file) : SubgraphReader(file) { }
    };
  }
}