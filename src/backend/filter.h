#pragma once

#include "index.h"
#include "features.h"
#include <map>
#include <unordered_map>

namespace inspire {
  namespace backend {
    class Filter {
      private:
      std::vector<size_t> MAPPING;

      public:
      Filter(std::string original, std::string filtered) {
        std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, size_t>>>> index;
        Index orig(original);
        if (orig.reset()) {
          do {
            // TODO: Not optimal, but this should not be critical
            index[orig.protein()][orig.model()][orig.chain()][orig.aminoacid()] = orig.index()-1;
            MAPPING.push_back(0);
            if (MAPPING.size() != orig.index()) {
              throw common::exception::TitledException("Unexpected format of index file");
            }
          } while (orig.next());
        }
        Index filt(filtered);
        if (filt.reset()) {
          do {
            // TODO: Not optimal, but this should not be critical
            auto protein_it = index.find(filt.protein());
            if (protein_it != index.end()) {
              auto model_it = protein_it->second.find(filt.model());
              if (model_it != protein_it->second.end()) {
                auto chain_it = model_it->second.find(filt.chain());
                if (chain_it != model_it->second.end()) {
                  auto aminoacid_it = chain_it->second.find(filt.aminoacid());
                  if (aminoacid_it != chain_it->second.end()) {
                    MAPPING[aminoacid_it->second] = filt.index();
                  }
                }
              }
            }
          } while (filt.next());
        }
      }

      void filter(std::string original, std::string filtered) {
        FeaturesReader features(original);
        std::ofstream output(common::filesystem::copy_filename(filtered, original, ".tur"));
        for (size_t i = 0; i < features.size(); i++) {
          if (i > 0) {
            output << '\t';
          }
          output << features.header(i);
        }
        output << '\n';
        size_t previous = 0;
        while (features.next_line()) {
          size_t index = MAPPING[features.index()-1];
          if (index > 0) {
            for (size_t i = 0; i < features.size(); i++) {
              if (i > 0) {
                output << '\t';
              }
              output << features.value(i);
            }
            if (++previous != index) {
              previous = index;
              output << '\t' << index;
            }
            output << '\n';
          }
        }
        output.flush();
        output.close();
      }
    };
  }
}