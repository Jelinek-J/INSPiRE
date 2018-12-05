#pragma once

#include "index.h"
#include "features.h"
#include <map>
#include <unordered_map>

namespace inspire {
  namespace backend {
    class Excluder {
      private:
      std::map<std::string, std::pair<std::vector<size_t>, std::vector<size_t>>> EXCLUDES;

      public:
      void add(std::string exclude) {
        EXCLUDES.insert({exclude, std::pair<std::vector<size_t>, std::vector<size_t>>()});
      }

      void generate(std::string knowledge_base, std::string queries, std::string output) {
        Index index(queries);
        if (index.reset()) {
          do {
            auto it = EXCLUDES.find(index.protein());
            if (it != EXCLUDES.end()) {
              it->second.first.push_back(index.index());
            }
          } while (index.next());
        }
        index = Index(knowledge_base);
        if (index.reset()) {
          do {
            auto it = EXCLUDES.find(index.protein());
            if (it != EXCLUDES.end()) {
              it->second.second.push_back(index.index());
            }
          } while (index.next());
        }

        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          output += "related.exc";
        }
        std::ofstream stream(output);
        for (auto excludes_it = EXCLUDES.begin(); excludes_it != EXCLUDES.end(); ++excludes_it) {
          if (excludes_it->second.first.size() > 0 && excludes_it->second.second.size() > 0) {
            for (auto q_it = excludes_it->second.first.begin(); q_it != excludes_it->second.first.end(); ++q_it) {
              if (q_it != excludes_it->second.first.begin()) {
                stream << ' ';
              }
              stream << *q_it;
            }
            stream << '\t';
            for (auto kb_it = excludes_it->second.second.begin(); kb_it != excludes_it->second.second.end(); ++kb_it) {
              if (kb_it != excludes_it->second.second.begin()) {
                stream << ' ';
              }
              stream << *kb_it;
            }
            stream << '\n';
          }
        }
        stream.flush();
        stream.close();
      }
    };
  }
}