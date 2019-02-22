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
        Index index_queries(queries);
        if (index_queries.reset()) {
          do {
            auto it = EXCLUDES.find(index_queries.protein());
            if (it != EXCLUDES.end()) {
              it->second.first.push_back(index_queries.index());
            }
          } while (index_queries.next());
        }
        Index index_kb(knowledge_base);
        if (index_kb.reset()) {
          do {
            auto it = EXCLUDES.find(index_kb.protein());
            if (it != EXCLUDES.end()) {
              it->second.second.push_back(index_kb.index());
            }
          } while (index_kb.next());
        }

        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          output += "related.exc";
        } else if (!common::string::ends_with(output, ".exc")) {
          output += ".exc";
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