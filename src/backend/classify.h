#pragma once

#include "features.h"
#include <map>

namespace inspire {
  namespace backend {
    class Classifier {
      private:
      std::set<std::string> headers;
      std::map<size_t, std::string> labels;
      // NOTE: Empty features skipped for memory reasons, however the right amount of separator is usefull for more efficient parsing
      std::string none;

      void load_labels(FeaturesReader &input, std::vector<int> features) {
        while (input.next_line()) {
          std::stringstream buffer;
          for (size_t i = 0; i < features.size(); i++) {
            if (i > 0) {
              buffer << elemental::filesystem::directory_separator;
              none.push_back(elemental::filesystem::directory_separator);
            }
            buffer << input.value(i);
          }
          const std::string &label = buffer.str();
          auto ins = labels.insert({input.index(), label});
          if (!ins.second) {
            throw elemental::exception::TitledException("Multiple features with the same index " + input.index());
          }
          headers.insert(label);
        }
      }

      public:
      Classifier(std::string input) {
        FeaturesReader reader(input);
        std::vector<int> features;
        for (auto i = 0; i < reader.size(); i++) {
          features.push_back(i);
        }
        load_labels(reader, features);
      }
      
      Classifier(std::string input, const std::vector<std::string> features) {
        FeaturesReader reader(input);
        std::map<std::string, size_t> map;
        for (size_t i = 0; i < reader.size(); i++) {
          map[reader.value(i)] = i;
        }
        std::vector<int> indices;
        for (auto features_it = features.begin(); features_it != features.end(); ++features_it) {
          auto it = map.find(*features_it);
          if (it == map.end()) {
            throw elemental::exception::TitledException("Feature '" + *features_it + "' is not present in feature file '" + input + "'");
          }
          indices.push_back(it->second);
        }
        load_labels(reader, indices);
      }

      void classify(std::string input, std::string output) {
        if (output.empty() || output.back() == elemental::filesystem::directory_separator) {
          size_t i = input.rfind(elemental::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (elemental::string::ends_with(tmp, ".med")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".sas";
        } else if (!elemental::string::ends_with(output, ".sas")) {
          output += ".sas";
        }
        std::ifstream input_file(input);
        std::ofstream output_file(output);
        for (auto headers_it = headers.begin(); headers_it != headers.end(); ++headers_it) {
          if (headers_it != headers.begin()) {
            output_file << ':';
          }
          output_file << *headers_it;
        }
        output_file << '\n';

        std::string line;
        while (!input_file.eof() && std::getline(input_file, line)) {
          output_file << line << '\n';
          std::map<int, std::map<std::string, int>> counts;
          while (!input_file.eof() && std::getline(input_file, line) && !line.empty()) {
            size_t index = line.find('\t');
            if (index == line.npos) {
              std::cerr << "Unexpected format of a line: '" << line << "'";
            } else {
              auto it = labels.find(std::stoi(line.substr(0, index)));
              counts[std::stoi(line.substr(index + 1))][it == labels.end() ? none : it->second]++;
            }
          }
          for (auto counts_it = counts.begin(); counts_it != counts.end(); ++counts_it) {
            output_file << counts_it->first << '\t';
            for (auto header_it = headers.begin(); header_it != headers.end(); ++header_it) {
              auto it = counts_it->second.find(*header_it);
              if (it == counts_it->second.end()) {
                output_file << 0;
              } else {
                output_file << it->second;
              }
              output_file << ':';
            }
            auto none_it = counts_it->second.find(none);
            if (none_it == counts_it->second.end()) {
              output_file << 0;
            } else {
              output_file << none;
            }
            output_file << '\n';
          }
          output_file << '\n';
        }
      }
    };
  }
}