#pragma once

#include "../elemental/filesystem.h"
#include "features.h"
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

namespace inspire {
  namespace backend {
    class Optimizer {
      public:
      virtual void optimize(std::string input) = 0;
      virtual void optimize(std::string input, std::string output) = 0;
    };

    class FractionalOptimizer : public Optimizer {
      private:
      size_t PRECISSION;
      std::map<size_t, std::string> LABELS;

      public:
      FractionalOptimizer(std::string labels, size_t threshold) {
        FeaturesReader reader(labels);
        std::vector<int> features;
        for (auto i = 0; i < reader.size(); i++) {
          features.push_back(i);
        }
        while (reader.next_line()) {
          std::stringstream buffer;
          for (size_t i = 0; i < features.size(); i++) {
            if (i > 0) {
              buffer << elemental::filesystem::directory_separator;
            }
            buffer << reader.value(i);
          }
          const std::string &label = buffer.str();
          auto ins = LABELS.insert({reader.index(), label});
          if (!ins.second) {
            throw elemental::exception::TitledException("Multiple features with the same index " + reader.index());
          }
        }
        PRECISSION = threshold;
      }

      void optimize(std::string input) override {
        if (elemental::string::ends_with(input, ".sas")) {
          optimize(input, input.substr(0, input.size()-4));
        } else if (elemental::string::ends_with(input, ".den")) {
          optimize(input, input+".den");
        } else {
          optimize(input, input);
        }
      }

      void optimize(std::string input, std::string output) override {
        if (output.empty() || output.back() == elemental::filesystem::directory_separator) {
          size_t i = input.rfind(elemental::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (elemental::string::ends_with(tmp, ".sas")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".den";
        } else if (!elemental::string::ends_with(output, ".den")) {
          output += ".den";
        }

        std::ifstream reader(input);
        
        std::string line;
        if (!std::getline(reader, line)) {
          throw elemental::exception::TitledException("Input file '" + input + "' is empty.");
        }

        std::vector<std::string> headers;
        {
          std::stringstream parts(line);
          while (std::getline(parts,line,':')) {
            headers.push_back(line);
          }
        }
        if (headers.size() > 2) {
          throw elemental::exception::TitledException("General n-dimensional optimisation is not supported yet.");
        }

        size_t dimension = 1;
        for (size_t i = 1; i < headers.size(); i++) {
          dimension *= PRECISSION;
        }
        std::vector<std::vector<long>> stats;
        for (size_t i = 0; i < dimension; i++) {
          stats.push_back({0,0,0,0});
        }
        std::string id;
        while (std::getline(reader, id)) {
          std::vector<size_t> counts(headers.size(), 0);
          while (std::getline(reader, line) && line.size() > 0) {
            std::stringstream parts(line);
            std::string count;
            if (!std::getline(parts, count, '\t')) {
              throw elemental::exception::TitledException("Unexpected exception during parsing line '" + line + "' for residue #" + id + ".");
            }
            for (size_t i = 0; i < counts.size(); i++) {
              if (!std::getline(parts, count, ':')) {
                throw elemental::exception::TitledException("Incomplete line '" + line + "' for residue #" + id + ": missing numbers.");
              }
              counts[i] += std::stol(count);
            }
            // Consider unclassified as optional
            if (std::getline(parts, count, ':') && std::getline(parts, count, ':')) {
              throw elemental::exception::TitledException("Invalid line '" + line + "' for residue #" + id + ": too many numbers.");
            }
          }

          int iface = counts[0];
          int nface = counts[1];
          int ratio = (PRECISSION*iface+0.5) / (iface + nface);
          if (LABELS[std::stol(id)] == headers[0]) {
            for (int i = 0; i < ratio; i++) {
              stats[i][0]++;
            }
            for (int i = ratio; i < dimension; i++) {
              stats[i][1]++;
            }
          } else {
            for (int i = 0; i < ratio; i++) {
              stats[i][2]++;
            }
            for (int i = ratio; i < dimension; i++) {
              stats[i][3]++;
            }
          }
        }
        reader.close();


        double best = -1;
        std::ofstream writer(output);
        writer << "threshold\tTP\tFP\tFN\tTN\tMCC\n";
        for (size_t i = 0; i < PRECISSION; i++) {
          writer << (0.5+i) / PRECISSION;
          for (size_t j = 0; j < 4; j++) {
            writer << '\t' << stats[i][j];
          }
          double mcc = stats[i][0] + stats[i][1];
          mcc *= stats[i][0] + stats[i][2];
          mcc *= stats[i][1] + stats[i][3];
          mcc *= stats[i][2] + stats[i][3];
          if (mcc > 0) {
            mcc = (stats[i][0] * stats[i][3] - stats[i][1] * stats[i][2]) / std::sqrt(mcc);
          }
          writer << '\t' << mcc << '\n';
          best = std::max(best, mcc);
        }
        writer.flush();
        std::cout << '\t' << best << std::endl;

        writer.flush();
        writer.close();
      }

    };
  }
}