#pragma once

#include "features.h"
#include <map>

namespace inspire {
  namespace backend {
    class Predictor {
      public:
      virtual void predict(std::string input, std::string output) = 0;
    };

    class FractionalPredictor : public Predictor {
      private:
      std::vector<std::vector<double>> THRESHOLDS;

      public:
      FractionalPredictor(std::string thresholds) {
        std::vector<double> tmp;
        std::stringstream parts(thresholds);
        std::string part;
        while (std::getline(parts, part, ' ')) {
          double threshold = std::stof(part);
          if (threshold < 0 || threshold > 1) {
            throw common::exception::TitledException("Threshold '" + part + "' is out of boundaries; all threshold must be within interval [0; 1].");
          }
          tmp.push_back(threshold);
        }

        size_t dimension = std::sqrt(tmp.size()*2)+1.5;
        if (tmp.size() != dimension*(dimension-1)/2) {
          throw common::exception::TitledException("Invalid number of threshold " + std::to_string(tmp.size()) + " instead of most probable " + std::to_string(dimension*(dimension-1)/2));
        }

        size_t counter = 0;
        for (size_t i = 0; i < dimension; ++i) {
          std::vector<double> row;
          for (size_t j = 0; j < dimension-i-1; ++j) {
            row.push_back(tmp[counter++]);
          }
          THRESHOLDS.push_back(row);
        }
      }

      void predict(std::string input, std::string output) override {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          size_t i = input.rfind(common::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (common::string::ends_with(tmp, ".sas")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".pec";
        } else if (!common::string::ends_with(output, ".pec")) {
          output += ".pec";
        }

        std::ifstream reader(input);
        std::ofstream writer(output);
        std::string line;
        if (!std::getline(reader, line)) {
          throw common::exception::TitledException("Input file '" + input + "' is empty.");
        }

        std::vector<std::string> headers;
        {
          std::stringstream parts(line);
          while (std::getline(parts,line,':')) {
            headers.push_back(line);
          }
        }
        if (headers.size() != THRESHOLDS.size()) {
          throw common::exception::TitledException("Input file '" + input + "' is not compatible with the selector: invalid number of classes.");
        }

        std::string id;
        while (std::getline(reader, id)) {
          std::vector<size_t> counts(THRESHOLDS.size(), 0);
          while (std::getline(reader, line) && line.size() > 0) {
            std::stringstream parts(line);
            std::string count;
            if (!std::getline(parts, count, '\t')) {
              throw common::exception::TitledException("Unexpected exception during parsing line '" + line + "' for residue #" + id + ".");
            }
            for (size_t i = 0; i < counts.size(); i++) {
              if (!std::getline(parts, count, ':')) {
                throw common::exception::TitledException("Incomplete line '" + line + "' for residue #" + id + ": missing numbers.");
              }
              counts[i] += std::stol(count);
            }
            // Consider unclassified as optional
            if (std::getline(parts, count, ':') && std::getline(parts, count, ':')) {
              throw common::exception::TitledException("Invalid line '" + line + "' for residue #" + id + ": too many numbers.");
            }
          }
          size_t best = 0;
          for (size_t i = 1; i < THRESHOLDS.size(); i++) {
            if (counts[i] > 0 && counts[best]/(double)(counts[best]+counts[i]) < THRESHOLDS[best][i-best-1]) {
              best = i;
            }
          }
          writer << id << '\t' << (counts[best] == 0 ? "" : headers[best]) << '\n';
        }
        reader.close();
        writer.flush();
        writer.close();
      }

    };
  }
}