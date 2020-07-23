#pragma once

#include "features.h"
#include <map>

namespace inspire {
  namespace backend {
    class Prioritizer {
      public:
      virtual void prioritize(std::string input, std::string output) = 0;
    };

    class LinearPrioritizer : public Prioritizer {
      private:
      double relative;
      double absolute;
      double distance;
      double coefficient;

      public:
      LinearPrioritizer(double relative, double absolute, double distance) {
        this->relative = relative;
        this->absolute = absolute;
        this->distance = distance;
        this->coefficient = 1/(absolute+relative+distance);
      }

      void prioritize(std::string input, std::string output) override {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          size_t i = input.rfind(common::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (common::string::ends_with(tmp, ".sas")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".pot";
        } else if (!common::string::ends_with(output, ".pot")) {
          output += ".pot";
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

        std::string id;
        while (std::getline(reader, id)) {
          double distance = 0;
          std::vector<size_t> counts(headers.size(), 0);
          while (std::getline(reader, line) && line.size() > 0) {
            std::stringstream parts(line);
            std::string level;
            if (!std::getline(parts, level, '\t')) {
              throw common::exception::TitledException("Unexpected exception during parsing line '" + line + "' for residue #" + id + ".");
            }
            std::string count;
            long long sum = 0;
            for (size_t i = 0; i < counts.size(); i++) {
              if (!std::getline(parts, count, ':')) {
                throw common::exception::TitledException("Incomplete line '" + line + "' for residue #" + id + ": missing numbers.");
              }
              long cnt = std::stol(count);
              sum += cnt;
              counts[i] += cnt;
            }
            // Consider unclassified as optional
            if (std::getline(parts, count, ':') && std::getline(parts, count, ':')) {
              throw common::exception::TitledException("Invalid line '" + line + "' for residue #" + id + ": too many numbers.");
            }
            distance += std::stol(level) * sum;
          }
          {
            long long count = 0;
            for (size_t i = 0; i < counts.size(); i++) {
              count += counts[i];
            }
            double score = count == 0 ? 0 : this->coefficient*(this->relative*counts[0]/(double)count + this->absolute*(1-1.0/(1+counts[0])) + this->distance/(1+distance/(double)count));
            writer << id << '\t' << score << '\n';
          }
        }
        reader.close();
        writer.flush();
        writer.close();
      }

    };
  }
}