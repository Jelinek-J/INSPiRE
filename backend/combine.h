#pragma once

#include "../common/filesystem.h"
#include "../common/exception.h"
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

namespace inspire {
  namespace backend {
    class MutuallyOptimizer {
      private:
      static const size_t HEADER = 5;
      size_t DIMENSION = -1;
      std::vector<std::ifstream*> INPUTS;
      std::ofstream OUTPUT;

      protected:
      virtual double combine_row(std::vector<double> &values) = 0;

      public:
      MutuallyOptimizer(std::string output) : OUTPUT(output) { }

      ~MutuallyOptimizer() {
        OUTPUT.flush();
        OUTPUT.close();
        for (size_t i = 0; i < INPUTS.size(); i++) {
          INPUTS[i]->close();
          delete INPUTS[i];
        }
      }

      void add_input(std::string header, std::string input) {
        size_t index = INPUTS.size();
        INPUTS.push_back(new std::ifstream(input));
        std::string line;
        if (!std::getline(*INPUTS[index], line)) {
          throw common::exception::TitledException("File '" + input + "' is empty or it is not possible to read from it.");
        }
        std::stringstream parts(line);
        std::string part;
        size_t columns = 0;
        while (std::getline(parts, part, '\t')) {
          ++columns;
        }
        if (DIMENSION == -1) {
          if (columns < HEADER) {
            throw common::exception::TitledException("Unexpected format of optimization file '" + input + "': header has less than 5 elements.");
          }
          DIMENSION = columns - HEADER;
          OUTPUT << "threshold";
          for (size_t i = 1; i < DIMENSION; i++) {
            OUTPUT << '\t';
          }
        } else {
          if (DIMENSION + HEADER != columns) {
            throw common::exception::TitledException("File '" + input + "' has a different dimension then previous file(s).");
          }
        }
        OUTPUT << '\t' << header;
      }

      void combine() {
        double max = -2;
        if (INPUTS.size() == 0) {
          throw common::exception::TitledException("No file was added");
        }
        OUTPUT << '\n';
        std::string line;
        while (std::getline(*INPUTS[0], line)) {
          std::vector<double> values;
          // Well, tests causes slow down, but this is not time nor memory critical point, so it is not a problem
          std::vector<std::string> thresholds;
          {
            std::stringstream parts(line);
            std::string part;
            for (size_t i = 0; i < DIMENSION; i++) {
              if (!std::getline(parts, part, '\t')) {
                throw common::exception::TitledException("The first file contains a line with an insufficient count of columns");
              }
              thresholds.push_back(part);
              OUTPUT << part << '\t';
            }
            for (size_t i = 0; i < HEADER; i++) {
              if (!std::getline(parts, part, '\t')) {
                throw common::exception::TitledException("The first file contains a line with an insufficient count of columns");
              }
            }
            OUTPUT << part << '\t';
            values.push_back(std::stod(part));
          }
          for (size_t input = 1; input < INPUTS.size(); input++) {
            if (!std::getline(*INPUTS[input], line)) {
              throw common::exception::TitledException("The " + std::to_string(input) + "th file contains a line with an insufficient count of colomns");
            }
            std::stringstream parts(line);
            std::string part;
            for (size_t i = 0; i < DIMENSION; i++) {
              if (!std::getline(parts, part, '\t')) {
                throw common::exception::TitledException("The " + std::to_string(input) + "th file contains a line with an insufficient count of colomns");
              }
              if (part != thresholds[i]) {
                throw common::exception::TitledException("The " + std::to_string(input) + "th file and " + std::to_string(input) + "th file don't have the same thresholds");
              }
            }
            for (size_t i = 0; i < HEADER; i++) {
              if (!std::getline(parts, part, '\t')) {
                throw common::exception::TitledException("The first file contains a line with an insufficient count of colomns");
              }
            }
            OUTPUT << part << '\t';
            values.push_back(std::stod(part));
          }
          double value = combine_row(values);
          OUTPUT << value << '\n';
          if (value > max) {
            max = value;
          }
        }
        OUTPUT.flush();
        std::cout << max << std::endl;
      }
    };

    class AverageMutuallyOptimize : public MutuallyOptimizer {
      protected:
      double combine_row(std::vector<double> &values) override {
        double sum = 0;
        for (size_t i = 0; i < values.size(); i++) {
          sum += values[i];
        }
        return sum / values.size();
      }

      public:
      AverageMutuallyOptimize(std::string output) : MutuallyOptimizer(output) { }
    };

    class MedianMutuallyOptimize : public MutuallyOptimizer {
      protected:
      double combine_row(std::vector<double> &values) override {
        std::sort(values.begin(), values.end());
        return values.size() %2 == 0 ? (values[values.size()/2-1] + values[values.size()/2])/2 : values[values.size()/2];
      }

      public:
      MedianMutuallyOptimize(std::string output) : MutuallyOptimizer(output) { }
    };
  }
}