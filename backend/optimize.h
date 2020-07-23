#pragma once

#include "../common/filesystem.h"
#include "features.h"
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

namespace inspire {
  namespace backend {
    class Optimizer {
      public:
      virtual void optimize(std::string input, std::string output) = 0;

      void optimize(std::string input) {
        if (common::string::ends_with(input, ".sas")) {
          optimize(input, input.substr(0, input.size()-4));
        } else if (common::string::ends_with(input, ".sed")) {
          optimize(input, input+".sed");
        } else {
          optimize(input, input);
        }
      }
    };

    class NaiveOptimizer : public Optimizer {
      protected:
      size_t PRECISSION;
      std::map<size_t, std::string> LABELS;

      public:
      NaiveOptimizer(std::string labels, size_t threshold) : PRECISSION(threshold) {
        FeaturesReader reader(labels);
        std::vector<int> features;
        for (auto i = 0; i < reader.size(); i++) {
          features.push_back(i);
        }
        while (reader.next_line()) {
          std::stringstream buffer;
          for (size_t i = 0; i < features.size(); i++) {
            if (i > 0) {
              buffer << common::filesystem::directory_separator;
            }
            buffer << reader.value(i);
          }
          const std::string &label = buffer.str();
          auto ins = LABELS.insert({reader.index(), label});
          if (!ins.second) {
            throw common::exception::TitledException("Multiple features with the same index " + reader.index());
          }
        }
      }
    };

    class FractionalOptimizer : public NaiveOptimizer {
      public:
      FractionalOptimizer(std::string labels, size_t threshold) : NaiveOptimizer(labels, threshold) { }

      void optimize(std::string input, std::string output) override {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          size_t i = input.rfind(common::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (common::string::ends_with(tmp, ".sas")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".sed";
        } else if (!common::string::ends_with(output, ".sed")) {
          output += ".sed";
        }

        std::ifstream reader(input);
        
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
        if (headers.size() > 2) {
          throw common::exception::TitledException("General n-dimensional optimisation is not supported yet.");
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
        writer << "threshold\tTP\tFN\tFP\tTN\tMCC\n";
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

    class WeightedOptimizer : public NaiveOptimizer {
      public:
      WeightedOptimizer(std::string labels, size_t threshold) : NaiveOptimizer(labels, threshold) { }

      void optimize(std::string input, std::string output) override {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          size_t i = input.rfind(common::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (common::string::ends_with(tmp, ".sas")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
          output += ".sed";
        } else if (!common::string::ends_with(output, ".sed")) {
          output += ".sed";
        }

        std::ifstream reader(input);

        std::string line;
        if (!std::getline(reader, line)) {
          throw common::exception::TitledException("Input file '" + input + "' is empty.");
        }

        std::map<std::string, std::pair<size_t, size_t>> indices;
        {
          std::stringstream parts(line);
          while (std::getline(parts, line, ':')) {
            indices[line] = std::pair<size_t, size_t>(indices.size(), 1);
          }
        }

        if (indices.size() > PRECISSION) {
          std::cerr << "The number of indices (" << indices.size() << ") id higher than the selected precission (" << PRECISSION << ")." << std::endl;
          return;
        }
        if (indices.size() == 0) {
          std::cerr << "There are no indices in the stats file." << std::endl;
          return;
        }
        if (indices.size() == 1) {
          std::cerr << "There is only one index in the stats file." << std::endl;
          return;
        }

        std::map<std::vector<size_t>, std::map<std::string, size_t>> stats;

        std::string id;
        while (std::getline(reader, id)) {
          std::vector<size_t> counts(indices.size(), 0);
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
          for (auto it = counts.begin(); it != counts.end(); ++it) {
            if (*it > 0) {
              ++stats[counts][LABELS[std::stol(id)]];
              break;
            }
          }
        }
        reader.close();

        // TODO: Consider empty stat - probably remove them

        std::ofstream writer(output);
        for (auto it = indices.begin(); it != indices.end(); ++it) {
          writer << it->first << '\t';
        }
        for (auto it1 = indices.begin(); it1 != indices.end(); ++it1) {
          for (auto it2 = indices.begin(); it2 != indices.end(); ++it2) {
            writer << it1->first << common::filesystem::directory_separator << it2->first << '\t';
          }
        }
        writer << "MCC\n";

        indices.begin()->second.second = PRECISSION - indices.size() + 1;
        double best_mcc = -1;
        std::vector<double> precissions;
        for (size_t i = 0; i <= indices.begin()->second.second; i++) {
          precissions.push_back(i/(double)PRECISSION);
        }
        do {
          // Calculate confusion matrix
          std::vector<std::vector<size_t>> results(indices.size(), std::vector<size_t>(indices.size(), 0));
          for (auto stats_it = stats.begin(); stats_it != stats.end(); ++stats_it) {
            double selected_score = 0;
            std::string selected_label;
            auto indices_it = indices.begin();
            for (size_t i = 0; indices_it != indices.end(); ++i, ++indices_it) {
              double score = stats_it->first[i]*indices_it->second.second;
              if (score > selected_score) {
                selected_score = score;
                selected_label = indices_it->first;
              }
            }
            for (auto stats_it_it = stats_it->second.begin(); stats_it_it != stats_it->second.end(); ++stats_it_it) {
              results[indices[stats_it_it->first].first][indices[selected_label].first] += stats_it_it->second;
            }
          }

          double a = 0;
          for (size_t k = 0; k < indices.size(); k++) {
            for (size_t l = 0; l < indices.size(); l++) {
              for (size_t m = 0; m < indices.size(); m++) {
                a += results[k][k]*results[l][m];
                a -= results[k][l]*results[m][k];
              }
            }
          }
          double b = 0;
          for (size_t k = 0; k < indices.size(); k++) {
            size_t b1 = 0;
            for (size_t l = 0; l < indices.size(); l++) {
              b1 += results[k][l];
            }
            size_t b2 = 0;
            for (size_t k2 = 0; k2 < indices.size(); k2++) {
              if (k != k2) {
                for (size_t l2 = 0; l2 < indices.size(); l2++) {
                  b2 += results[k2][l2];
                }
              }
            }
            b += b1*(double)b2;
          }
          double c = 0;
          for (size_t k = 0; k < indices.size(); k++) {
            size_t c1 = 0;
            for (size_t l = 0; l < indices.size(); l++) {
              c1 += results[l][k];
            }
            size_t c2 = 0;
            for (size_t k2 = 0; k2 < indices.size(); k2++) {
              if (k != k2) {
                for (size_t l2 = 0; l2 < indices.size(); l2++) {
                  c2 += results[l2][k2];
                }
              }
            }
            c += c1*(double)c2;
          }
          double mcc = a / std::sqrt(b*c);
          for (auto it = indices.begin(); it != indices.end(); ++it) {
            writer << precissions[it->second.second] << '\t';
          }
          for (auto it1 = indices.begin(); it1 != indices.end(); ++it1) {
            for (auto it2 = indices.begin(); it2 != indices.end(); ++it2) {
              writer << results[it1->second.first][it2->second.first] << '\t';
            }
          }
          writer << mcc << '\n';
          if (mcc > best_mcc) {
            best_mcc = mcc;
          }

          // Increment indices
          if (indices.rbegin()->second.second == PRECISSION - indices.size() + 1) {
            break;
          } else {
            auto it = indices.begin();
            while (it != indices.end() && it->second.second == 1) {
              ++it;
            }
            if (it == indices.begin()) {
              --(it->second.second);
            } else {
              indices.begin()->second.second = it->second.second-1;
              it->second.second = 1;
            }
            ++it;
            ++(it->second.second);
          }
        } while (true);
        writer.flush();
        std::cout << '\t' << best_mcc << std::endl;

        writer.flush();
        writer.close();
      }

    };
  }
}