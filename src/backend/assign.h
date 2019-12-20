#pragma once

#include "index.h"

namespace inspire {
  namespace backend {
    class Assignator {
      private:
      // [Protein; model; chain; aminoacid]
      std::vector<std::tuple<std::string, std::string, std::string, std::string>> RESIDUES;
      // *.pec file
      std::ifstream INPUT;

      inline static std::string validate_name(std::string output, const std::string extension, const std::string input) {
        if (output.size() == 0 || output.back() == common::filesystem::directory_separator) {
          size_t sep = input.rfind(common::filesystem::directory_separator);
          if (sep == input.npos) {
            output += input;
          } else {
            output += input.substr(sep+1);
          }
        }
        if (!common::string::ends_with(output, extension)) {
          if (common::string::ends_with(output, ".pec")) {
            output = output.substr(0, output.size()-4);
          }
          output += extension;
        }
        return output;
      }

      public:
      Assignator(const std::string &index, const std::string &results) : INPUT(results) {
        Index indices(index);
        if (!indices.reset()) {
          throw common::exception::TitledException("The index file is empty or it is not possible to read it.");
        }
        // Residues file starts from 0;
        RESIDUES.push_back(std::make_tuple("Protein", "Model", "Chain", "Aminoacid"));
        do {
          if (RESIDUES.size() != indices.index()) {
            throw common::exception::TitledException("Unexpected format of index file: expected continuous aritmetic sequence starting at 1 with step 1");
          }
          RESIDUES.push_back(std::make_tuple(indices.protein(), indices.model(), indices.chain(), indices.aminoacid()));
        } while (indices.next());
      }

      bool next(std::string &protein, std::string &model, std::string &chain, std::string &aminoacid, std::string &value) {
        std::string line;
        if (std::getline(INPUT, line)) {
          size_t tab = line.find('\t');
          if (tab == line.npos || tab == line.size()-1) {
            throw common::exception::TitledException("Unexpected line format in results file: " + line);
          }
          auto &labels = RESIDUES[std::stoll(line.substr(0, tab))];
          protein = std::get<0>(labels);
          model = std::get<1>(labels);
          chain = std::get<2>(labels);
          aminoacid = std::get<3>(labels);
          value = line.substr(tab+1);
          return true;
        }
        return false;
      }

      // Delimiter-separated value file format
      static std::string Csv(char delimiter, const std::string &index, const std::string &results, const std::string &output) {
        if (delimiter == '"') {
          throw common::exception::TitledException("Double quotes are not allowed as a delimiter");
        }
        std::string output_path = validate_name(output, ".csv", results);
        std::ofstream stream(output_path);
        Assignator assignator(index, results);
        std::string protein;
        std::string model;
        std::string chain;
        std::string aminoacid;
        std::string value;
        stream << "\"Protein\"" << delimiter << "\"Model\"" << delimiter << "\"Chain\"" << delimiter << "\"Residue\"" << delimiter << "\"value\"\n";
        while (assignator.next(protein, model, chain, aminoacid, value)) {
          protein = common::string::replace_all(protein, '"', "\"\"");
          model = common::string::replace_all(model, '"', "\"\"");
          chain = common::string::replace_all(chain, '"', "\"\"");
          aminoacid = common::string::replace_all(aminoacid, '"', "\"\"");
          value = common::string::replace_all(value, '"', "\"\"");
          stream << '"' << protein << '"' << delimiter << '"' << model << '"' << delimiter << '"' << chain << '"' << delimiter << '"' << aminoacid << '"' << delimiter << '"' << value << "\"\n";
        }
        stream.flush();
        stream.close();
        return output_path;
      }

      // Aligned file format
      static std::string List(const std::string &index, const std::string &results, const std::string &output) {
        std::string output_path = validate_name(output, ".les", results);
        std::ofstream stream(output_path);
        Assignator assignator(index, results);
        std::string protein;
        std::string last_protein;
        std::string model;
        std::string last_model;
        std::string chain;
        std::string last_chain;
        std::string aminoacid;
        std::string value;
        bool first = true;
        while (assignator.next(protein, model, chain, aminoacid, value)) {
          if (first) {
            first = false;
            stream << protein << '\t' << model << '\t' << chain;
          } else {
            if (protein != last_protein) {
              stream << protein << '\t' << model << '\t' << chain;
            } else if (model != last_model) {
              stream << '\t' << model << '\t' << chain;
            } else if (chain != last_chain) {
              stream << "\t\t" << chain;
            } else {
              stream << "\t\t";
            }
          }
          stream << '\t' << aminoacid << '\t' << value << '\n';
          last_protein = protein;
          last_model = model;
          last_chain = chain;
        }
        stream.flush();
        stream.close();
        return output_path;
      }

      // Xml file format
      static std::string Xml(const std::string &index, const std::string &results, const std::string &output) {
        std::string output_path = validate_name(output, ".xml", results);
        std::ofstream stream(output_path);
        Assignator assignator(index, results);
        std::string protein;
        std::string last_protein;
        std::string model;
        std::string last_model;
        std::string chain;
        std::string last_chain;
        std::string aminoacid;
        std::string value;
        bool first = true;
        while (assignator.next(protein, model, chain, aminoacid, value)) {
          value = common::string::replace_all(common::string::replace_all(common::string::replace_all(value, '&', "&amp;"), '<', "&lt;"), '>', "&gt;");
          if (first) {
            first = false;
            stream << "<protein id=\"" << protein << "\">\n"
                   << "  <model id=\"" << model << "\">\n"
                   << "    <chain id=\"" << chain << "\">\n";
          } else {
            if (protein != last_protein) {
              stream << "    </chain>\n"
                     << "  </model>\n"
                     << "</protein>\n"
                     << "<protein id=\"" << protein << "\">\n"
                     << "  <model id=\"" << model << "\">\n"
                     << "    <chain id=\"" << chain << "\">\n";
            } else if (model != last_model) {
              stream << "    </chain>\n"
                     << "  </model>\n"
                     << "  <model id=\"" << model << "\">\n"
                     << "    <chain id=\"" << chain << "\">\n";
            } else if (chain != last_chain) {
              stream << "    </chain>\n"
                     << "    <chain id=\"" << chain << "\">\n";
            }
          }
          stream << "      <aminoacid id=\"" << aminoacid << "\">" << value << "</aminoacid>\n";
          last_protein = protein;
          last_model = model;
          last_chain = chain;
        }
        if (!first) {
          stream << "    </chain>\n"
                 << "  </model>\n"
                 << "</protein>\n";
        }
        stream.flush();
        stream.close();
        return output_path;
      }
    };
  }
}