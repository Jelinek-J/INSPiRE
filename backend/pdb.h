#pragma once

#include "protein.h"
#include "filters.h"
#include "../common/string.h"
#include "../common/exception.h"
#include <vector>
#include <set>
#include <istream>
#include <sstream>
#include <iostream>

// This file specializes on reading PDB files http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

namespace inspire {
  namespace backend {
    // A class grouping methods for parsing from pdb file
    struct Pdb {
      private:

      // Abstraction for parsing transformation matrices compatible with both REMARK 290 and REMARK350
      // lines: a list of lines to parse
      // shift: start index where to begin with parsing in the list
      // prefix: prefix for a format check, e.g. 'REMARK 290   SMTRY' or 'REMARK 350   BIOMT'
      // error_id: string for "localization of errors", e.g. 'crystallographic transformation' or 'biomolecule transformation'
      // returns: parsed transformation matrix
      private: static TransformationMatrix parseTransformation(const std::vector<std::string> &lines, size_t shift, std::string prefix, size_t id_end, std::string error_id) {
        const size_t end = shift+3;
        if (lines.size() < end) {
          throw common::exception::TitledException("Unexpected format of " + error_id + " matrix: missing line(s).");
        }
        for (size_t i = shift; i < end; i++) {
          if (lines[i].size() < 68 || !common::string::starts_with(lines[i], prefix) ||
              lines[i][18] != '1'+i-shift) {
            throw common::exception::TitledException("Unexpected format of " + error_id + " line: '"+ lines[i] + "'");
          }
        }
        size_t length = id_end-19;
        std::string id = lines[shift].substr(19, length);
        for (size_t i = shift+1; i < end; i++) {
          if (id != lines[i].substr(19, length)) {
            throw common::exception::TitledException("Unexpected format of " + error_id + " matrix: ambiguous identifiers: '"+ lines[i] + "'");
          }
        }
        std::vector<TransformationVector> transformations;
        for (size_t i = shift; i < end; i++) {
          transformations.push_back(std::make_tuple(std::stod(lines[i].substr(id_end, 10)), std::stod(lines[i].substr(id_end+10, 10)), std::stod(lines[i].substr(id_end+20, 10)), std::stod(lines[i].substr(id_end+30, 15))));
        }
        return std::make_tuple(transformations[0], transformations[1], transformations[2]);
      }

      // Finds end if transformation identifier for REMARK 290 and REMARK 350
      static size_t find_transformation_id_end(const std::string &line, const std::string &title) {
        size_t stop = line.find_first_not_of(' ', 19);
        if (stop == line.npos) {
          throw common::exception::TitledException("Unexpected format of " + title + " line, missing transformation identifier: '" + line + "'");
        }
        stop = line.find_first_of(' ', stop+1);
        if (stop == line.npos) {
          throw common::exception::TitledException("Unexpected format of " + title + " line, missing transformation coefficients: '" + line + "'");
        }
        return stop;
      }

      // Parse a single biomolecule from REMARK 350
      // lines: a list of lines to parse
      // protein: a protein where to insert parsed biomolecule
      // allowed_chains: a set of identifiers that are defined in a CHAIN line of COMPND section
      static void parseBiologicalTransformation(const std::vector<std::string> &lines, Protein &protein, std::map<char, std::pair<size_t, std::vector<std::string>>> &allowed_chains) {
        if (lines.size() == 0) {
          return;
        }
        std::stringstream comment;
        size_t i;
        for (i = 0; i < lines.size() && !common::string::starts_with(lines[i], "REMARK 350 APPLY THE FOLLOWING TO CHAINS:"); i++) {
          // TODO: Realy is test in each iteration quicker than a single one deletion of a single char?
          if (i > 0) {
            comment << '\n';
          }
          comment << lines[i];
        }
        auto biomolecule = protein.BIOMOLECULES.insert({std::stoi(lines[0].substr(23)), Biomolecule(comment.str())});
        if (!biomolecule.second) {
          throw common::exception::TitledException("Unexpected format of REMARK 350 section: multiple biomolecules with the same id: '" + lines[0] + "'");
        }
        std::set<std::string> chains;
        do {
          if (common::string::starts_with(lines[i], "REMARK 350 APPLY THE FOLLOWING TO CHAINS:")) {
            chains.clear();
            size_t j = 41;
            do {
              for (; j < lines[i].size(); ++j) {
                if (((lines[i][j] >= 'A' && lines[i][j] <= 'Z') || (lines[i][j] >= 'a' && lines[i][j] <= 'z') || (lines[i][j] >= '0' && lines[i][j] <= '9'))) {
                  if (allowed_chains.count(lines[i][j])) {
                    chains.emplace(std::string(1, lines[i][j]));
                  }
                } else if (lines[i][j] != ' ' && lines[i][j] != ',' &&lines[i][j] != ';') {

                }
              }
              ++i;
              if (common::string::starts_with(lines[i], "REMARK 350")) {
                j = 10;
                while (j < lines[i].size() && lines[i][j] == ' ') {
                  ++j;
                }
                if (common::string::contains_at(lines[i], "AND CHAINS:", j)) {
                  j += 12;
                } else {
                  j = -1;
                }
              } else {
                j = -1;
              }
            } while (i < lines.size() && j != -1);
          } else if (common::string::starts_with(lines[i], "REMARK 350   BIOMT1")) {
            size_t stop = find_transformation_id_end(lines[i], "biomolecule");
            // TODO: Does not it needlessly recreate objects?
            if (!biomolecule.first->second.TRANSFORMATIONS.insert({common::string::trim(lines[i].substr(19, stop-19)), {parseTransformation(lines, i, "REMARK 350   BIOMT", stop, "biomolecule"), chains}}).second) {
              throw common::exception::TitledException("Unexpected format of biomolecule: multiple matrices with the same id: " + lines[i] + "'");
            }
            i += 3;
          } else {
            for (size_t j = (common::string::starts_with(lines[i], "REMARK 350") ? 10 : 0); j < lines[i].size(); j++) {
              if (lines[i][j] != ' ') {
                throw common::exception::TitledException("Unexpected format of biomolecule: expected x-coordinate: " + lines[i] + "'");
              }
            }
            ++i;
          }
        } while (i < lines.size());
      }

      // Parses ATOM / HETATM line and adds corresponding Characteristic, Atom, Aminoacid, Chain and Model (if new) to the given protein
      // line: ATOM / HETATM line to parse
      // protein: protein where to add the result
      // model, chain, aminoacid: the currently processed model, chain and aminoacid respectivelly
      //
      // NOTE: The last three arguments serves to save time during searchin the corresponding element within maps
      // NOTO: There is not also an atom as alternative location are not common and thus it does not save time
      static void parse_atom_line(const std::string line, Protein& protein, std::pair<const int, Model>*& model, std::pair<const std::string, Chain>*& chain, std::pair<const std::pair<int, char>, Aminoacid>*& aminoacid) {
        if (model == nullptr) {
          if (protein.MODELS.empty()) {
            model = &*protein.MODELS.insert({1, Model()}).first;
          } else {
            throw common::exception::TitledException("'MODEL' line expected instead of '" + line + "'");
          }
          chain = nullptr;
          aminoacid = nullptr;
        }
        if (chain == nullptr) {
          auto ins = model->second.CHAINS.insert({std::string(1, line[21]), Chain()});
          if (!ins.second) {
            throw common::exception::TitledException("Multiple chains with the same identifier: '" + line + "'");
          }
          chain = &*ins.first;
          aminoacid = nullptr;
        } else if (chain->first[0] != line[21]) {
          throw common::exception::TitledException("Missing 'TER' line separating multiple chains: '" + line + "'");
        }
        const std::pair<int, char> aa_id = {std::stoi(line.substr(22, 4)), line[26]};
        if (aminoacid == nullptr || aa_id != aminoacid->first) {
          std::pair<std::map<std::pair<int, char>, Aminoacid>::iterator, bool> ins;
          if (common::string::starts_with(line, "ATOM  ")) {
            ins = chain->second.AMINOACIDS.insert({aa_id, Aminoacid(common::string::trim(line.substr(17,3)))});
          } else if (common::string::starts_with(line, "HETATM")) {
            ins = chain->second.AMINOACIDS.insert({aa_id, ModifiedAminoacid(common::string::trim(line.substr(17,3)))});
          } else {
            throw common::exception::TitledException("UNEXPECTED format of atomic coordinates line: '" + line + "'");
          }
          if (!ins.second) {
            throw common::exception::TitledException("Multiple aminoacids with the same identifier: '" + line + "'");
          }
          aminoacid = &*ins.first;
        }
        Atom& atom = aminoacid->second.ATOMS[common::string::trim(line.substr(12, 4))];
        if (line[16] == ' ') {
          if (atom.ALTERNATIVE_LOCATIONS.size() > 0) {
            throw common::exception::TitledException("Multiple ATOM lines for the same atom without altLoc specifier: '" + line + "'");
          }
          atom.ELEMENT = common::string::trim(line.substr(76, 2));
        } else {
          if (atom.ALTERNATIVE_LOCATIONS.empty()) {
            atom.ELEMENT = common::string::trim(line.substr(76, 2));
          } else if (atom.ELEMENT !=  common::string::trim(line.substr(76, 2))) {
            throw common::exception::TitledException("Multiple ATOM lines for the same atom contains differents elements: '" + line + "'");
          }
        }
        auto coordinate = atom.ALTERNATIVE_LOCATIONS.insert({line[16], Characteristic()});
        if (!coordinate.second) {
          throw common::exception::TitledException("Multiple alternative locations with the same identifier: '" + line + "'");
        }
        coordinate.first->second.SERIAL_NUMBER = std::stoi(line.substr(6, 5));
        coordinate.first->second.X = std::stod(line.substr(30, 8));
        coordinate.first->second.Y = std::stod(line.substr(38, 8));
        coordinate.first->second.Z = std::stod(line.substr(46, 8));
        coordinate.first->second.OCCUPANCY = std::stof(line.substr(54, 6));
        coordinate.first->second.TEMPERATURE = std::stof(line.substr(60, 6));
        if (line[78] == ' ' && line[79] == ' ') {
          coordinate.first->second.CHARGE = 0;
        } else {
          if (line[78] >= '0' && line[78] <= '9') {
            coordinate.first->second.CHARGE = line[78]-'0';
            if (line[79] == '-') {
              coordinate.first->second.CHARGE *= -1;
            } else if (line[79] != '+' && line[79] != ' ') {
              throw common::exception::TitledException("Unexpected format of charge: '" + line + "'");
            }
          } else if (line[79] >= '0' && line[79] <= '9') {
            coordinate.first->second.CHARGE = line[79]-'0';
            if (line[78] == '-') {
              coordinate.first->second.CHARGE *= -1;
            } else if (line[78] != '+' && line[78] != ' ') {
              throw common::exception::TitledException("Unexpected format of charge: '" + line + "'");
            }
          } else {
            throw common::exception::TitledException("Unexpected format of charge: '" + line + "'");
          }
        }
      }

      static inline std::string parse_id(std::string line) {
        if (line.size() < 66) {
          throw common::exception::TitledException("Unexpected format of HEADER line, identifier should be at positions [62; 66): '" + line + "'");
        }
        return line.substr(62, 4);
      }

      public:
      static std::string parse_id(std::istream& input) {
        std::string line;
        while (!input.eof() && std::getline(input, line) && !(common::string::starts_with(line, "END") && (line.size() == 3 || line[3] == ' '))) {
          if (common::string::starts_with(line, "HEADER")) {
            return parse_id(line);
          }
        }
        throw common::exception::TitledException("No header line");
      }

      // Read pdb file from input filtered by filters and returns the corresponding protein
      //
      // TODO: Currently, there are parsed only aminoacids with ATOM/HETATM line, in the future, it will be better to use SEQRES for it.
      // However, if an aminoacid has no coordinates, it is totally unusefull for us. I.E. it should be coordinated with the validation part of processing.
      static Protein parse_protein(std::istream& input, BasicFilter* filter) {
        Protein protein;
        std::pair<const int, Model>* model = nullptr;
        std::pair<const std::string, Chain>* chain = nullptr;
        std::pair<const std::pair<int, char>, Aminoacid>* aminoacid = nullptr;
        std::map<char, std::pair<size_t, std::vector<std::string>>> chains;
        bool parse_chains = false;
        std::string line;
        std::vector<std::string> crystallographic;
        std::vector<std::string> biomolecule;
        bool missing_model = false;

        // NOTE: input.eof() check is necessary for the case that the file does not ends with an empty line
        while (!input.eof() && std::getline(input, line) && !(common::string::starts_with(line, "END") && (line.size() == 3 || line[3] == ' '))) {
          if (filter->keep(line)) {
            if (common::string::starts_with(line, "HEADER")) {
              protein.ID_CODE = parse_id(line);
            } else if (common::string::starts_with(line, "COMPND")) {
              if (parse_chains || common::string::contains_at(line, "CHAIN: ", 11)) {
                size_t start = parse_chains ? 10 : 17;
                size_t i = start-1;
                size_t state = 0; // 0=space; 1=chain; 2=separator;3=end;
                while (state < 3 && ++i < line.size()) {
                  switch (line[i]) {
                    case ' ':
                      switch (state) {
                        case 0:
                          ++state;
                          break;
                        case 1:
                          while (i+1 < line.size() && line[i+1] == ' ') {
                            ++i;
                          }
                          if (i == line.size()) {
                            parse_chains = true;
                            state = 3;
                          } else {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN contains multiple consecutive spaces: '" << line << "'" << std::endl;
                          }
                          break;
                        case 2:
                          while (i+1 < line.size() && line[i+1] == ' ') {
                            ++i;
                          }
                          if (i == line.size()) {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line does not ends with a semicolon: '" << line << "'" << std::endl;
                            parse_chains = false;
                            state = 3;
                          } else {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN contains a space between a chain identifier and a separator: '" << line << "'" << std::endl;
                          }
                          break;
                        default:
                          throw common::exception::TitledException("Unexpected state occured within line: '" + line + "'");
                      }
                      break;
                    case ',':
                      switch (state) {
                        case 0:
                          std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN contains multiple consecutive separators: '" << line << "'" << std::endl;
                          break;
                        case 1:
                          std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN miss a chain identifier between two consecutive separators: '" << line << "'" << std::endl;
                          break;
                        case 2:
                          state=0;
                          break;
                        default:
                          throw common::exception::TitledException("Unexpected state occured within line: '" + line + "'");
                      }
                      break;
                    case ';':
                      switch (state) {
                        case 0:
                          if (i <= start) {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line contains no chain identifier: '" << line << "'" << std::endl;
                          } else {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line contains an extra comma in front of the semicolon: '" << line << "'" << std::endl;
                          }
                          break;
                        case 1:
                          if (i <= start+1) {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line contains no chain identifier: '" << line << "'" << std::endl;
                          } else {
                            std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN miss a chain identifier in front of the semicolon: '" << line << "'" << std::endl;
                          }
                          break;
                        case 2:
                          break;
                        default:
                          throw common::exception::TitledException("Unexpected state occured within line: '" + line + "'");
                      }
                      while (++i < line.size()) {
                        if (line[i] != ' ') {
                          std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line contains a non-space character after the semicolon: '" << line << "'" << std::endl;
                          break;
                        }
                      }
                      parse_chains = false;
                      state = 3;
                      break;
                    default:
                      switch (state) {
                        case 0:
                          std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line does not contain a space in front of a chain identifier: '" << line << "'" << std::endl;
                          break;
                        case 1:
                          break;
                        default:
                          throw common::exception::TitledException("Unexpected state occured within COMPOUND's CHAIN line: '" + line + "'");
                      }
                      chains[line[i]];
                      state = 2;
                      break;
                  }
                }
                switch (state) {
                  case 0:
                  case 1:
                    parse_chains = true;
                    break;
                  case 2:
                    parse_chains = false;
                    std::cerr << "WARNING [" << protein.ID_CODE << "]:    COMPOUND's CHAIN line does not ends with a semicolon: '" << line << "'" << std::endl;
                    break;
                  case 3:
                    break;
                  default:
                    throw common::exception::TitledException("Unexpected state occured within line: '" + line + "'");
                }
              }
            } else if (common::string::starts_with(line, "REMARK 290")) {
              if (common::string::starts_with(line, "REMARK 290   SMTRY")) {
                crystallographic.push_back(line);
                if (line.size() > 19 && line[18] == '3') {
                  switch (crystallographic.size()) {
                    case 1:
                    case 2:
                      throw common::exception::TitledException("Unexpected format of crystallographic transformation matrix: missing line(s) prior line: '" + line + "'");
                      break;
                    case 3:
                    {
                      size_t stop = find_transformation_id_end(crystallographic[0], "crystallographic transformation");
                      if (!protein.CRYSTALLOGRAPHIC_SYMMETRIES.insert({std::stoi(crystallographic[0].substr(19, stop-19)), parseTransformation(crystallographic, 0, "REMARK 290   SMTRY", stop, "crystallographic transformation")}).second) {
                        throw common::exception::TitledException("Unexpected format of crystallographic transformations: multiple matrices with the same id: '" + crystallographic[0] + "'");
                      }
                      crystallographic.clear();
                    }
                    break;
                    default:
                      throw common::exception::TitledException("Unexpected format of crystallographic transformation matrix: excess line(s): '" + crystallographic[0] + "'");
                      break;
                  }
                }
              }
            } else if (common::string::starts_with(line, "REMARK 350")) {
              if (common::string::starts_with(line, "REMARK 350 BIOMOLECULE:")) {
                parseBiologicalTransformation(biomolecule, protein, chains);
                biomolecule.clear();
                biomolecule.push_back(line);
              } else if (biomolecule.size() > 0) {
                biomolecule.push_back(line);
              }
            } else if (common::string::starts_with(line, "SEQRES ")) {
              size_t count = std::stol(common::string::trim(line.substr(13, 4)));
              auto chains_it = chains.find(line[11]);
              if (chains_it == chains.end()) {
                if (line[11] == ' ') {
                  switch (chains.size()) {
                    case 0:
                      throw common::exception::TitledException("Allowed chains are not definned prior SEQRES section: '" + line + "'");
                      break;
                    case 1:
                      chains_it = chains.begin();
                      break;
                    default:
                      throw common::exception::TitledException("Blank chain identifier, but multiple chains: '" + line + "'");
                  }
                } else {
                  throw common::exception::TitledException("Unknown chain identifier: '" + line + "'");
                }
              }
              if (std::stol(common::string::trim(line.substr(7, 3))) == 1) {
                if (chains_it->second.first != 0) {
                  throw common::exception::TitledException("Multiple lines with the same serial numbers corresponding to the same chain: '" + line + "'");
                }
                chains_it->second.first = count;
              } else {
                if (chains_it->second.first != count) {
                  throw common::exception::TitledException("Ambiguity in number of residues : '" + std::string(1, chains_it->second.first) + "' vs. '" + line + "'");
                }
              }
              for (size_t i = 19; i < line.size(); i+=4) {
                std::string residue = common::string::trim(line.substr(i, 3));
                if (residue.empty()) {
                  break;
                }
                chains_it->second.second.push_back(residue);
              }
            } else if (common::string::starts_with(line, "MODEL")) {
              int serial = std::stoi(line.substr(10, 4));
              auto ins = protein.MODELS.insert({serial, Model()});
              if (!ins.second) {
                throw common::exception::TitledException("Unexpected format of pdb file: multiple models with the same serial number: '" + line + "'");
              }
              if (!missing_model && serial != protein.MODELS.size()) {
                missing_model = true;
                std::cerr << "WARNING [" << protein.ID_CODE << "]:    Missing models with serial numbers from " << protein.MODELS.size() << " to " << serial-1 << std::endl;
              }
              model = &*ins.first;
            } else if (common::string::starts_with(line, "ENDMDL")) {
              model = nullptr;
              chain = nullptr;
              aminoacid = nullptr;
            } else if (common::string::starts_with(line, "TER")) {
              chain = nullptr;
              aminoacid = nullptr;
            } else if (common::string::starts_with(line, "ATOM")) {
              parse_atom_line(line, protein, model, chain, aminoacid);
            } else if (common::string::starts_with(line, "HETATM")) {
              if (line.size() < 22) {
                throw common::exception::TitledException("Unexpected format of HETATM line, too short: '" + line + "'");
              }
              // NOTE: if chain is not nullptr, than it is should be the same chain, else it is corrupted and throw an exception
              // If a model is nullptr, than it is corrupted and should throw an exception
              if (chain == nullptr && (chains.find(line[21]) == chains.end() || (model != nullptr && model->second.CHAINS.count(std::string(1, line[21]))))) {
                // TODO: If additional heteroatoms will be required too
              } else {
                parse_atom_line(line, protein, model, chain, aminoacid);
              }
            }
          }
        }
        if (parse_chains) {
          throw common::exception::TitledException("Unexpected format of COMPND section, CHAIN sections is not terminated.");
        }
        if (crystallographic.size() > 0) {
          throw common::exception::TitledException("Unexpected format of crystallographic transformation matrix: missing line with z-coordinate: '" + crystallographic[0] + "'");
        }
        parseBiologicalTransformation(biomolecule, protein, chains);

        if (protein.BIOMOLECULES.empty()) {
          std::cerr << "WARNING [" << protein.ID_CODE << "]:    Protein does not contain a definition of biomolecules (REMARK 350)" << std::endl;
          Biomolecule &biomolecule = protein.BIOMOLECULES[1];
          auto &transformation = biomolecule.TRANSFORMATIONS["1"];
          std::get<0>(std::get<0>(transformation.first)) = 1;
          std::get<1>(std::get<1>(transformation.first)) = 1;
          std::get<2>(std::get<2>(transformation.first)) = 1;
          for (auto it = chains.begin(); it != chains.end(); ++it) {
            transformation.second.emplace(std::string(1, it->first));
          }
        }

        if (protein.CRYSTALLOGRAPHIC_SYMMETRIES.empty()) {
          //std::cerr << "WARNING [" << protein.ID_CODE << "]:    Protein does not contain a definition of crystallographic transformation (REMARK 290)" << std::endl;
          TransformationMatrix &transformation = protein.CRYSTALLOGRAPHIC_SYMMETRIES[1];
          std::get<0>(std::get<0>(transformation)) = 1;
          std::get<1>(std::get<1>(transformation)) = 1;
          std::get<2>(std::get<2>(transformation)) = 1;
        }

        for (auto validation_it = chains.begin(); validation_it != chains.end(); ++validation_it) {
          if (validation_it->second.first != validation_it->second.second.size()) {
            throw common::exception::TitledException("Ambiguity between declared and real number of residues in the SEQRES record for chain '" + std::string(1, validation_it->first) + "'");
          }
          for (auto models_it = protein.MODELS.begin(); models_it != protein.MODELS.end(); ++models_it) {
            auto chains_it = models_it->second.CHAINS.find(std::string(1, validation_it->first));
            if (chains_it == models_it->second.CHAINS.end()) {
              std::cerr << "WARNING [" << protein.ID_CODE << "]:    Chain '" << validation_it->first << "' in model '" << models_it->first << "' does not contain information about " << validation_it->second.first << " residues" << std::endl;
            } else if (chains_it->second.AMINOACIDS.size() != validation_it->second.first) {
              if (chains_it->second.AMINOACIDS.size() < validation_it->second.first) {
                std::cerr << "WARNING [" << protein.ID_CODE << "]:    Chain '" << validation_it->first << "' in model '" << models_it->first << "' does not contain information about " << validation_it->second.first - chains_it->second.AMINOACIDS.size() << " residues" << std::endl;
              } else {
                std::cerr << "WARNING [" << protein.ID_CODE << "]:    Chain '" << validation_it->first << "' in model '" << models_it->first << "' has defined more residues in ATOM/HETATM section than in SEQRES section" << std::endl;
              }
            }
          }
        }

        return protein;
      }

      static bool contains_hetatoms(std::istream& input) {
        std::string line;
        while (std::getline(input, line)) {
          if (common::string::starts_with(line, "HETATM") && line.substr(17, 3) != "HOH" && line.substr(17, 3) != "DOD" && line.substr(17, 3) != "DIS" && line.substr(17, 3) != "MTO") {
            return true;
          }
        }
        return false;
      }
    };
  }
}