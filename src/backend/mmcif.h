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
    struct Mmcif {
      private:

      // Parse value from the given row; to make it faster, the key is expected
      // NOTE: Is . the same as empty value ('', "", resp. ;\n;)?
      // NOTE: Suspection that no value is equal to ? ('?', "?", resp ;?\n;) and thus no ambiguation with ? as unknown value.
      static inline std::string parse_item(std::string &line, std::istream &input, const std::string &key) {
        size_t nonspace = line.find_first_not_of(' ', key.size());
        if (nonspace == line.npos) {
          if (!std::getline(input, line)) {
            throw common::exception::TitledException("Unexpected end of file, missing value for " + line);
          }
          if (line.size() < 1 || line[0] != ';') {
            return line;
          }
          std::stringstream buffer(line.substr(1));
          while (std::getline(input, line)) {
            if (line.size() > 0 || line[0] == ';') {
              return buffer.str();
            }
          }
          throw common::exception::TitledException("Unexpected end of file, uncompleted value for key " + key + ": " + buffer.str());
        } else {
          std::string value = common::string::trim(line.substr(nonspace));
          if (value[0] == '\'' || value[0] == '"') {
            if (value.size() < 2 || value[value.size()-1] != value[0]) {
              throw common::exception::TitledException("Unexpected format of line: " + line);
            }
            return value.substr(1, value.size()-2);
          } else if (value == "." || value == "?") {
            return "";
          } else {
            return value;
          }
        }
      }

      // Extract value from a vector and skip reading input of the last line of the row
      // NOTE: Is . the same as empty value ('', "", resp. ;\n;)?
      // NOTE: Suspection that no value is equal to ? ('?', "?", resp ;?\n;) and thus no ambiguation with ? as unknown value.
      static inline std::map<size_t, std::string> extract_item(std::string &line, std::istream &input, const size_t count, std::set<size_t> indices) {
        std::map<size_t, std::string> ret;
        size_t position = line.find_first_not_of(' ', 0);
        for (size_t i = 0; i < count; i++) {
          while (position == line.npos) {
            if (!std::getline(input, line)) {
              throw common::exception::TitledException("Unexpected end of file");
            }
            position = line.find_first_not_of(' ', 0);
          }
          switch (line[position]) {
            case '\'':
            case '"':
              {
                ++position;
                size_t next = line.find(line[position-1], position);
                if (next == line.npos) {
                  throw common::exception::TitledException("Missing closing quotation marks on line " + line);
                }
                if (indices.find(i) != indices.end()) {
                  ret[i] = line.substr(position, next-position);
                }
                if (next < line.size()-1 && line[next+1] != ' ') {
                  throw common::exception::TitledException("Unexpected character after closing quotation marks on line " + line);
                }
                position = next+1;
              }
              break;
            case ';':
              if (position == 0) {
                if (indices.find(i) != indices.end()) {
                  ret[i] = line.substr(1);
                }
                if (!std::getline(input, line)) {
                  throw common::exception::TitledException("Unexpected end of file within a matrix " + line);
                }
                while (line.size() == 0 || line[0] != ';') {
                  if (indices.find(i) != indices.end()) {
                    ret[i] += line;
                  }
                  if (!std::getline(input, line)) {
                    throw common::exception::TitledException("Unexpected end of file within a matrix " + line);
                  }
                }
                if (line.size() > 1 && line[1] != ' ') {
                  throw common::exception::TitledException("Unexpected character after closing semicolon of multiline value " + line);
                }
                position = 1;
              } else {
                size_t next = line.find(' ', position);
                if (indices.find(i) != indices.end()) {
                  if (next == line.npos) {
                    ret[i] = line.substr(position);
                  } else {
                    ret[i] = line.substr(position, next-position);
                  }
                }
                position = next;
              }
              break;
            default:
              size_t next = line.find(' ', position);
              if (indices.find(i) != indices.end()) {
                std::string tmp;
                if (next == line.npos) {
                  tmp = line.substr(position);
                } else {
                  tmp = line.substr(position, next-position);
                }
                ret[i] = (tmp == "." || tmp == "?") ? "" : tmp;
              }
              position = next;
              break;
          }
          position = line.find_first_not_of(' ', position);
        }
        return ret;
      }

      // Add biomolecule header to list of biomolecules
      static std::set<std::string> parseChains(std::string &affected) {
        std::set<std::string> ret;
        std::stringstream parts_chains(affected);
        std::string part;
        while (std::getline(parts_chains, part, ',')) {
          ret.emplace(part);
        }
        return ret;
      }

      // Parse what transformations should be used
      static std::map<std::string, TransformationMatrix> parseTransformations(const std::string &ids, const std::map<std::string, TransformationMatrix> &transformations) {
        if (ids.find('*') != ids.npos || ids.find('+') != ids.npos || ids.find('.') != ids.npos) {
          throw common::exception::TitledException("Star '*', plus '+' and dot '.' are not currently supported in identifiers of pdbx_struct_oper_list");
        }
        std::map<std::string, TransformationMatrix> ret;
        auto transformations_it = transformations.find(ids);
        if (transformations_it != transformations.end()) {
          ret[ids] = transformations_it->second;
          return ret;
        }
        std::stringstream items(ids);
        for (std::string item; std::getline(items, item, ','); ) {
          transformations_it = transformations.find(item);
          if (transformations_it == transformations.end()) {
            ret.clear();
            break;
          } else {
            ret[item] = transformations_it->second;
          }
        }
        if (ret.size() > 0) {
          return ret;
        }
        size_t minus = ids.find('-');
        if (minus != ids.npos) {
          long long from = std::stoll(ids.substr(0, minus));
          long long to = std::stoll(ids.substr(minus+1));
          for (long long i = from; i <= to; ++i) {
            transformations_it = transformations.find(std::to_string(i));
            if (transformations_it == transformations.end()) {
              ret.clear();
              break;
            } else {
              ret[transformations_it->first] = transformations_it->second;
            }
          }
        }
        if (ret.size() > 0) {
          return ret;
        }
        throw common::exception::TitledException("Operation expression expression '_pdbx_struct_assembly_gen.oper_expression' or the used format is not currently supported: " + ids);
      }

      public:
      static std::string parse_id(std::istream& input) {
        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          if (line.size() > 10 && common::string::starts_with(line, "_entry.id ")) {
            return parse_item(line, input, "_entry.id");
          }
        }
        throw common::exception::TitledException("No id line");
      }

      // Read pdb file from input filtered by filters and returns the corresponding protein
      //
      // TODO: Currently, there are parsed only aminoacids with ATOM/HETATM line, in the future, it will be better to use SEQRES for it.
      // However, if an aminoacid has no coordinates, it is totally unusefull for us. I.E. it should be coordinated with the validation part of processing.
      static Protein parse_protein(std::istream& input, BasicFilter* filter) {
        // Parsed protein to be returned
        Protein protein;
        // List of protein chains (to ignore other chains
        std::set<std::string> chains;
        // Crystallographic symmetry group to generate crystallographic transformations
        std::string space_group;
        // Assignation of biomolecule transformations to chains <biomolecule_id, <transformation, chains>>
        std::map<int, std::map<std::string, std::set<std::string>>> generation;
        // Whether we are in a loop block
        bool loop = false;
        // Index of item within a loop block
        size_t index = 0;
        // Currently parsed item of the file
        std::string line;

        // NOTE: input.eof() check is necessary for the case that the file does not ends with an empty line
        while (!input.eof() && std::getline(input, line)) {
          if (filter->keep(line)) {
            // Type of block
            if (common::string::starts_with(line, "loop_")) {
              loop = true;
              index = 0;
            // Protein identifier
            } else if (common::string::starts_with(line, "_entry.id ")) {
              if (loop) {
                throw common::exception::TitledException("Unexpected protein ID within a loop");
              }
              protein.ID_CODE = parse_item(line, input, "_entry.id ");
            // Identifiers of chains that should be parsed
            } else if (common::string::starts_with(line, "_entity_poly.pdbx_strand_id")) {
              if (loop) {
                size_t i = index;
                do {
                  ++index;
                } while (std::getline(input, line) && line.size() > 0 && line[0] == '_');
                while (!common::string::starts_with(line, "#")) {
                  std::stringstream values(extract_item(line, input, index, {i})[i]);
                  std::string value;
                  while (std::getline(values, value, ',')) {
                    chains.emplace(value);
                  }
                  if (!std::getline(input, line)) {
                    throw common::exception::TitledException("Unexpected end of file, missing closing sharp");
                  }
                }
              } else {
                std::stringstream values(parse_item(line, input, "_entity_poly.pdbx_strand_id"));
                std::string value;
                while (std::getline(values, value, ',')) {
                  chains.emplace(value);
                }
              }
            // Crystallographic symmetry
            } else if (common::string::starts_with(line, "_symmetry.space_group_name_H-M")) {
              if (loop) {
                throw common::exception::TitledException("Unexpected crystallographic space group within a loop");
              }
              space_group = parse_item(line, input, "_symmetry.space_group_name_H-M");
              // TODO: Dictionary of space groups and generate transformations
            // Biomolecules composition
            } else if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.")) {
              if (loop) {
                size_t id;
                size_t transformations;
                size_t affected;
                do {
                  if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.assembly_id")) {
                    id = index;
                  } else if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.oper_expression")) {
                    transformations = index;
                  } else if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.asym_id_list")) {
                    affected = index;
                  }
                  ++index;
                } while (std::getline(input, line) && line.size() > 0 && line[0] == '_');
                while (!common::string::starts_with(line, "#")) {
                  auto values = extract_item(line, input, index, {id, transformations, affected});
                  auto &subgeneration = generation[std::stoi(values[id])];
                  auto subgeneration_it = subgeneration.find(values[transformations]);
                  if (subgeneration_it == subgeneration.end()) {
                    subgeneration[values[transformations]] = parseChains(values[affected]);
                  } else {
                    auto parsed = parseChains(values[affected]);
                    subgeneration_it->second.insert(parsed.begin(), parsed.end());
                  }
                  if (!std::getline(input, line)) {
                    throw common::exception::TitledException("Unexpected end of file, missing closing sharp");
                  }
                }
              } else {
                int id;
                std::string transformations;
                std::string affected;
                do {
                  if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.assembly_id")) {
                    id = std::stoi(parse_item(line, input, "_pdbx_struct_assembly_gen.assembly_id"));
                  } else if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.oper_expression")) {
                    transformations = parse_item(line, input, "_pdbx_struct_assembly_gen.oper_expression");
                  } else if (common::string::starts_with(line, "_pdbx_struct_assembly_gen.asym_id_list")) {
                    affected = parse_item(line, input, "_pdbx_struct_assembly_gen.asym_id_list");
                  }
                } while (std::getline(input, line) && line.size() > 0 && line[0] != '#');
                if (id <= 0) {
                  throw common::exception::TitledException("Missing ID of biomolecule transformations (_pdbx_struct_assembly_gen.assembly_id)");
                }
                generation[id][transformations] = parseChains(affected);
              }
            // Biomolecules transformations
            } else if (common::string::starts_with(line, "_pdbx_struct_oper_list.")) {
              std::map<std::string, TransformationMatrix> matrices;
              if (loop) {
                size_t id_index;
                size_t matrix_indices[3][4];
                std::set<size_t> indices;
                do {
                  if (common::string::starts_with(line, "_pdbx_struct_oper_list.id")) {
                    id_index = index;
                    indices.emplace(index);
                  } else if (common::string::starts_with(line, "_pdbx_struct_oper_list.matrix[") && line[31] == ']' && line[32] == '[' && line[34] == ']') {
                    matrix_indices[line[30]-'1'][line[33]-'1'] = index;
                    indices.emplace(index);
                  } else if (common::string::starts_with(line, "_pdbx_struct_oper_list.vector[") && line[31] == ']') {
                    matrix_indices[line[30]-'1'][3] = index;
                    indices.emplace(index);
                  }
                  ++index;
                } while (std::getline(input, line) && line.size() > 0 && line[0] == '_');
                if (indices.size() != 13 ) {
                  throw common::exception::TitledException("Incorrectly deffined biomolecolue transformation matrix: _pdbx_struct_oper_list");
                }
                while (!common::string::starts_with(line, "#")) {
                  auto values = extract_item(line, input, index, indices);
                  TransformationMatrix matrix(std::make_tuple(std::make_tuple(std::stod(values[matrix_indices[0][0]]), std::stod(values[matrix_indices[0][1]]), std::stod(values[matrix_indices[0][2]]), std::stod(values[matrix_indices[0][3]])),
                                                              std::make_tuple(std::stod(values[matrix_indices[1][0]]), std::stod(values[matrix_indices[1][1]]), std::stod(values[matrix_indices[1][2]]), std::stod(values[matrix_indices[1][3]])),
                                                              std::make_tuple(std::stod(values[matrix_indices[2][0]]), std::stod(values[matrix_indices[2][1]]), std::stod(values[matrix_indices[2][2]]), std::stod(values[matrix_indices[2][3]]))));
                  matrices[values[id_index]] = matrix;
                  if (!std::getline(input, line)) {
                    throw common::exception::TitledException("Unexpected end of file, missing closing sharp");
                  }
                }
              } else {
                std::string id;
                TransformationMatrix matrix;
                do {
                  if (common::string::starts_with(line, "_pdbx_struct_oper_list.id")) {
                    id = parse_item(line, input, "_pdbx_struct_oper_list.id");
                  } else if (common::string::starts_with(line, "_pdbx_struct_oper_list.matrix[")) { // It is not possible to use variable in std::get<i>(tuple)
                    switch (line[30]) {
                      case '1':
                        switch (line[33]) {
                          case '1':
                            std::get<0>(std::get<0>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[1][1]"));
                            break;
                          case '2':
                            std::get<1>(std::get<0>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[1][2]"));
                            break;
                          case '3':
                            std::get<2>(std::get<0>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[1][3]"));
                            break;
                          default:
                            throw common::exception::TitledException("Unknown row index in biomolecule transformation amtrix: " + line);
                            break;
                        }
                        break;
                      case '2':
                        switch (line[33]) {
                          case '1':
                            std::get<0>(std::get<1>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[2][1]"));
                            break;
                          case '2':
                            std::get<1>(std::get<1>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[2][2]"));
                            break;
                          case '3':
                            std::get<2>(std::get<1>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[2][3]"));
                            break;
                          default:
                            throw common::exception::TitledException("Unknown row index in biomolecule transformation amtrix: " + line);
                            break;
                        }
                        break;
                      case '3':
                        switch (line[33]) {
                          case '1':
                            std::get<0>(std::get<2>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[3][1]"));
                            break;
                          case '2':
                            std::get<1>(std::get<2>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[3][2]"));
                            break;
                          case '3':
                            std::get<2>(std::get<2>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.matrix[3][3]"));
                            break;
                          default:
                            throw common::exception::TitledException("Unknown row index in biomolecule transformation amtrix: " + line);
                            break;
                        }
                        break;
                      default:
                        throw common::exception::TitledException("Unknown row index in biomolecule transformation amtrix: " + line);
                        break;
                    }
                  } else if (common::string::starts_with(line, "_pdbx_struct_oper_list.vector[")) {
                    switch (line[30]) {
                      case '1':
                        std::get<3>(std::get<0>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.vector[1]"));
                        break;
                      case '2':
                        std::get<3>(std::get<1>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.vector[2]"));
                        break;
                      case '3':
                        std::get<3>(std::get<2>(matrix)) = std::stod(parse_item(line, input, "_pdbx_struct_oper_list.vector[3]"));
                        break;
                      default:
                        throw common::exception::TitledException("Unknown row index in biomolecule transformation amtrix: " + line);
                        break;
                    }
                  }
                } while (std::getline(input, line) && line.size() > 0 && line[0] != '#');
                if (id.size() == 0) {
                  throw common::exception::TitledException("Missing ID of biomolecule transformations (_pdbx_struct_assembly_gen.assembly_id)");
                }
                matrices[id] = matrix;
              }
              for (auto generation_it = generation.begin(); generation_it != generation.end(); ++generation_it) {
                for (auto subgeneration_it = generation_it->second.begin(); subgeneration_it != generation_it->second.end(); ++subgeneration_it) {
                  std::string matrix_ids = subgeneration_it->first;
                  std::map<std::string, TransformationMatrix> transformations;
                  if (matrix_ids.size() == 0 || matrix_ids[0] != '(' || matrix_ids.back() != ')') {
                    transformations = parseTransformations(matrix_ids, matrices);
                  } else {
                    size_t end = matrix_ids.find(")(", 1);
                    if (end == matrix_ids.npos) {
                      transformations = parseTransformations(matrix_ids.substr(1, matrix_ids.size()-2), matrices);
                    } else {
                      size_t start = 1;
                      transformations = parseTransformations(matrix_ids.substr(1, end-1), matrices);
                      do {
                        start = end+2;
                        end = matrix_ids.find(")(", start);
                        if (end == matrix_ids.npos) {
                          end = matrix_ids.size()-1;
                        }
                        auto transformations2 = parseTransformations(matrix_ids.substr(start, end-start), matrices);
                        std::map<std::string, TransformationMatrix> transformations_next;
                        for (auto transformations_it = transformations.begin(); transformations_it != transformations.end(); ++transformations_it) {
                          for (auto transformations2_it = transformations2.begin(); transformations2_it != transformations2.end(); ++transformations2_it) {
                            TransformationMatrix matrix;
                            std::get<0>(std::get<0>(matrix)) = std::get<0>(std::get<0>(transformations_it->second))*std::get<0>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<0>(transformations_it->second))*std::get<0>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<0>(transformations_it->second))*std::get<0>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<0>(transformations_it->second));
                            std::get<1>(std::get<0>(matrix)) = std::get<0>(std::get<0>(transformations_it->second))*std::get<1>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<0>(transformations_it->second))*std::get<1>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<0>(transformations_it->second))*std::get<1>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<0>(transformations_it->second));
                            std::get<2>(std::get<0>(matrix)) = std::get<0>(std::get<0>(transformations_it->second))*std::get<2>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<0>(transformations_it->second))*std::get<2>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<0>(transformations_it->second))*std::get<2>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<0>(transformations_it->second));
                            std::get<3>(std::get<0>(matrix)) = std::get<0>(std::get<0>(transformations_it->second))*std::get<3>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<0>(transformations_it->second))*std::get<3>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<0>(transformations_it->second))*std::get<3>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<0>(transformations_it->second));
                            std::get<0>(std::get<1>(matrix)) = std::get<0>(std::get<1>(transformations_it->second))*std::get<0>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<1>(transformations_it->second))*std::get<0>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<1>(transformations_it->second))*std::get<0>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<1>(transformations_it->second));
                            std::get<1>(std::get<1>(matrix)) = std::get<0>(std::get<1>(transformations_it->second))*std::get<1>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<1>(transformations_it->second))*std::get<1>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<1>(transformations_it->second))*std::get<1>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<1>(transformations_it->second));
                            std::get<2>(std::get<1>(matrix)) = std::get<0>(std::get<1>(transformations_it->second))*std::get<2>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<1>(transformations_it->second))*std::get<2>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<1>(transformations_it->second))*std::get<2>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<1>(transformations_it->second));
                            std::get<3>(std::get<1>(matrix)) = std::get<0>(std::get<1>(transformations_it->second))*std::get<3>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<1>(transformations_it->second))*std::get<3>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<1>(transformations_it->second))*std::get<3>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<1>(transformations_it->second));
                            std::get<0>(std::get<2>(matrix)) = std::get<0>(std::get<2>(transformations_it->second))*std::get<0>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<2>(transformations_it->second))*std::get<0>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<2>(transformations_it->second))*std::get<0>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<2>(transformations_it->second));
                            std::get<1>(std::get<2>(matrix)) = std::get<0>(std::get<2>(transformations_it->second))*std::get<1>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<2>(transformations_it->second))*std::get<1>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<2>(transformations_it->second))*std::get<1>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<2>(transformations_it->second));
                            std::get<2>(std::get<2>(matrix)) = std::get<0>(std::get<2>(transformations_it->second))*std::get<2>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<2>(transformations_it->second))*std::get<2>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<2>(transformations_it->second))*std::get<2>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<2>(transformations_it->second));
                            std::get<3>(std::get<2>(matrix)) = std::get<0>(std::get<2>(transformations_it->second))*std::get<3>(std::get<0>(transformations2_it->second))+
                              std::get<1>(std::get<2>(transformations_it->second))*std::get<3>(std::get<1>(transformations2_it->second))+
                              std::get<2>(std::get<2>(transformations_it->second))*std::get<3>(std::get<2>(transformations2_it->second))+
                              std::get<3>(std::get<2>(transformations_it->second));
                            transformations_next[transformations_it->first + "." + transformations2_it->first] = matrix;
                          }
                        }
                        transformations = transformations_next;
                      } while (end < matrix_ids.size()-1);
                    }
                  }
                  auto &biomolecule = protein.BIOMOLECULES[generation_it->first];
                  for (auto transformations_it = transformations.begin(); transformations_it != transformations.end(); ++transformations_it) {
                    auto trans_it = biomolecule.TRANSFORMATIONS.find(transformations_it->first);
                    if (trans_it == biomolecule.TRANSFORMATIONS.end()) {
                      biomolecule.TRANSFORMATIONS[transformations_it->first] = {transformations_it->second, subgeneration_it->second};
                    } else {
                      trans_it->second.second.insert(subgeneration_it->second.begin(), subgeneration_it->second.end());
                    }
                  }
                }
              }
            // Atoms
            } else if (common::string::starts_with(line, "_atom_site.")) {
              if (chains.size() == 0) {
                throw common::exception::TitledException("Block '_entity_poly' missing or occures after block '_atom_site' which is not currently supported");
              }
              if (!loop) {
                throw common::exception::TitledException("Single atom definition is not supported");
              }
              size_t model_index = -1;
              size_t chain_index = -1;
              size_t chain_index2 = -1;
              size_t residue_index = -1;
              size_t residue_index2 = -1;
              size_t insertion_index = -1;
              size_t residue_type_index = -1;
              size_t residue_type_index2 = -1;
              size_t residue_group_index = -1;
              size_t atom_index = -1;
              size_t atom_name_index = -1;
              size_t atom_name_index2 = -1;
              size_t atom_type_index = -1;
              size_t altloc_index = -1;
              size_t x_index = -1;
              size_t y_index = -1;
              size_t z_index = -1;
              size_t occupancy_index = -1;
              size_t temperature_index = -1;
              size_t charge_index = -1;
              std::set<size_t> indices;
              do {
                if (common::string::starts_with(line, "_atom_site.id")) {
                  atom_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.auth_asym_id")) {
                  chain_index2 = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.auth_atom_id")) {
                  atom_name_index2 = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.auth_comp_id")) {
                  residue_type_index2 = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.auth_seq_id")) {
                  residue_index2 = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.B_iso_or_equiv")) {
                  temperature_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.Cartn_x")) {
                  x_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.Cartn_y")) {
                  y_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.Cartn_z")) {
                  z_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.group_PDB")) {
                  residue_group_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.label_alt_id")) {
                  altloc_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.label_asym_id")) {
                  chain_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.label_atom_id")) {
                  atom_name_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.label_comp_id")) {
                  residue_type_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.label_seq_id")) {
                  residue_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.occupancy")) {
                  occupancy_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.pdbx_formal_charge")) {
                  charge_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.pdbx_PDB_ins_code")) {
                  insertion_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.pdbx_PDB_model_num")) {
                  model_index = index;
                  indices.emplace(index);
                } else if (common::string::starts_with(line, "_atom_site.type_symbol")) {
                  atom_type_index = index;
                  indices.emplace(index);
                }
                ++index;
              } while (std::getline(input, line) && line.size() > 0 && line[0] == '_');
              if (indices.size() != 20) {
                throw common::exception::TitledException("Incorrectly deffined atom headers: _atom_sites");
              }
              if (chain_index == -1) {
                if (chain_index2 == -1) {
                  throw common::exception::TitledException("Chain identifier is not specified");
                } else {
                  chain_index = chain_index2;
                }
              } else {
                if (chain_index2 == -1) {
                  chain_index2 = chain_index;
                }
              }
              if (residue_index == -1) {
                if (residue_index2 == -1) {
                  throw common::exception::TitledException("Residue identifier is not specified");
                } else {
                  residue_index = residue_index2;
                }
              } else {
                if (residue_index2 == -1) {
                  residue_index2 = residue_index;
                }
              }
              if (residue_type_index == -1) {
                if (residue_type_index2 == -1) {
                  throw common::exception::TitledException("Residue type identifier identifier is not specified");
                } else {
                  residue_type_index = residue_type_index2;
                }
              } else {
                if (residue_type_index2 == -1) {
                  residue_type_index2 = residue_type_index;
                }
              }
              if (atom_name_index == -1) {
                if (atom_name_index2 == -1) {
                  throw common::exception::TitledException("Residue type identifier identifier is not specified");
                } else {
                  residue_type_index = atom_name_index2;
                }
              } else {
                if (atom_name_index2 == -1) {
                  atom_name_index2 = atom_name_index;
                }
              }

              std::pair<const int, Model>* model = nullptr;
              std::pair<const std::string, Chain>* chain = nullptr;
              std::pair<const std::pair<int, char>, Aminoacid>* aminoacid = nullptr;
              bool missing_model = false;
              std::set<std::string> warnings;
              while (!common::string::starts_with(line, "#")) {
                auto values = extract_item(line, input, index, indices);
                std::string chain_id2 = values[chain_index2];
                std::string residue = values[residue_index];
                if (chains.find(chain_id2) != chains.end() && residue.size() > 0) {
                  // Set the current model
                  size_t model_id = (model_index == -1) ? 1 : std::stol(values[model_index]);
                  if (model == nullptr || model_id != model->first) {
                    auto ins = protein.MODELS.insert({model_id, Model()});
                    if (!ins.second) {
                      throw common::exception::TitledException("Unexpected format of mmCIF file: multiple models with the same serial number: '" + std::to_string(model_id) + "'");
                    }
                    if (!missing_model && model_id != protein.MODELS.size()) {
                      missing_model = true;
                      std::cerr << "WARNING [" << protein.ID_CODE << "]:    Missing models with serial numbers from " << protein.MODELS.size() << " to " << model_id-1 << std::endl;
                    }
                    model = &*ins.first;
                    chain = nullptr;
                    aminoacid = nullptr;
                  }
                  // Set the current chain
                  std::string chain_id = values[chain_index];
                  //It seems, mmcif require each item have a different label, while in pdb files water and ionts have the same label
                  if (chain_id2 != chain_id) {
                    warnings.emplace("Different values of '_atom_site.auth_asym_id' and '_atom_site.label_asym_id' are not currently supported");
                  }
                  if (chain == nullptr || chain_id != chain->first) {
                    auto ins = model->second.CHAINS.insert({chain_id, Chain()});
                    if (!ins.second) {
                      throw common::exception::TitledException("Multiple chains with the same identifier: '" + chain_id + "'");
                    }
                    chain = &*ins.first;
                    aminoacid = nullptr;
                  }
                  // Set the current aminoacid
                  if (residue != values[residue_index2]) {
                    warnings.emplace("Different values of '_atom_site.auth_seq_id' and '_atom_site.label_seq_id' are not currently supported");
                  }
                  if (values[residue_type_index] != values[residue_type_index2]) {
                    throw common::exception::TitledException("Different values of '_atom_site.auth_comp_id' and '_atom_site.label_comp_id' are not currently supported");
                  }
                  std::string insertion = insertion_index == -1 ? " " : values[insertion_index];
                  if (insertion == "?" || insertion.size() == 0) {
                    insertion = " ";
                  } else if (insertion.size() > 1) {
                    throw common::exception::TitledException("Alternative location identifier should be a character, not a string: '" + insertion + "'");
                  }
                  const std::pair<int, char> aa_id = {std::stoi(values[residue_index]), insertion[0]};
                  if (aminoacid == nullptr || aa_id != aminoacid->first) {
                    std::pair<std::map<std::pair<int, char>, Aminoacid>::iterator, bool> ins;
                    if (values[residue_group_index] == "ATOM") {
                      ins = chain->second.AMINOACIDS.insert({aa_id, Aminoacid(values[residue_type_index])});
                    } else if (values[residue_group_index] == "HETATM") {
                      ins = chain->second.AMINOACIDS.insert({aa_id, ModifiedAminoacid(values[residue_type_index])});
                    } else {
                      throw common::exception::TitledException("UNEXPECTED type of atom group: '" + values[residue_type_index] + "'");
                    }
                    if (!ins.second) {
                      throw common::exception::TitledException("Multiple aminoacids with the same identifier: '" + std::to_string(aa_id.first) + aa_id.second + "'");
                    }
                    aminoacid = &*ins.first;
                  }
                  // Set the current atom
                  std::string atom_name = values[atom_name_index];
                  if (atom_name != values[atom_name_index2]) {
                    throw common::exception::TitledException("Different values of '_atom_site.auth_atom_id' and '_atom_site.label_atom_id' are not currently supported");
                  }
                  std::string altloc = altloc_index == -1 ? " " : values[altloc_index];
                  if (altloc == "." || altloc.size() == 0) {
                    altloc = " ";
                  } else if (altloc.size() > 1) {
                    throw common::exception::TitledException("Alternative location identifier should be a character, not a string: '" + altloc + "'");
                  }
                  Atom& atom = aminoacid->second.ATOMS[atom_name];
                  if (altloc == " ") {
                    if (atom.ALTERNATIVE_LOCATIONS.size() > 0) {
                      throw common::exception::TitledException("Multiple ATOM lines for the same atom without altLoc specifier: '" + altloc + "'");
                    }
                    atom.ELEMENT = values[atom_type_index];
                  } else {
                    if (atom.ALTERNATIVE_LOCATIONS.empty()) {
                      atom.ELEMENT = values[atom_type_index];
                    } else if (atom.ELEMENT != values[atom_type_index]) {
                      throw common::exception::TitledException("Multiple ATOM lines for the same atom contains differents elements: '" + atom.ELEMENT + " vs. " + values[atom_type_index] + "'");
                    }
                  }
                  // Set the current coordinate
                  auto coordinate = atom.ALTERNATIVE_LOCATIONS.insert({altloc[0], Characteristic()});
                  if (!coordinate.second) {
                    throw common::exception::TitledException("Multiple alternative locations with the same identifier: '" + line + "'");
                  }
                  coordinate.first->second.SERIAL_NUMBER = std::stoi(values[atom_index]);
                  coordinate.first->second.X = std::stod(values[x_index]);
                  coordinate.first->second.Y = std::stod(values[y_index]);
                  coordinate.first->second.Z = std::stod(values[z_index]);
                  coordinate.first->second.OCCUPANCY = std::stof(values[occupancy_index]);
                  coordinate.first->second.TEMPERATURE = std::stof(values[temperature_index]);
                  std::string charge = values[charge_index];
                  coordinate.first->second.CHARGE = charge.size() == 0 ? 0 : std::stoi(charge);
                }
                if (!std::getline(input, line)) {
                  throw common::exception::TitledException("Unexpected end of file, missing closing sharp");
                }
              }
              for (auto warnings_it = warnings.begin(); warnings_it != warnings.end(); ++warnings_it) {
                //std::cerr << "WARNING: " << *warnings_it << std::endl;
              }
            }
            if (loop && line.size() > 0 && line[0] == '_') {
              ++index;
            }
            if (common::string::starts_with(line, "#")) {
              loop = false;
              index = 0;
            }
          }
        }

        if (protein.BIOMOLECULES.empty()) {
          std::cerr << "WARNING [" << protein.ID_CODE << "]:    Protein does not contain a definition of biomolecules (REMARK 350)" << std::endl;
          Biomolecule &biomolecule = protein.BIOMOLECULES[1];
          auto &transformation = biomolecule.TRANSFORMATIONS["1"];
          std::get<0>(std::get<0>(transformation.first)) = 1;
          std::get<1>(std::get<1>(transformation.first)) = 1;
          std::get<2>(std::get<2>(transformation.first)) = 1;
          transformation.second = chains;
        }

        if (protein.CRYSTALLOGRAPHIC_SYMMETRIES.empty()) {
          //std::cerr << "WARNING [" << protein.ID_CODE << "]:    Protein does not contain a definition of crystallographic transformation (REMARK 290)" << std::endl;
          TransformationMatrix &transformation = protein.CRYSTALLOGRAPHIC_SYMMETRIES[1];
          std::get<0>(std::get<0>(transformation)) = 1;
          std::get<1>(std::get<1>(transformation)) = 1;
          std::get<2>(std::get<2>(transformation)) = 1;
        }

        return protein;
      }
    };
  }
}