#pragma once

#include "protein.h"
#include "filters.h"
#include "../common/string.h"
#include "../common/exception.h"
#include "../common/xml.h"
#include <vector>
#include <set>
#include <istream>
#include <sstream>
#include <iostream>

// This file specializes on reading PDB files http://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html

namespace inspire {
  namespace backend {
    // A class grouping methods for parsing from pdb file
    struct Xml {
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
        // Load xml
        common::xml::Xml xml(input);
        // id
        std::string id;
        // Protein identifier
        if (!xml.get_attribute("PDBx:datablock.PDBx:entryCategory.PDBx:entry", "id", id)) {
          throw common::exception::TitledException("The file does not contain any identifier");
        }
        return id;
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
        // Alternative names of protein chains
        std::map<std::string, std::string> alternatives;
        // List of residues to make validation of defined atoms
        std::map<std::string, std::vector<std::string>> entities;
        // List of assemblies
        std::set<size_t> assemblies;
        // Crystallographic symmetry group to generate crystallographic transformations
        std::string space_group;
        // Assignation of biomolecule transformations to chains <biomolecule_id, <transformation, chains>>
        std::map<int, std::map<std::string, std::set<std::string>>> generation;
        // Load xml
        common::xml::Xml xml(input);

        // Protein identifier
        if (!xml.get_attribute("PDBx:datablock.PDBx:entryCategory.PDBx:entry", "id", protein.ID_CODE)) {
          throw common::exception::TitledException("The file does not contain any identifier");
        }

        // Identifiers of chains that should be parsed
        std::vector<common::xml::Xml> nodes;
        if (xml.get_nodes("PDBx:datablock.PDBx:entity_polyCategory", "PDBx:entity_poly", nodes)) {
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string value;
            if (!nodes_it->get_value("PDBx:pdbx_strand_id", value)) {
              throw common::exception::TitledException("The file does not contain any identifier");
            }
            std::stringstream parts(value);
            std::string part;
            while (std::getline(parts, part, ',')) {
              chains.emplace(part);
            }
          }
          if (chains.size() == 0) {
            throw common::exception::TitledException("Block 'entity_poly' does not contain any valid chain identifier");
          }
        } else {
          throw common::exception::TitledException("Block 'entity_poly' is missing");
        }

        // Alternative identifiers of chains that should be parsed
        if (xml.get_nodes("PDBx:datablock.PDBx:struct_asymCategory", "PDBx:struct_asym", nodes)) {
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string id;
            std::string value;
            if (!(nodes_it->get_attribute("id", id) && nodes_it->get_value("PDBx:entity_id", value))) {
              throw common::exception::TitledException("Incomplete definition of '_struct_asym'");
            }
            alternatives[id] = value;
          }
          if (chains.size() == 0) {
            throw common::exception::TitledException("Block 'struct_asym' does not contain any valid chain identifier");
          }
        } else {
          throw common::exception::TitledException("Block 'struct_asym' is missing");
        }

        // Composition of parsed entities
        if (xml.get_nodes("PDBx:datablock.PDBx:entity_poly_seqCategory", "PDBx:entity_poly_seq", nodes)) {
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string id;
            std::string monomer;
            std::string number;
            if (!(nodes_it->get_attribute("entity_id", id) && nodes_it->get_attribute("mon_id", monomer) && nodes_it->get_attribute("num", number))) {
              throw common::exception::TitledException("Incomplete definition of '_entity_poly_seq'");
            }
            auto &chain = entities[id];
            chain.push_back(monomer);
            if (chain.size() != std::stol(number)) {
              throw common::exception::TitledException("Entity polymer numbers '_entity_poly_seq.num' are not sequential");
            }
          }
          if (chains.size() == 0) {
            throw common::exception::TitledException("Block 'entity_poly_seq' does not contain any valid chain identifier");
          }
        } else {
          throw common::exception::TitledException("Block 'entity_poly_seq' is missing");
        }

        // Crystallographic symmetry
        if (!xml.get_value("PDBx:datablock.PDBx:symmetryCategory.PDBx:symmetry.PDBx:space_group_name_H-M", space_group)) {
          protein.CRYSTALLOGRAPHIC_SYMMETRIES[1] = TransformationMatrix();
        }
        // TODO: Dictionary of space groups and generate transformations;

        // Biomolecules composition
        if (xml.get_nodes("PDBx:datablock.PDBx:pdbx_struct_assemblyCategory", "PDBx:pdbx_struct_assembly", nodes)) {
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string id;
            std::string details;
            if (!(nodes_it->get_attribute("id", id) && nodes_it->get_value("PDBx:details", details))) {
              throw common::exception::TitledException("Incomplete definition of 'pdbx_struct_assembly_gen'");
            }
            if (common::string::ends_with(common::string::to_upper(details), "ASSEMBLY")) {
              assemblies.emplace(std::stol(id));
            }
          }
        }
        if (xml.get_nodes("PDBx:datablock.PDBx:pdbx_struct_assembly_genCategory", "PDBx:pdbx_struct_assembly_gen", nodes)) {
          if (assemblies.empty()) {
            throw common::exception::TitledException("Data category '_pdbx_struct_assembly' is missing, or it occurs after '_pdbx_struct_assembly_gen', or no valid assembly is defined (~ missing 'REMARKS 350')");
          }
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string id;
            std::string affected;
            std::string transformations;
            if (!(nodes_it->get_attribute("assembly_id", id) && nodes_it->get_attribute("asym_id_list", affected) && nodes_it->get_attribute("oper_expression", transformations))) {
              throw common::exception::TitledException("Incomplete definition of 'pdbx_struct_assembly_gen'");
            }
            size_t assembly = std::stoi(id);
            if (assemblies.find(assembly) != assemblies.end()) {
              auto &subgeneration = generation[assembly];
              auto subgeneration_it = subgeneration.find(transformations);
              if (subgeneration_it == subgeneration.end()) {
                subgeneration[transformations] = parseChains(affected);
              } else {
                auto parsed = parseChains(affected);
                subgeneration_it->second.insert(parsed.begin(), parsed.end());
              }
            }
          }
        }

        // Biomolecules transformations
        if (xml.get_nodes("PDBx:datablock.PDBx:pdbx_struct_oper_listCategory", "PDBx:pdbx_struct_oper_list", nodes)) {
          std::map<std::string, TransformationMatrix> matrices;
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string id;
            std::string cell;
            double matrix[3][4];
            if (!nodes_it->get_attribute("id", id)) {
              throw common::exception::TitledException("Missing identifier of 'pdbx_struct_oper_list'");
            }
            for (size_t i = 0; i < 3; i++) {
              for (size_t j = 0; j < 3; j++) {
                if (!nodes_it->get_value("PDBx:matrix" + std::to_string(i+1) + std::to_string(j+1), cell)) {
                  throw common::exception::TitledException("Incomplete definition of 'pdbx_struct_oper_list', missing 'matrix" + std::to_string(i+1) + std::to_string(j+1) + "'");
                }
                matrix[i][j] = std::stod(cell);
              }
              if (!nodes_it->get_value("PDBx:vector" + std::to_string(i+1), cell)) {
                throw common::exception::TitledException("Incomplete definition of 'pdbx_struct_oper_list', missing 'vector" + std::to_string(i+1) + "'");
              }
              matrix[i][3] = std::stod(cell);
            }
            TransformationMatrix m(std::make_tuple(std::make_tuple(matrix[0][0], matrix[0][1], matrix[0][2], matrix[0][3]),
                                                   std::make_tuple(matrix[1][0], matrix[1][1], matrix[1][2], matrix[1][3]),
                                                   std::make_tuple(matrix[2][0], matrix[2][1], matrix[2][2], matrix[2][3])));
            matrices[id] = m;
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
        }

        // Atoms etc.
        if (xml.get_nodes("PDBx:datablock.PDBx:atom_siteCategory", "PDBx:atom_site", nodes)) {
          std::pair<const int, Model>* model = nullptr;
          std::pair<const std::string, Chain>* chain = nullptr;
          std::pair<const std::pair<int, char>, Aminoacid>* aminoacid = nullptr;
          bool missing_model = false;
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string chain_index;
            std::string chain_index2;
            std::string residue_index;
            std::string residue_index2;
            if (!((nodes_it->get_value("PDBx:label_asym_id", chain_index)  | nodes_it->get_value("PDBx:auth_asym_id", chain_index2)) &&
                  (nodes_it->get_value("PDBx:label_seq_id", residue_index) | nodes_it->get_value("PDBx:auth_seq_id", residue_index2)))) {
              throw common::exception::TitledException("Incomplete definition of 'atom_site'");
            }
            if (chain_index == "") {
              chain_index = chain_index2;
            } else if (chain_index2 == "") {
              chain_index2 = chain_index;
            }
            if (residue_index2 == "") {
              residue_index2 = residue_index;
            }
            if (chains.find(chain_index2) != chains.end() && residue_index.size() > 0) {
              std::string model_index;
              std::string insertion;
              std::string residue_type;
              std::string residue_type2;
              std::string residue_group;
              std::string atom_index;
              std::string atom_name;
              std::string atom_name2;
              std::string atom_type;
              std::string altloc;
              std::string x;
              std::string y;
              std::string z;
              std::string occupancy;
              std::string temperature;
              std::string charge;
              if (!(nodes_it->get_attribute("id", atom_index) && (nodes_it->get_value("PDBx:label_comp_id", residue_type) | nodes_it->get_value("PDBx:auth_comp_id", residue_type2)) &&
                                                                 (nodes_it->get_value("PDBx:label_atom_id", atom_name)    | nodes_it->get_value("PDBx:auth_atom_id", atom_name2)))) {
                throw common::exception::TitledException("Incomplete definition of 'atom_site'");
              }
              if (residue_type == "") {
                residue_type = residue_type2;
              } else if (residue_type2 == "") {
                residue_type2 = residue_type;
              }
              if (atom_name == "") {
                atom_name = atom_name2;
              } else if (atom_name2 == "") {
                atom_name2 = atom_name;
              }
              nodes_it->get_value("PDBx:pdbx_formal_charge", charge);
              nodes_it->get_value("PDBx:B_iso_or_equiv", temperature);
              nodes_it->get_value("PDBx:Cartn_x", x);
              nodes_it->get_value("PDBx:Cartn_y", y);
              nodes_it->get_value("PDBx:Cartn_z", z);
              nodes_it->get_value("PDBx:group_PDB", residue_group);
              nodes_it->get_value("PDBx:label_alt_id", altloc);
              nodes_it->get_value("PDBx:occupancy", occupancy);
              nodes_it->get_value("PDBx:pdbx_PDB_ins_code", insertion);
              nodes_it->get_value("PDBx:pdbx_PDB_model_num", model_index);
              nodes_it->get_value("PDBx:type_symbol", atom_type);

              // Set the current model
              size_t model_id = (model_index == "") ? 1 : std::stol(model_index);
              if (model == nullptr || model_id != model->first) {
                auto ins = protein.MODELS.insert({model_id, Model()});
                if (!ins.second) {
                  throw common::exception::TitledException("Unexpected format of XML file: multiple models with the same serial number: '" + std::to_string(model_id) + "'");
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
              //It seems, xml require each item have a different label, while in pdb files water and ionts have the same label
              if (chain == nullptr || chain_index != chain->first) {
                auto ins = model->second.CHAINS.insert({chain_index, Chain()});
                if (!ins.second) {
                  throw common::exception::TitledException("Multiple chains with the same identifier: '" + chain_index + "'");
                }
                chain = &*ins.first;
                aminoacid = nullptr;
              }

              // Set the current aminoacid
              if (residue_type != residue_type2) {
                throw common::exception::TitledException("Different values of '_atom_site.auth_comp_id' and '_atom_site.label_comp_id' are not currently supported");
              }
              if (insertion == "?" || insertion.size() == 0) {
                insertion = " ";
              } else if (insertion.size() > 1) {
                throw common::exception::TitledException("Alternative location identifier should be a character, not a string: '" + insertion + "'");
              }
              const std::pair<int, char> aa_id = {std::stoi(residue_index), insertion[0]};
              if (aminoacid == nullptr || aa_id != aminoacid->first) {
                std::pair<std::map<std::pair<int, char>, Aminoacid>::iterator, bool> ins;
                if (residue_group == "ATOM") {
                  ins = chain->second.AMINOACIDS.insert({aa_id, Aminoacid(residue_type)});
                } else if (residue_group == "HETATM") {
                  ins = chain->second.AMINOACIDS.insert({aa_id, ModifiedAminoacid(residue_type)});
                } else {
                  throw common::exception::TitledException("UNEXPECTED type of atom group: '" + residue_type + "'");
                }
                if (!ins.second) {
                  throw common::exception::TitledException("Multiple aminoacids with the same identifier: '" + std::to_string(aa_id.first) + aa_id.second + "'");
                }
                aminoacid = &*ins.first;
              }

              // Set the current atom
              if (atom_name != atom_name2) {
                throw common::exception::TitledException("Different values of '_atom_site.auth_atom_id' and '_atom_site.label_atom_id' are not currently supported");
              }
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
                atom.ELEMENT = atom_type;
              } else {
                if (atom.ALTERNATIVE_LOCATIONS.empty()) {
                  atom.ELEMENT = atom_type;
                } else if (atom.ELEMENT != atom_type) {
                  throw common::exception::TitledException("Multiple ATOM lines for the same atom contains differents elements: '" + atom.ELEMENT + " vs. " + atom_type + "'");
                }
              }
              // Set the current coordinate
              auto coordinate = atom.ALTERNATIVE_LOCATIONS.insert({altloc[0], Characteristic()});
              if (!coordinate.second) {
                throw common::exception::TitledException("Multiple alternative locations with the same identifier: '<PDBx:atom_site id=\"" + atom_index + "\">'");
              }
              coordinate.first->second.SERIAL_NUMBER = std::stoi(atom_index);
              coordinate.first->second.X = std::stod(x);
              coordinate.first->second.Y = std::stod(y);
              coordinate.first->second.Z = std::stod(z);
              coordinate.first->second.OCCUPANCY = std::stof(occupancy);
              coordinate.first->second.TEMPERATURE = std::stof(temperature);
              coordinate.first->second.CHARGE = charge.size() == 0 ? 0 : std::stoi(charge);
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

        for (auto models_it = protein.MODELS.begin(); models_it != protein.MODELS.end(); ++models_it) {
          for (auto chains_it = models_it->second.CHAINS.begin(); chains_it != models_it->second.CHAINS.end(); ++chains_it) {
            auto alternatives_it = alternatives.find(chains_it->first);
            if (alternatives_it == alternatives.end()) {
              std::cerr << "WARNING [" << protein.ID_CODE << "]:    Protein does not contain chain '" << chains_it->first << "' in model '" << models_it->first << "'" << std::endl;
            } else {
              auto entity = entities.find(alternatives_it->second);
              if (entity == entities.end()) {
                throw common::exception::TitledException("Entity '" + alternatives_it->second + "' is declared in '_entity_poly' section but not defined in '_entity_poly_seq' section'");
              }
              if (chains_it->second.AMINOACIDS.size() != entity->second.size()) {
                if (chains_it->second.AMINOACIDS.size() < entity->second.size()) {
                  std::cerr << "WARNING [" << protein.ID_CODE << "]:    Chain '" << chains_it->first << "' in model '" << models_it->first << "' does not contain information about " << entity->second.size() - chains_it->second.AMINOACIDS.size() << " residues" << std::endl;
                } else {
                  std::cerr << "WARNING [" << protein.ID_CODE << "]:    Chain '" << chains_it->first << "' in model '" << models_it->first << "' has defined more residues in ATOM/HETATM section than in SEQRES section" << std::endl;
                }
              }
            }
          }
        }

        return protein;
      }

      static bool contains_hetatoms(std::istream& input) {
        // Load xml
        common::xml::Xml xml(input);
        std::vector<common::xml::Xml> nodes;
        if (xml.get_nodes("PDBx:datablock.PDBx:atom_siteCategory", "PDBx:atom_site", nodes)) {
          for (auto nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it) {
            std::string residue_group;
            if (nodes_it->get_value("PDBx:group_PDB", residue_group) && residue_group == "HETATM") {
              std::string residues[2];
              if (nodes_it->get_value("PDBx:label_comp_id", residues[0]) && nodes_it->get_value("PDBx:auth_comp_id", residues[1])) {
                char situation = 0;
                for (size_t i = 0; i < 2; i++) {
                  if (residues[i] == "HOH" || residues[i] == "DOD" || residues[i] == "DIS" || residues[i] == "MTO") {
                    ++situation;
                  }
                }
                switch (situation) {
                  case 0: // Contains non-water heteroatom
                    return true;
                    break;
                  case 1: // One label claims it is a water, the other one claims it is not a water
                    throw common::exception::TitledException("Ambiguity between heteroatom labels: '" + residues[0] + "' vs. '" + residues[1]);
                  case 2:
                    break;
                  default:
                    throw common::exception::TitledException("This should not happen: ++(++(0)) <= 2");
                }
              }
            }
          }
        }
        return false;
      }
    };
  }
}
