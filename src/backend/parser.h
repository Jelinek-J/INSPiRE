#pragma once

#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include "../common/string.h"
#include "protein.h"
#include "pdb.h"
#include "mmcif.h"
#include "xml.h"

namespace inspire {
  namespace backend {
    struct ProteinParser {
      static std::string parse_id(const std::string& input) {
        std::ifstream stream(input);
        try {
          if (common::string::ends_with(input, ".pdb") || common::string::ends_with(input, ".ent")) {
            return Pdb::parse_id(stream);
          } else if (common::string::ends_with(input, ".cif")) {
            return Mmcif::parse_id(stream);
          } else if (common::string::ends_with(input, ".xml")) {
            return Xml::parse_id(stream);
          } else {
  // TODO: Try to recognize format based on the first line
            throw common::exception::TitledException("Unknown file format is not currently supported");
          }
        } catch (const common::exception::TitledException& e) {
          throw common::exception::TitledException("[" + input + "] " + e.what());
        }
      }

      static bool contains_hetatoms(const std::string& input) {
        std::ifstream stream(input);
        try {
          if (common::string::ends_with(input, ".pdb") || common::string::ends_with(input, ".ent")) {
            return Pdb::contains_hetatoms(stream);
          } else if (common::string::ends_with(input, ".cif")) {
            return Mmcif::contains_hetatoms(stream);
          } else if (common::string::ends_with(input, ".xml")) {
            return Xml::contains_hetatoms(stream);
          } else {
  // TODO: Try to recognize format based on the first line
            throw common::exception::TitledException("Unknown file format is not currently supported");
          }
        } catch (const common::exception::TitledException& e) {
          throw common::exception::TitledException("[" + input + "] " + e.what());
        }
      }

      static Protein parse_protein(std::string& input, BasicFilter* filter) {
        std::ifstream stream(input);
        try {
          if (common::string::ends_with(input, ".pdb") || common::string::ends_with(input, ".ent")) {
            return Pdb::parse_protein(stream, filter);
          } else if (common::string::ends_with(input, ".cif")) {
            return Mmcif::parse_protein(stream, filter);
          } else if (common::string::ends_with(input, ".xml")) {
            return Xml::parse_protein(stream, filter);
          } else {
  // TODO: Try to recognize format based on the first line
            throw common::exception::TitledException("Unknown file format is not currently supported");
          }
        } catch (const common::exception::TitledException& e) {
          throw common::exception::TitledException("[" + input + "] " + e.what());
        }
      }
    };
  }
}
