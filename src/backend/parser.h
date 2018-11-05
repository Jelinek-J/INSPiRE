#pragma once

#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include "../elemental/string.h"
#include "protein.h"
#include "pdb.h"
#include "mmcif.h"
#include "xml.h"

namespace inspire {
  namespace backend {
    struct ProteinParser {
      static std::string parse_id(const std::string& input) {
        std::ifstream stream(input);
        if (elemental::string::ends_with(input, ".pdb")) {
          return Pdb::parse_id(stream);
        } else if (elemental::string::ends_with(input, ".cif")) {
          return Mmcif::parse_id(stream);
        } else if (elemental::string::ends_with(input, ".xml")) {
#if defined(_DEBUG) && !defined(BIGOBJ)
          throw elemental::exception::TitledException(input + ": PDBML file format is not currently supported"); 
#else
          return Xml::parse_id(stream);
#endif
        } else {
// TODO: Try to recognize format based on the first line
          throw elemental::exception::TitledException(input + ": Unknown file format is not currently supported");
        }
      }

      static Protein parse_protein(std::string& input, BasicFilter* filter) {
        std::ifstream stream(input);
        if (elemental::string::ends_with(input, ".pdb")) {
          return Pdb::parse_protein(stream, filter);
        } else if (elemental::string::ends_with(input, ".cif")) {
          return Mmcif::parse_protein(stream, filter);
        } else if (elemental::string::ends_with(input, ".xml")) {
#if defined(_DEBUG) && !defined(BIGOBJ)
          throw elemental::exception::TitledException(input + ": PDBML file format is not currently supported");
#else
          return Xml::parse_protein(stream, filter);
#endif
        } else {
// TODO: Try to recognize format based on the first line
          throw elemental::exception::TitledException(input + ": Unknown file format is not currently supported");
        }
      }
    };
  }
}