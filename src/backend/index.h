#pragma once

#include "iterators.h"
#include "../elemental/filesystem.h"
#include "../elemental/string.h"
#include "protein.h"
#include <string>
#include <sstream>
#include <iostream>

namespace inspire {
  namespace backend {
    class Index {
      private:
      std::ifstream FILE;
      size_t INDEX;
      std::string PROTEIN;
      std::string MODEL;
      std::string CHAIN;
      std::string AMINOACID;

      public:
      Index(std::string filename) : FILE(filename) { }

      // TODO: I expect this is more effective than create and return clumping object
      bool next() {
        std::string line;
        if (FILE.eof() || !std::getline(FILE, line)) {
          return false;
        }
        ++INDEX;
        std::stringstream elements(line);
        std::string element;
        if (!std::getline(elements, element, '\t')) {
          throw elemental::exception::TitledException("Unexpected empty line in the index file.");
        }
        AMINOACID = element;
        if (std::getline(elements, element, '\t')) {
          CHAIN = element;
          if (std::getline(elements, element, '\t')) {
            MODEL = element;
            if (std::getline(elements, element, '\t')) {
              PROTEIN = element;
            }
          }
        }
        return true;
      }
      bool reset() {
        // TODO: Is it better to needlessly reset in the first reading, or test whether it is necessary in every reading?
        FILE.clear();
        FILE.seekg(0, std::ios::beg);
        INDEX = 0;
        return next();
      }

      size_t index() { return INDEX; }
      std::string protein() { return PROTEIN; }
      std::string model() { return MODEL; }
      std::string chain() { return CHAIN; }
      std::string aminoacid() { return AMINOACID; }
    };

    // This class iterate through proteins and create an index for a knowledge base
    // NOTE: Due to efficiency reasons it is necessary to reset each stage before iterating it. 
    class Indexer {
      private:
      // File where to write index
      std::ofstream OUTPUT;
      // How to iterate proteins
      ProteinIterator* ITERATOR;

      public:
      // Create an index file and initialize iterator
      // filepath: path where to create index file
      // iterator: how to iterate proteins and what model(s), biomolecule(s) and crystallographic transformation(s) to use
      Indexer(std::string filepath, ProteinIterator* iterator) {
        if (filepath.empty() || filepath.back() == elemental::filesystem::directory_separator) {
          filepath += "residue.ind";
        } else if (!elemental::string::ends_with(filepath, ".ind")) {
          filepath += ".ind";
        }
        if (elemental::filesystem::exists(filepath)) {
          if (elemental::filesystem::is_directory(filepath)) {
            throw elemental::exception::TitledException("'" + filepath + "' exists and is a directory, please select another filename.");
          }
          std::cerr << "'" << filepath << "' exists and will be overriden." << std::endl;
        }
        OUTPUT.open(filepath);
        ITERATOR = iterator;
      }

      // Flush and close the file with created index
      ~Indexer() {
        OUTPUT.flush();
        OUTPUT.close();
      }

      // Adds 'protein' to the index file
      void index(Protein *protein) {
        if (ITERATOR->init(protein)) {
          if (ITERATOR->resetModel()) {
            bool first_model = true;
            do {
              if (ITERATOR->resetChain()) {
                bool first_chain = true;
                do {
                  if (ITERATOR->resetAminoacid()) {
                    bool first_aminoacid = true;
                    do {
                      OUTPUT << ITERATOR->getAminoacidName();
                      if (first_aminoacid) {
                        OUTPUT << '\t' << ITERATOR->getChainName();
                        first_aminoacid = false;
                        if (first_chain) {
                          OUTPUT << '\t' << ITERATOR->getModelName();
                          first_chain = false;
                          if (first_model) {
                            OUTPUT << '\t' << protein->ID_CODE;
                            first_model = false;
                          }
                        }
                      }
                      OUTPUT << std::endl;
                    } while (ITERATOR->nextAminoacid());
                  }
                } while (ITERATOR->nextChain());
              }
            } while (ITERATOR->nextModel());
          }
        }
      }

    };
  }
}