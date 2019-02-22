#pragma once

#include "../common/exception.h"
#include "protein.h"
#include <string>
#include <iostream>

namespace inspire {
  namespace backend {
    // Beware that iterators are lazy, i.e. each level must be reset, if a parent level is changed
    class ProteinIterator {
      public:
      virtual bool init(Protein* protein) = 0;

      virtual bool setModel(std::string name) = 0;
      virtual bool setChain(std::string name) = 0;
      virtual bool setAminoacid(std::string name) = 0;
      virtual bool setAtom(std::string name) = 0;

      virtual bool resetModel() { return true; }
      virtual bool resetChain() { return true; }
      virtual bool resetAminoacid() { return true; }
      virtual bool resetAtom() { return true; }
      
      virtual bool nextModel() { return false; };
      virtual bool nextChain() { return false; };
      virtual bool nextAminoacid() { return false; };
      virtual bool nextAtom() { return false; }
      
      virtual std::string getProteinName() = 0;
      virtual std::string getModelName() = 0;
      virtual std::string getChainName() = 0;
      virtual std::string getAminoacidName() = 0;
      virtual std::string getAtomName() = 0;
      // Default implementation expects chain name in format where chainID is the first followed by other informations and where chainID is a single character
      virtual std::string chainSymmetryGroup(const std::string &chain_name) {
        if (chain_name.empty()) {
          return "";
        }
        return chain_name.substr(0, 1);
      }

      virtual bool computeCharacteristics() = 0;

      // Residue type of the current residue
      virtual std::string residue_name() = 0;
      // Atom name
      virtual std::string atom_name() = 0;
      // Element type of the atom
      virtual std::string element() = 0;
      // Coordinates of the current atom, dependings on intern Characteristics selector or aggregator
      // NOTE: This function is preferred prior triplet x(), y() and z() due to repetitive intermediate calculations
      virtual Coordinate coordinates() = 0;
      // X coordinates of the current atom, dependings on intern Characteristics selector or aggregator
      virtual double x() = 0;
      // Y coordinates of the current atom, dependings on intern Characteristics selector or aggregator
      virtual double y() = 0;
      // Z coordinates of the current atom, dependings on intern Characteristics selector or aggregator
      virtual double z() = 0;
      // Occupancy of the current atom, dependings on intern Characteristics selector or aggregator
      virtual float occupancy() = 0;
      // Temperature of the current atom, dependings on intern Characteristics selector or aggregator
      virtual float temperature() = 0;
      // Charge of the current atom, dependings on intern Characteristics selector or aggregator
      virtual signed char charge() = 0;
    };

    // This iterator use only the first model, the first biomolecule, the first crystallographic transformation and the first alternative location
    // and iterates through biomolecule transformations, chains and aminoacids only
    class ExplicitIterator : public ProteinIterator {
      private:
      Protein* PROTEIN;
      std::map<std::string, Chain>::const_iterator CHAIN;
      std::map<std::pair<int, char>, Aminoacid>::const_iterator AMINOACID;
      std::map<std::string, Atom>::const_iterator ATOM;
      std::map<char, Characteristic>::const_iterator CHARACTERISTIC;

      // Set chain iterator on the first valid chain id
      // NOTE: Expects that there is at least one valid chain in each model
      bool set_first_chain() {
        CHAIN = PROTEIN->MODELS.begin()->second.CHAINS.begin();
        return CHAIN != PROTEIN->MODELS.begin()->second.CHAINS.end();
      }
      // Set chain iterator on the next valid chain id
      bool set_next_chain() {
        return ++CHAIN != PROTEIN->MODELS.begin()->second.CHAINS.end();
      }

      // NOTE: CHAIN is not validate as that is validated in the set_(first|next)_chain method and thus it will be redundant
      inline bool valid_chain() {
        return true;
      }
      inline bool valid_aminoacid() {
        return AMINOACID != CHAIN->second.AMINOACIDS.end();
      }
      inline bool valid_atom() {
        return ATOM != AMINOACID->second.ATOMS.end();
      }

      public:
      bool init(Protein* protein) override {
        if (protein->BIOMOLECULES.empty()) {
          throw common::exception::TitledException("Protein does contains no biomolecule definition");
        }
        if (protein->CRYSTALLOGRAPHIC_SYMMETRIES.empty()) {
          throw common::exception::TitledException("Protein does contains no crystallographic transformation");
        }
        PROTEIN = protein;
        return true;
      };

      bool resetChain() override {
        return set_first_chain() && valid_chain();
      }
      bool resetAminoacid() override {
        if (CHAIN == PROTEIN->MODELS.begin()->second.CHAINS.end()) {
          return false;
        }
        AMINOACID = CHAIN->second.AMINOACIDS.begin();
        // NOTE: Due to efficiency reasons, this test is used instead of valid_aminoacid() as it saves one lookup in std::map
        return AMINOACID != CHAIN->second.AMINOACIDS.end();
      }
      bool resetAtom() override {
        ATOM = AMINOACID->second.ATOMS.begin();
        return valid_atom();
      }

      bool nextChain() override {
        return set_next_chain() && valid_chain();
      };
      bool nextAminoacid() override {
        ++AMINOACID;
        return valid_aminoacid();
      };
      bool nextAtom() override {
        ++ATOM;
        return valid_atom();
      };

      bool setModel(std::string name) override {
        return resetModel();
      }
      bool setChain(std::string name) override {
        if (name.size() < 1 || !valid_chain()) {
          return false;
        }
        CHAIN = PROTEIN->MODELS.begin()->second.CHAINS.find(name);
        return CHAIN != PROTEIN->MODELS.begin()->second.CHAINS.end();
      }
      bool setAminoacid(std::string name) override {
        if (name.size() == 0) {
          return false;
        }
        int id;
        char alt_loc;
        if (name[name.size()-1] < '0' || name[name.size()-1] > '9') {
          id = std::stoi(name.substr(0, name.size()-1));
          alt_loc = name[name.size()-1];
        } else {
          id = std::stoi(name);
          alt_loc = ' ';
        }
        AMINOACID = CHAIN->second.AMINOACIDS.find({id, alt_loc});
        return valid_aminoacid();
      }
      bool setAtom(std::string name) override {
        ATOM = AMINOACID->second.ATOMS.find(name);
        return valid_atom();
      }

      // Returns protein's id_code
      std::string getProteinName() override {
        return PROTEIN->ID_CODE;
      };
      // Returns ""
      std::string getModelName() override {
        return "";
      };
      // Returns 'chainID'+'assemblyTransformationID'
      std::string getChainName() override {
        return CHAIN->first;
      }
      // Returns 'resSeq''iCode'
      std::string getAminoacidName() override {
        std::string ret = std::to_string(AMINOACID->first.first);
        if (AMINOACID->first.second != ' ') {
          ret.push_back(AMINOACID->first.second);
        }
        return ret;
      };
      // Returns 'name'
      std::string getAtomName() override {
        return ATOM->first;
      };

      // Set Characteristics with the highest occupancy.
      // If multiple Characteristics have the same occupancy, the first one is used.
      bool computeCharacteristics() override {
        CHARACTERISTIC = ATOM->second.ALTERNATIVE_LOCATIONS.begin();
        if (CHARACTERISTIC == ATOM->second.ALTERNATIVE_LOCATIONS.end()) {
          return false;
        }
        for (auto it = ATOM->second.ALTERNATIVE_LOCATIONS.begin(); ++it != ATOM->second.ALTERNATIVE_LOCATIONS.end(); ) {
          if (CHARACTERISTIC->second.OCCUPANCY < it->second.OCCUPANCY) {
            CHARACTERISTIC = it;
          }
        }
        return true;
      }

      std::string residue_name() override {
        return AMINOACID->second.RESIDUE_NAME;
      }
      std::string atom_name() override {
        return ATOM->first;
      }
      std::string element() override {
        return ATOM->second.ELEMENT;
      }
      Coordinate coordinates() override {
        return std::make_tuple(CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
      }
      double x() override {
        return CHARACTERISTIC->second.X;
      }
      double y() override {
        return CHARACTERISTIC->second.Y;
      }
      double z() override {
        return CHARACTERISTIC->second.Z;
      }
      float occupancy() override {
        return CHARACTERISTIC->second.OCCUPANCY;
      }
      float temperature() override {
        return CHARACTERISTIC->second.TEMPERATURE;
      }
      signed char charge() override {
        return CHARACTERISTIC->second.CHARGE;
      }
    };

    // This iterator use only the first model, the first biomolecule, the first crystallographic transformation and the first alternative location
    // and iterates through biomolecule transformations, chains and aminoacids only
    class FirstModelIterator : public ProteinIterator {
      private:
      Protein* PROTEIN;
      std::map<std::string, std::pair<TransformationMatrix, std::set<std::string>>>::const_iterator BIOMOLECULE_TRANSFORMATION;
      std::set<std::string>::const_iterator CHAIN;
      std::map<std::pair<int, char>, Aminoacid>::const_iterator AMINOACID;
      std::map<std::string, Atom>::const_iterator ATOM;
      std::map<char, Characteristic>::const_iterator CHARACTERISTIC;

      // Set chain iterator on the first valid chain id
      // NOTE: Expects that there is at least one valid chain in each model
      bool set_first_chain() {
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.begin();
        while (CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN)) {
            return true;
          }
          ++CHAIN;
        }
        return false;
      }
      // Set chain iterator on the next valid chain id
      bool set_next_chain() {
        while (++CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN)) {
            return true;
          }
        }
        return false;
      }

      // NOTE: CHAIN is not validate as that is validated in the set_(first|next)_chain method and thus it will be redundant
      inline bool valid_chain() {
        return BIOMOLECULE_TRANSFORMATION != PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end();
      }
      inline bool valid_aminoacid() {
        return AMINOACID != PROTEIN->MODELS.begin()->second.CHAINS.at(*CHAIN).AMINOACIDS.end();
      }
      inline bool valid_atom() {
        return ATOM != AMINOACID->second.ATOMS.end();
      }

      public:
      bool init(Protein* protein) override {
        if (protein->BIOMOLECULES.empty()) {
          throw common::exception::TitledException("Protein does contains no biomolecule definition");
        }
        if (protein->CRYSTALLOGRAPHIC_SYMMETRIES.empty()) {
          throw common::exception::TitledException("Protein does contains no crystallographic transformation");
        }
        PROTEIN = protein;
        return true;
      };

      bool resetChain() override {
        BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.begin();
        while (!set_first_chain()) {
          if (++BIOMOLECULE_TRANSFORMATION == PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end()) {
            return false;
          }
        }
        return valid_chain();
      }
      bool resetAminoacid() override {
        auto chain = PROTEIN->MODELS.begin()->second.CHAINS.find(*CHAIN);
        if (chain == PROTEIN->MODELS.begin()->second.CHAINS.end()) {
          return false;
        }
        AMINOACID = chain->second.AMINOACIDS.begin();
        // NOTE: Due to efficiency reasons, this test is used instead of valid_aminoacid() as it saves one lookup in std::map
        return AMINOACID != chain->second.AMINOACIDS.end();
      }
      bool resetAtom() override {
        ATOM = AMINOACID->second.ATOMS.begin();
        return valid_atom();
      }

      bool nextChain() override {
        if (!set_next_chain()) {
          do {
            ++BIOMOLECULE_TRANSFORMATION;
            if (BIOMOLECULE_TRANSFORMATION == PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end()) {
              return false;
            }
          } while (!set_first_chain());
        }
        return valid_chain();
      };
      bool nextAminoacid() override {
        ++AMINOACID;
        return valid_aminoacid();
      };
      bool nextAtom() override {
        ++ATOM;
        return valid_atom();
      };

      bool setModel(std::string name) override {
        return resetModel();
      }
      bool setChain(std::string name) override {
        size_t plus = name.find('+');
        if (plus == name.npos) {
          BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.begin();
        } else {
          BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.find(name.substr(plus+1));
        }
        if (!valid_chain()) {
          return false;
        }
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.find(name.substr(0, plus));
        return CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end() && PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN);
      }
      bool setAminoacid(std::string name) override {
        if (name.size() == 0) {
          return false;
        }
        int id;
        char alt_loc;
        if (name[name.size()-1] < '0' || name[name.size()-1] > '9') {
          id = std::stoi(name.substr(0, name.size()-1));
          alt_loc = name[name.size()-1];
        } else {
          id = std::stoi(name);
          alt_loc = ' ';
        }
        AMINOACID = PROTEIN->MODELS.begin()->second.CHAINS.at(*CHAIN).AMINOACIDS.find({id, alt_loc});
        return valid_aminoacid();
      }
      bool setAtom(std::string name) override {
        ATOM = AMINOACID->second.ATOMS.find(name);
        return valid_atom();
      }

      // Returns protein's id_code
      std::string getProteinName() override {
        return PROTEIN->ID_CODE;
      };
      // Returns ""
      std::string getModelName() override {
        return "";
      };
      // Returns 'chainID'+'assemblyTransformationID'
      std::string getChainName() override {
        std::string ret(*CHAIN);
        if (BIOMOLECULE_TRANSFORMATION->first != "1") {
          ret.push_back('+');
          ret.append(BIOMOLECULE_TRANSFORMATION->first);
        }
        return ret;
      }
      // Returns 'resSeq''iCode'
      std::string getAminoacidName() override {
        std::string ret = std::to_string(AMINOACID->first.first);
        if (AMINOACID->first.second != ' ') {
          ret.push_back(AMINOACID->first.second);
        }
        return ret;
      };
      // Returns 'name'
      std::string getAtomName() override {
        return ATOM->first;
      };

      // Set Characteristics with the highest occupancy.
      // If multiple Characteristics have the same occupancy, the first one is used.
      bool computeCharacteristics() override {
        CHARACTERISTIC = ATOM->second.ALTERNATIVE_LOCATIONS.begin();
        if (CHARACTERISTIC == ATOM->second.ALTERNATIVE_LOCATIONS.end()) {
          return false;
        }
        for (auto it = ATOM->second.ALTERNATIVE_LOCATIONS.begin(); ++it != ATOM->second.ALTERNATIVE_LOCATIONS.end(); ) {
          if (CHARACTERISTIC->second.OCCUPANCY < it->second.OCCUPANCY) {
            CHARACTERISTIC = it;
          }
        }
        return true;
      }

      std::string residue_name() override {
        return AMINOACID->second.RESIDUE_NAME;
      }
      std::string atom_name() override {
        return ATOM->first;
      }
      std::string element() override {
        return ATOM->second.ELEMENT;
      }
      Coordinate coordinates() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return std::make_tuple(PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z));
      }
      double x() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      double y() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      double z() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      float occupancy() override {
        return CHARACTERISTIC->second.OCCUPANCY;
      }
      float temperature() override {
        return CHARACTERISTIC->second.TEMPERATURE;
      }
      signed char charge() override {
        return CHARACTERISTIC->second.CHARGE;
      }
    };

    // This iterator use only the first model, the first biomolecule and the first alternative location
    // and iterates through biomolecule transformations, crystallographic transformations, chains and aminoacids only
    class FirstModelCrystallographicIterator : public ProteinIterator {
      private:
      Protein* PROTEIN;
      std::map<int, TransformationMatrix>::const_iterator CRYSTALLOGRAPHIC_TRANSFORMATION;
      std::map<std::string, std::pair<TransformationMatrix, std::set<std::string>>>::const_iterator BIOMOLECULE_TRANSFORMATION;
      std::set<std::string>::const_iterator CHAIN;
      std::map<std::pair<int, char>, Aminoacid>::const_iterator AMINOACID;
      std::map<std::string, Atom>::const_iterator ATOM;
      std::map<char, Characteristic>::const_iterator CHARACTERISTIC;

      // Set chain iterator on the first valid chain id
      // NOTE: Expects that there is at least one valid chain in each model
      bool set_first_chain() {
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.begin();
        while (CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN)) {
            return true;
          }
          ++CHAIN;
        }
        return false;
      }
      // Set chain iterator on the next valid chain id
      bool set_next_chain() {
        while (++CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN)) {
            return true;
          }
        }
        return false;
      }

      // NOTE: CHAIN is not validate as that is validated in the set_(first|next)_chain method and thus it will be redundant
      inline bool valid_chain() {
        return CRYSTALLOGRAPHIC_TRANSFORMATION != PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.end()
            && BIOMOLECULE_TRANSFORMATION != PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end();
      }
      inline bool valid_aminoacid() {
        return AMINOACID != PROTEIN->MODELS.begin()->second.CHAINS.at(*CHAIN).AMINOACIDS.end();
      }
      inline bool valid_atom() {
        return ATOM != AMINOACID->second.ATOMS.end();
      }

      public:
      bool init(Protein* protein) override {
        if (protein->BIOMOLECULES.empty()) {
          throw common::exception::TitledException("Protein does contains no biomolecule definition");
        }
        PROTEIN = protein;
        return true;
      };

      bool resetChain() override {
        CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
        BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.begin();
        while (!set_first_chain()) {
          if (++BIOMOLECULE_TRANSFORMATION == PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end()) {
            return false;
          }
        }
        return valid_chain();
      }
      bool resetAminoacid() override {
        auto chain = PROTEIN->MODELS.begin()->second.CHAINS.find(*CHAIN);
        if (chain == PROTEIN->MODELS.begin()->second.CHAINS.end()) {
          return false;
        }
        AMINOACID = chain->second.AMINOACIDS.begin();
        // NOTE: Due to efficiency reasons, this test is used instead of valid_aminoacid() as it saves one lookup in std::map
        return AMINOACID != chain->second.AMINOACIDS.end();
      }
      bool resetAtom() override {
        ATOM = AMINOACID->second.ATOMS.begin();
        return valid_atom();
      }

      bool nextChain() override {
        ++CRYSTALLOGRAPHIC_TRANSFORMATION;
        if (CRYSTALLOGRAPHIC_TRANSFORMATION == PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.end()) {
          if (!set_next_chain()) {
            do {
              ++BIOMOLECULE_TRANSFORMATION;
              if (BIOMOLECULE_TRANSFORMATION == PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.end()) {
                return false;
              }
            } while (!set_first_chain());
          }
          CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
        }
        return valid_chain();
      };
      bool nextAminoacid() override {
        ++AMINOACID;
        return valid_aminoacid();
      };
      bool nextAtom() override {
        ++ATOM;
        return valid_atom();
      };

      bool setModel(std::string name) override {
        return resetModel();
      }
      bool setChain(std::string name) override {
        size_t plus = name.find('+');
        size_t star = name.find('*');
        std::string chain;
        if (plus == name.npos) {
          if (star == name.npos) {
            BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.begin();
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
            chain = name;
          } else {
            BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.begin();
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1)));
            chain = name.substr(0, star);
          }
        } else {
          if (star == name.npos) {
            BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.find(name.substr(plus+1));
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
            chain = name.substr(0, plus);
          } else {
            if (plus < star) {
              BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.find(name.substr(plus+1, star-plus-1));
              CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1)));
              chain = name.substr(0, plus);
            } else {
              BIOMOLECULE_TRANSFORMATION = PROTEIN->BIOMOLECULES.begin()->second.TRANSFORMATIONS.find(name.substr(plus+1));
              CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1, plus-star-1)));
              chain = name.substr(0, star);
            }
          }
        }
        if (!valid_chain()) {
          return false;
        }
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.find(chain);
        return CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end() && PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN);
      }
      bool setAminoacid(std::string name) override {
        if (name.size() == 0) {
          return false;
        }
        int id;
        char alt_loc;
        if (name[name.size()-1] < '0' || name[name.size()-1] > '9') {
          id = std::stoi(name.substr(0, name.size()-1));
          alt_loc = name[name.size()-1];
        } else {
          id = std::stoi(name);
          alt_loc = ' ';
        }
        AMINOACID = PROTEIN->MODELS.begin()->second.CHAINS.at(*CHAIN).AMINOACIDS.find({id, alt_loc});
        return valid_aminoacid();
      }
      bool setAtom(std::string name) override {
        ATOM = AMINOACID->second.ATOMS.find(name);
        return valid_atom();
      }

      // Returns protein's id_code
      std::string getProteinName() override {
        return PROTEIN->ID_CODE;
      };
      // Returns ""
      std::string getModelName() override {
        return "";
      };
      // Returns 'chainID'+'assemblyTransformationID'*'crystallographicTransformationID'
      std::string getChainName() override {
        std::string ret(*CHAIN);
        if (BIOMOLECULE_TRANSFORMATION->first != "1") {
          ret.push_back('+');
          ret.append(BIOMOLECULE_TRANSFORMATION->first);
        }
        if (CRYSTALLOGRAPHIC_TRANSFORMATION->first > 1) {
          ret.push_back('*');
          ret.append(std::to_string(CRYSTALLOGRAPHIC_TRANSFORMATION->first));
        }
        return ret;
      }
      // Returns 'resSeq''iCode'
      std::string getAminoacidName() override {
        std::string ret = std::to_string(AMINOACID->first.first);
        if (AMINOACID->first.second != ' ') {
          ret.push_back(AMINOACID->first.second);
        }
        return ret;
      };
      // Returns 'name'
      std::string getAtomName() override {
        return ATOM->first;
      };

      // Set Characteristics with the highest occupancy.
      // If multiple Characteristics have the same occupancy, the first one is used.
      bool computeCharacteristics() override {
        CHARACTERISTIC = ATOM->second.ALTERNATIVE_LOCATIONS.begin();
        if (CHARACTERISTIC == ATOM->second.ALTERNATIVE_LOCATIONS.end()) {
          return false;
        }
        for (auto it = ATOM->second.ALTERNATIVE_LOCATIONS.begin(); ++it != ATOM->second.ALTERNATIVE_LOCATIONS.end(); ) {
          if (CHARACTERISTIC->second.OCCUPANCY < it->second.OCCUPANCY) {
            CHARACTERISTIC = it;
          }
        }
        return true;
      }

      std::string residue_name() override {
        return AMINOACID->second.RESIDUE_NAME;
      }
      std::string atom_name() override {
        return ATOM->first;
      }
      std::string element() override {
        return ATOM->second.ELEMENT;
      }
      Coordinate coordinates() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return std::make_tuple(PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z));
      }
      double x() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<0>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      double y() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<1>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      double z() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<2>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      float occupancy() override {
        return CHARACTERISTIC->second.OCCUPANCY;
      }
      float temperature() override {
        return CHARACTERISTIC->second.TEMPERATURE;
      }
      signed char charge() override {
        return CHARACTERISTIC->second.CHARGE;
      }
    };

    // This iterator use only the first crystallographic transformations and the first alternative location
    // and iterates through models, biomolecules biomolecule transformations, chains and aminoacids only
    class BiomoleculesIterator : public ProteinIterator {
      private:
      Protein* PROTEIN;
      std::map<int, Model>::const_iterator MODEL;
      std::map<int, Biomolecule>::const_iterator BIOMOLECULE;
      std::map<std::string, std::pair<TransformationMatrix, std::set<std::string>>>::const_iterator BIOMOLECULE_TRANSFORMATION;
      std::set<std::string>::const_iterator CHAIN;
      std::map<std::pair<int, char>, Aminoacid>::const_iterator AMINOACID;
      std::map<std::string, Atom>::const_iterator ATOM;
      std::map<char, Characteristic>::const_iterator CHARACTERISTIC;

      // Set chain iterator on the first valid chain id
      // NOTE: Expects that there is at least one valid chain in each model
      bool set_first_chain() {
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.begin();
        while (CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (MODEL->second.CHAINS.count(*CHAIN)) {
            return true;
          }
          ++CHAIN;
        }
        return false;
      }
      // Set chain iterator on the next valid chain id
      bool set_next_chain() {
        while (++CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (MODEL->second.CHAINS.count(*CHAIN)) {
            return true;
          }
        }
        return false;
      }

      inline bool valid_model() {
        return MODEL != PROTEIN->MODELS.end() && BIOMOLECULE != PROTEIN->BIOMOLECULES.end();
      }
      // NOTE: CHAIN is not validate as that is validated in the set_(first|next)_chain method and thus it will be redundant
      inline bool valid_chain() {
        return BIOMOLECULE_TRANSFORMATION != BIOMOLECULE->second.TRANSFORMATIONS.end();
      }
      inline bool valid_aminoacid() {
        return AMINOACID != MODEL->second.CHAINS.at(*CHAIN).AMINOACIDS.end();
      }
      inline bool valid_atom() {
        return ATOM != AMINOACID->second.ATOMS.end();
      }

      public:
      bool init(Protein* protein) override {
        if (protein->CRYSTALLOGRAPHIC_SYMMETRIES.empty()) {
          throw common::exception::TitledException("Protein does contains no crystallographic transformation");
        }
        PROTEIN = protein;
        return true;
      };

      bool resetModel() override {
        MODEL = PROTEIN->MODELS.begin();
        BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
        return valid_model();
      }
      bool resetChain() override {
        BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.begin();
        while (!set_first_chain()) {
          if (++BIOMOLECULE_TRANSFORMATION == BIOMOLECULE->second.TRANSFORMATIONS.end()) {
            return false;
          }
        }
        return valid_chain();
      }
      bool resetAminoacid() override {
        auto chain = MODEL->second.CHAINS.find(*CHAIN);
        if (chain == MODEL->second.CHAINS.end()) {
          return false;
        }
        AMINOACID = chain->second.AMINOACIDS.begin();
        // NOTE: Due to efficiency reasons, this test is used instead of valid_aminoacid() as it saves one lookup in std::map
        return AMINOACID != chain->second.AMINOACIDS.end();
      }
      bool resetAtom() override {
        ATOM = AMINOACID->second.ATOMS.begin();
        return valid_atom();
      }

      bool nextModel() override {
        ++BIOMOLECULE;
        if (BIOMOLECULE == PROTEIN->BIOMOLECULES.end()) {
          ++MODEL;
          BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
        }
        return valid_model();
      };
      bool nextChain() override {
        if (!set_next_chain()) {
          do {
            ++BIOMOLECULE_TRANSFORMATION;
            if (BIOMOLECULE_TRANSFORMATION == BIOMOLECULE->second.TRANSFORMATIONS.end()) {
              return false;
            }
          } while (!set_first_chain());
        }
        return valid_chain();
      };
      bool nextAminoacid() override {
        ++AMINOACID;
        return valid_aminoacid();
      };
      bool nextAtom() override {
        ++ATOM;
        return valid_atom();
      };

      bool setModel(std::string name) override {
        if (name.empty()) {
          return resetModel();
        }
        if (name[0] == ';') {
          size_t comma = name.find(',');
          if (comma == name.npos) {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(1)));
            MODEL = PROTEIN->MODELS.begin();
          } else {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(1, comma-1)));
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(comma+1)));
          }
        } else if (name[0] == ',') {
          size_t semicolon = name.find(';');
          if (semicolon == name.npos) {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(1)));
          } else {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(semicolon+1)));
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(1, semicolon-1)));
          }
        } else {
          throw common::exception::TitledException("Invalid format of model code, 'AllExceptAltLocIterator' does not recognize first character of '" + name + "'");
        }
        return valid_model();
      }
      bool setChain(std::string name) override {
        size_t plus = name.find('+');
        if (plus == name.npos) {
          BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.begin();
        } else {
          BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.find(name.substr(plus+1));
        }
        if (!valid_chain()) {
          return false;
        }
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.find(name.substr(0, plus));
        return CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end() && PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN);
      }
      bool setAminoacid(std::string name) override {
        if (name.size() == 0) {
          return false;
        }
        int id;
        char alt_loc;
        if (name[name.size()-1] < '0' || name[name.size()-1] > '9') {
          id = std::stoi(name.substr(0, name.size()-1));
          alt_loc = name[name.size()-1];
        } else {
          id = std::stoi(name);
          alt_loc = ' ';
        }
        AMINOACID = MODEL->second.CHAINS.at(*CHAIN).AMINOACIDS.find({id, alt_loc});
        return valid_aminoacid();
      }
      bool setAtom(std::string name) override {
        ATOM = AMINOACID->second.ATOMS.find(name);
        return valid_atom();
      }

      // Returns protein's id_code
      std::string getProteinName() override {
        return PROTEIN->ID_CODE;
      };
      // Returns ;'biomoleculeID','modelID'
      std::string getModelName() override {
        std::string ret = "";
        if (BIOMOLECULE->first > 1) {
          ret.push_back(';');
          ret.append(std::to_string(BIOMOLECULE->first));
        }
        if (MODEL->first > 1) {
          ret.push_back(',');
          ret.append(std::to_string(MODEL->first));
        }
        return ret;
      };
      // Returns 'chainID'+'assemblyTransformationID'
      std::string getChainName() override {
        std::string ret(*CHAIN);
        if (BIOMOLECULE_TRANSFORMATION->first != "1") {
          ret.push_back('+');
          ret.append(BIOMOLECULE_TRANSFORMATION->first);
        }
        return ret;
      }
      // Returns 'resSeq''iCode'
      std::string getAminoacidName() override {
        std::string ret = std::to_string(AMINOACID->first.first);
        if (AMINOACID->first.second != ' ') {
          ret.push_back(AMINOACID->first.second);
        }
        return ret;
      };
      // Returns 'name'
      std::string getAtomName() override {
        return ATOM->first;
      };

      // Set Characteristics with the highest occupancy.
      // If multiple Characteristics have the same occupancy, the first one is used.
      bool computeCharacteristics() override {
        CHARACTERISTIC = ATOM->second.ALTERNATIVE_LOCATIONS.begin();
        if (CHARACTERISTIC == ATOM->second.ALTERNATIVE_LOCATIONS.end()) {
          return false;
        }
        for (auto it = ATOM->second.ALTERNATIVE_LOCATIONS.begin(); ++it != ATOM->second.ALTERNATIVE_LOCATIONS.end(); ) {
          if (CHARACTERISTIC->second.OCCUPANCY < it->second.OCCUPANCY) {
            CHARACTERISTIC = it;
          }
        }
        return true;
      }

      std::string residue_name() override {
        return AMINOACID->second.RESIDUE_NAME;
      }
      std::string atom_name() override {
        return ATOM->first;
      }
      std::string element() override {
        return ATOM->second.ELEMENT;
      }
      Coordinate coordinates() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return std::make_tuple(PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z));
      }
      double x() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      double y() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      double z() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z);
      }
      float occupancy() override {
        return CHARACTERISTIC->second.OCCUPANCY;
      }
      float temperature() override {
        return CHARACTERISTIC->second.TEMPERATURE;
      }
      signed char charge() override {
        return CHARACTERISTIC->second.CHARGE;
      }
    };

    // This iterator use only the first alternative location
    // and iterates through models, biomolecules biomolecule transformations, crystallographic transformations, chains and aminoacids only
    class AllExceptAltLocIterator : public ProteinIterator {
      private:
      Protein* PROTEIN;
      std::map<int, Model>::const_iterator MODEL;
      std::map<int, Biomolecule>::const_iterator BIOMOLECULE;
      std::map<int, TransformationMatrix>::const_iterator CRYSTALLOGRAPHIC_TRANSFORMATION;
      std::map<std::string, std::pair<TransformationMatrix, std::set<std::string>>>::const_iterator BIOMOLECULE_TRANSFORMATION;
      std::set<std::string>::const_iterator CHAIN;
      std::map<std::pair<int, char>, Aminoacid>::const_iterator AMINOACID;
      std::map<std::string, Atom>::const_iterator ATOM;
      std::map<char, Characteristic>::const_iterator CHARACTERISTIC;

      // Set chain iterator on the first valid chain id
      // NOTE: Expects that there is at least one valid chain in each model
      bool set_first_chain() {
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.begin();
        while (CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (MODEL->second.CHAINS.count(*CHAIN)) {
            return true;
          }
          ++CHAIN;
        }
        return false;
      }
      // Set chain iterator on the next valid chain id
      bool set_next_chain() {
        while (++CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end()) {
          if (MODEL->second.CHAINS.count(*CHAIN)) {
            return true;
          }
        }
        return false;
      }

      inline bool valid_model() {
        return MODEL != PROTEIN->MODELS.end() && BIOMOLECULE != PROTEIN->BIOMOLECULES.end();
      }
      // NOTE: CHAIN is not validate as that is validated in the set_(first|next)_chain method and thus it will be redundant
      inline bool valid_chain() {
        return CRYSTALLOGRAPHIC_TRANSFORMATION != PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.end() && BIOMOLECULE_TRANSFORMATION != BIOMOLECULE->second.TRANSFORMATIONS.end();
      }
      inline bool valid_aminoacid() {
        return AMINOACID != MODEL->second.CHAINS.at(*CHAIN).AMINOACIDS.end();
      }
      inline bool valid_atom() {
        return ATOM != AMINOACID->second.ATOMS.end();
      }

      public:
      bool init(Protein* protein) override {
        PROTEIN = protein;
        return true;
      };

      bool resetModel() override {
        MODEL = PROTEIN->MODELS.begin();
        BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
        return valid_model();
      }
      bool resetChain() override {
        CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
        BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.begin();
        while (!set_first_chain()) {
          if (++BIOMOLECULE_TRANSFORMATION == BIOMOLECULE->second.TRANSFORMATIONS.end()) {
            return false;
          }
        }
        return valid_chain();
      }
      bool resetAminoacid() override {
        auto chain = MODEL->second.CHAINS.find(*CHAIN);
        if (chain == MODEL->second.CHAINS.end()) {
          return false;
        }
        AMINOACID = chain->second.AMINOACIDS.begin();
        // NOTE: Due to efficiency reasons, this test is used instead of valid_aminoacid() as it saves one lookup in std::map
        return AMINOACID != chain->second.AMINOACIDS.end();
      }
      bool resetAtom() override {
        ATOM = AMINOACID->second.ATOMS.begin();
        return valid_atom();
      }

      bool nextModel() override {
        ++BIOMOLECULE;
        if (BIOMOLECULE == PROTEIN->BIOMOLECULES.end()) {
          ++MODEL;
          BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
        }
        return valid_model();
      };
      bool nextChain() override {
        ++CRYSTALLOGRAPHIC_TRANSFORMATION;
        if (CRYSTALLOGRAPHIC_TRANSFORMATION == PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.end()) {
          if (!set_next_chain()) {
            do {
              ++BIOMOLECULE_TRANSFORMATION;
              if (BIOMOLECULE_TRANSFORMATION == BIOMOLECULE->second.TRANSFORMATIONS.end()) {
                return false;
              }
            } while (!set_first_chain());
          }
          CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
        }
        return valid_chain();
      };
      bool nextAminoacid() override {
        ++AMINOACID;
        return valid_aminoacid();
      };
      bool nextAtom() override {
        ++ATOM;
        return valid_atom();
      };

      bool setModel(std::string name) override {
        if (name.empty()) {
          return resetModel();
        }
        if (name[0] == ';') {
          size_t comma = name.find(',');
          if (comma == name.npos) {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(1)));
            MODEL = PROTEIN->MODELS.begin();
          } else {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(1, comma-1)));
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(comma+1)));
          }
        } else if (name[0] == ',') {
          size_t semicolon = name.find(';');
          if (semicolon == name.npos) {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.begin();
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(1)));
          } else {
            BIOMOLECULE = PROTEIN->BIOMOLECULES.find(std::stoi(name.substr(semicolon+1)));
            MODEL = PROTEIN->MODELS.find(std::stoi(name.substr(1, semicolon-1)));
          }
        } else {
          throw common::exception::TitledException("Invalid format of model code, 'AllExceptAltLocIterator' does not recognize first character of '" + name + "'");
        }
        return valid_model();
      }
      bool setChain(std::string name) override {
        size_t plus = name.find('+');
        size_t star = name.find('*');
        std::string chain;
        if (plus == name.npos) {
          if (star == name.npos) {
            BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.begin();
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
            chain = name;
          } else {
            BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.begin();
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1)));
            chain = name.substr(0, star);
          }
        } else {
          if (star == name.npos) {
            BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.find(name.substr(plus+1));
            CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin();
            chain = name.substr(0, plus);
          } else {
            if (plus < star) {
              BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.find(name.substr(plus+1, star-plus-1));
              CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1)));
              chain = name.substr(0, plus);
            } else {
              BIOMOLECULE_TRANSFORMATION = BIOMOLECULE->second.TRANSFORMATIONS.find(name.substr(plus+1));
              CRYSTALLOGRAPHIC_TRANSFORMATION = PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.find(std::stoi(name.substr(star+1, plus-star-1)));
              chain = name.substr(0, star);
            }
          }
        }
        if (!valid_chain()) {
          return false;
        }
        CHAIN = BIOMOLECULE_TRANSFORMATION->second.second.find(chain);
        return CHAIN != BIOMOLECULE_TRANSFORMATION->second.second.end() && PROTEIN->MODELS.begin()->second.CHAINS.count(*CHAIN);
      }
      bool setAminoacid(std::string name) override {
        if (name.size() == 0) {
          return false;
        }
        int id;
        char alt_loc;
        if (name[name.size()-1] < '0' || name[name.size()-1] > '9') {
          id = std::stoi(name.substr(0, name.size()-1));
          alt_loc = name[name.size()-1];
        } else {
          id = std::stoi(name);
          alt_loc = ' ';
        }
        AMINOACID = MODEL->second.CHAINS.at(*CHAIN).AMINOACIDS.find({id, alt_loc});
        return valid_aminoacid();
      }
      bool setAtom(std::string name) override {
        ATOM = AMINOACID->second.ATOMS.find(name);
        return valid_atom();
      }

      // Returns protein's id_code
      std::string getProteinName() override {
        return PROTEIN->ID_CODE;
      };
      // Returns ;'biomoleculeID','modelID'
      std::string getModelName() override {
        std::string ret = "";
        if (BIOMOLECULE->first > 1) {
          ret.push_back(';');
          ret.append(std::to_string(BIOMOLECULE->first));
        }
        if (MODEL->first > 1) {
          ret.push_back(',');
          ret.append(std::to_string(MODEL->first));
        }
        return ret;
      };
      // Returns 'chainID'+'assemblyTransformationID'*'crystallographicTransformationID'
      std::string getChainName() override {
        std::string ret(*CHAIN);
        if (BIOMOLECULE_TRANSFORMATION->first != "1") {
          ret.push_back('+');
          ret.append(BIOMOLECULE_TRANSFORMATION->first);
        }
        if (CRYSTALLOGRAPHIC_TRANSFORMATION->first > 1) {
          ret.push_back('*');
          ret.append(std::to_string(CRYSTALLOGRAPHIC_TRANSFORMATION->first));
        }
        return ret;
      }
      // Returns 'resSeq''iCode'
      std::string getAminoacidName() override {
        std::string ret = std::to_string(AMINOACID->first.first);
        if (AMINOACID->first.second != ' ') {
          ret.push_back(AMINOACID->first.second);
        }
        return ret;
      };
      // Returns 'name'
      std::string getAtomName() override {
        return ATOM->first;
      };

      // Set Characteristics with the highest occupancy.
      // If multiple Characteristics have the same occupancy, the first one is used.
      bool computeCharacteristics() override {
        CHARACTERISTIC = ATOM->second.ALTERNATIVE_LOCATIONS.begin();
        if (CHARACTERISTIC == ATOM->second.ALTERNATIVE_LOCATIONS.end()) {
          return false;
        }
        for (auto it = ATOM->second.ALTERNATIVE_LOCATIONS.begin(); ++it != ATOM->second.ALTERNATIVE_LOCATIONS.end(); ) {
          if (CHARACTERISTIC->second.OCCUPANCY < it->second.OCCUPANCY) {
            CHARACTERISTIC = it;
          }
        }
        return true;
      }

      std::string residue_name() override {
        return AMINOACID->second.RESIDUE_NAME;
      }
      std::string atom_name() override {
        return ATOM->first;
      }
      std::string element() override {
        return ATOM->second.ELEMENT;
      }
      Coordinate coordinates() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return std::make_tuple(PROTEIN->transform(std::get<0>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<1>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z),
                               PROTEIN->transform(std::get<2>(PROTEIN->CRYSTALLOGRAPHIC_SYMMETRIES.begin()->second), x, y, z));
      }
      double x() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<0>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      double y() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<1>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      double z() override {
        double x = PROTEIN->transform(std::get<0>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double y = PROTEIN->transform(std::get<1>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        double z = PROTEIN->transform(std::get<2>(BIOMOLECULE_TRANSFORMATION->second.first), CHARACTERISTIC->second.X, CHARACTERISTIC->second.Y, CHARACTERISTIC->second.Z);
        return PROTEIN->transform(std::get<2>(CRYSTALLOGRAPHIC_TRANSFORMATION->second), x, y, z);
      }
      float occupancy() override {
        return CHARACTERISTIC->second.OCCUPANCY;
      }
      float temperature() override {
        return CHARACTERISTIC->second.TEMPERATURE;
      }
      signed char charge() override {
        return CHARACTERISTIC->second.CHARGE;
      }
    };
  }
}