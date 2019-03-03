#pragma once

#include "../common/filesystem.h"
#include "index.h"
#include "octree.h"
#include "parser.h"
#include "protein.h"
#include <map>
#include <cmath>
#include <iostream>
#include <limits>

namespace inspire {
  namespace backend {
    // Extract features from protein
    class IFeature {
      public:
      // Returns th name of this feature
      virtual std::string title() = 0;

      // Inform about the current protein for the case that something joint needs to be precomputed
      virtual void init(Protein* protein) { }
    };

    template <class T> class Feature : public IFeature {
      protected:
      // How to select required aminoacids from proteins
      ProteinIterator* ITERATOR;

      public:
      static const T UNDEFINED;

      // Initializes Feature extractor with how to select required aminoacids and their features from proteins
      Feature(ProteinIterator* iterator) : ITERATOR(iterator) { }

      // Returns the feature of the given aminoacid within the given chain within the given model from the given protein
      virtual T feature(const std::string &model, const std::string &chain, const std::string &aminoacid) = 0;
    };
    template <> const size_t Feature<size_t>::UNDEFINED = -1;
    template <> const float Feature<float>::UNDEFINED = std::numeric_limits<float>::infinity();
    template <> const std::string Feature<std::string>::UNDEFINED = "";

    template <class T> class ToStringFeature : public Feature<std::string> {
      private:
      Feature<T>* FEATURE;

      public:
      // *.bod
      ToStringFeature(ProteinIterator* iterator, Feature<T>* feature) : Feature(iterator), FEATURE(feature) { }

      std::string title() override { return FEATURE->title(); }

      void init(Protein* protein) override {
        FEATURE->init(protein);
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        T ret = FEATURE->feature(model, chain, aminoacid);
        return ret == FEATURE->UNDEFINED ? UNDEFINED : std::to_string(ret);
      }
    };

    class BinFeature : Feature<size_t> {
      private:
      Feature<float>* FEATURE;
      std::map<float, size_t> BINS;

      public:
      // *.bod
      BinFeature(ProteinIterator* iterator, Feature<float>* feature, std::string boundaries) : Feature(iterator), FEATURE(feature) {
        std::ifstream input(boundaries);
        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          BINS[std::stof(line)] = BINS.size();
        }
        input.close();
      }

      std::string title() override { return FEATURE->title() + std::to_string(BINS.size()+1); }

      void init(Protein* protein) override {
        FEATURE->init(protein);
      }

      size_t feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        float ret = FEATURE->feature(model, chain, aminoacid);
        if (ret == FEATURE->UNDEFINED) {
          return UNDEFINED;
        } else {
          auto it = BINS.upper_bound(ret);
          return it == BINS.end() ? BINS.size() : it->second;
        }
      }
    };

    class RelativeAminoacidFeature : public Feature<float> {
      private:
      Feature<float>* FEATURE;
      std::map<std::string, float> MAXIMAL;

      public:
      RelativeAminoacidFeature(ProteinIterator* iterator, Feature<float>* feature, std::string maximals) : Feature(iterator), FEATURE(feature) {
        std::ifstream input(maximals);
        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          size_t tab = line.find('\t');
          if (tab == line.npos) {
            std::cerr << "File with van der Waals radiuses '" << maximals << "' contains invalid line '" << line << "'" << std::endl;
          } else {
            auto ins = MAXIMAL.insert({line.substr(0, tab), std::stof(line.substr(tab+1))});
            if (!ins.second) {
              std::cerr << "File with van der Waals radiuses '" << maximals << "' contains multiple definitions of element '" << line.substr(0, tab) << "'" << std::endl;
            }
          }
        }
        input.close();
      }

      // NOTE: Not save in general approach
      //~RelativeAminoacidFeature() {
      //  if (FEATURE != nullptr) {
      //    delete FEATURE;
      //  }
      //}

      std::string title() override { return "r" + FEATURE->title(); }

      void init(Protein* protein) override {
        FEATURE->init(protein);
      }

      float feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        // NOTE: In this stage, it is expected that each query is on different aminoacid
        // TODO: Consider the first aminoacid
        // TODO: That it is useless to set it in this and inner feature too
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto it = MAXIMAL.find(ITERATOR->residue_name());
        if (it == MAXIMAL.end()) {
          std::cerr << ITERATOR->getProteinName() << ": Residue type '" << ITERATOR->residue_name() << "' is not known" << std::endl;
          return UNDEFINED;
        }
        float ret = FEATURE->feature(model, chain, aminoacid);
        if (ret == FEATURE->UNDEFINED) {
          return UNDEFINED;
        } else {
          return ret/it->second;
        }
      }
    };

    class SingleAminoacidFeature : public Feature<std::string> {
      private:
      std::map<std::string, std::string> TRANSFORMATIONS;

      public:
      SingleAminoacidFeature(ProteinIterator* iterator, std::string table) : Feature(iterator) {
        std::ifstream transformations(table);
        std::string line;
        while (!transformations.eof() && std::getline(transformations, line)) {
          size_t tab = line.find('\t');
          if (tab == line.npos) {
            std::cerr << "'" << line << "' is not a valid transformation line." << std::endl;
          } else {
            auto ins = TRANSFORMATIONS.insert({line.substr(0, tab), line.substr(tab+1)});
            if (!ins.second) {
              std::cerr << "Multiple lines have a same key: '" << ins.first->first << '\t' << ins.first->second << "' vs. '" << line << "'" << std::endl;
            }
          }
        }
      }

      std::string title() override { return "aminoacid"; }

      void init(Protein* protein) override {
        ITERATOR->init(protein);
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        // NOTE: In this stage, it is expected that each query is on different aminoacid
        // TODO: Consider the first aminoacid
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto it = TRANSFORMATIONS.find(ITERATOR->residue_name());
        if (it == TRANSFORMATIONS.end()) {
          std::cerr << ITERATOR->getProteinName() << ": Residue type '" << ITERATOR->residue_name() << "' is not known" << std::endl;
          return UNDEFINED;
        }
        return it->second;
      }
    };

    class BasicAminoacidFeature : public Feature<std::string> {
      public:
      BasicAminoacidFeature(ProteinIterator* iterator) : Feature(iterator) { }

      std::string title() override { return "aminoacid"; }

      void init(Protein* protein) override {
        ITERATOR->init(protein);
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        // NOTE: In this stage, it is expected that each query is on different aminoacid
        // TODO: Consider the first aminoacid
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        return ITERATOR->residue_name();
      }
    };

    class CoordinateFeature : public Feature<std::string> {
      public:
      static Coordinate parse(std::string coordinates) {
        size_t first = coordinates.find(' ');
        if (first == coordinates.npos) {
          throw common::exception::TitledException("Invalid format of coordinate feature: '" + coordinates + "'");
        }
        size_t last = coordinates.rfind(' ');
        if (last == first) {
          throw common::exception::TitledException("Invalid format of coordinate feature, missing third coordinate: '" + coordinates + "'");
        }
        return std::make_tuple(std::stod(coordinates.substr(0, first)), std::stod(coordinates.substr(first+1, last-first-1)), std::stod(coordinates.substr(last+1)));
      }

      CoordinateFeature(ProteinIterator* iterator) : Feature(iterator) { }

      std::string title() override { return "coordinate"; }

      void init(Protein* protein) override {
        ITERATOR->init(protein);
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAtom()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAtom()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        // NOTE: In this stage, it is expected that each query is on different aminoacid
        // TODO: Consider the first aminoacid
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        if (!ITERATOR->resetAtom()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for aminoacid '" << aminoacid << "' chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        do {
          if (ITERATOR->atom_name() == "CA") {
            if (!ITERATOR->computeCharacteristics()) {
              std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid characteristics for Carbon_alpha atom of aminoacid '" << aminoacid << "' chain '" << chain << "' in model '" << model << "'" << std::endl;
              return UNDEFINED;
            }
            Coordinate coordinates = ITERATOR->coordinates();
            // NOTE: Theoretically, uniqueness of Carbon_alpha should be tested. However, in a valid file it will be a vaste of time as C_alpha should be one of the first atoms
            return std::to_string(std::get<0>(coordinates)) + " "+ std::to_string(std::get<1>(coordinates)) + " "+std::to_string(std::get<2>(coordinates));
          }
        } while (ITERATOR->nextAtom());
        // TODO: Consider some estimation of Carbon_alpha position based on other atoms
        std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid Carbon_alpha atom for aminoacid '" << aminoacid << "' chain '" << chain << "' in model '" << model << "'" << std::endl;
        return UNDEFINED;
      }
    };

    class CompositionFeature : public Feature<std::string> {
      public:
      CompositionFeature(ProteinIterator* iterator) : Feature(iterator) { }

      std::string title() override { return "composition"; }

      void init(Protein* protein) override {
        ITERATOR->init(protein);
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        std::map<std::string, size_t> composition;
        if (ITERATOR->resetAtom()) {
          do {
            if (!(ITERATOR->element() == "H" || ITERATOR->element() == "D")) {
              auto it = composition.find(ITERATOR->getAtomName());
              if (it == composition.end()) {
                composition[ITERATOR->getAtomName()] = 1;
              } else {
                std::cerr << "Protein '" << ITERATOR->getProteinName() << "' contains multiple atoms with the same name '" << ITERATOR->getAtomName() << "' in aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
                ++(it->second);
              }
            }
          } while (ITERATOR->nextAtom());
        }
        std::stringstream buffer;
        for (auto composition_it = composition.begin(); composition_it != composition.end(); ++composition_it) {
          for (size_t i = 0; i < composition_it->second; ++i) {
            buffer << composition_it->first << ' ';
          }
        }
        std::string ret = buffer.str();
        if (ret.empty()) {
          return "";
        } else {
          ret.pop_back();
          return ret;
        }
      }
    };


    class InterfaceFeature : public Feature<std::string> {
      protected:
      std::map<std::string, double> ELEMENTS;
      // NOTE: The opposite approach with suming element distances during testing allows default distance for unknown elements,
      // however the current approach should be quicker and in the case of unknown elements, user can rerun extraction with modified radiuses table
      std::map<std::string, std::map<std::string, double>> DISTANCES;
      // Model_name => [Chain_name => [Aminoacid_name => interface]]
      std::map<std::string, std::map<std::string, std::map<std::string, bool>>> INTERFACES;

      public:
      InterfaceFeature(ProteinIterator* iterator, std::string distances, double distance) : Feature(iterator) {
        std::ifstream input(distances);
        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          size_t tab = line.find('\t');
          if (tab == line.npos) {
            std::cerr << "File with van der Waals radiuses '" << distances << "' contains invalid line '" << line << "'" << std::endl;
          } else {
            auto ins = ELEMENTS.insert({line.substr(0, tab), std::stod(line.substr(tab+1))});
            if (!ins.second) {
              std::cerr << "File with van der Waals radiuses '" << distances << "' contains multiple definitions of element '" << line.substr(0, tab) << "'" << std::endl;
            }
          }
        }
        input.close();
        for (auto it_1 = ELEMENTS.begin(); it_1 != ELEMENTS.end(); ++it_1) {
          auto &dist = DISTANCES[it_1->first];
          for (auto it_2 = ELEMENTS.begin(); it_2 != ELEMENTS.end(); ++it_2) {
            dist[it_2->first] = (it_1->second < 0 || it_2->second < 0) ? -1 : std::pow(it_1->second + distance + it_2->second, 2);
          }
        }
      }

      std::string title() override { return "interface"; }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        auto model_it = INTERFACES.find(model);
        if (model_it == INTERFACES.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto chain_it = model_it->second.find(chain);
        if (chain_it == model_it->second.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        auto aminoacid_it = chain_it->second.find(aminoacid);
        if (aminoacid_it == chain_it->second.end()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" + aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        return aminoacid_it->second ? "I" : "N";
      }
    };

    class InterfaceShortFeature : public InterfaceFeature {
      public:
      InterfaceShortFeature(ProteinIterator* iterator, std::string distances, double distance) : InterfaceFeature(iterator, distances, distance) { }

      void init(Protein* protein) override {
        INTERFACES.clear();

        std::set<std::string> warned;
        ITERATOR->init(protein);
        if (ITERATOR->resetModel()) {
          do {
            std::map<std::string, std::map<std::string, std::list<std::pair<std::string, Coordinate>>>> coordinates;
            if (ITERATOR->resetChain()) {
              do {
                auto &chain = coordinates[ITERATOR->getChainName()];
                if (ITERATOR->resetAminoacid()) {
                  do {
                    auto &aminoacid = chain[ITERATOR->getAminoacidName()];
                    if (ITERATOR->resetAtom()) {
                      do {
                        if (DISTANCES.find(ITERATOR->element()) == DISTANCES.end()) {
                          if (warned.insert(ITERATOR->element()).second) {
                            std::cerr << ITERATOR->getProteinName() << ": chemical element '" << ITERATOR->element() << "' does not have specified van der Waals radius" << std::endl;
                          }
                        } else if (ITERATOR->computeCharacteristics()) {
                          aminoacid.push_back({ITERATOR->element(), ITERATOR->coordinates()});
                        } else {
                          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not contains a valid characteristic for an atom '" << ITERATOR->getAtomName() << "' in aminoacid '"
                                    << ITERATOR->getAminoacidName() << "' in chain '" << ITERATOR->getChainName() << "' in model '" << ITERATOR->getModelName() << "'" << std::endl;
                        }
                      } while (ITERATOR->nextAtom());
                    }
                  } while (ITERATOR->nextAminoacid());
                }
              } while (ITERATOR->nextChain());
            }

            // TODO: Consider diagonal calculation, it could save up to 1/2 time, however it requires squared memory
            auto &model = INTERFACES[ITERATOR->getModelName()];
            for (auto chain_1_it = coordinates.begin(); chain_1_it != coordinates.end(); ++chain_1_it) {
              auto &chain = model[chain_1_it->first];
              for (auto aminoacid_1_it = chain_1_it->second.begin(); aminoacid_1_it != chain_1_it->second.end(); ++aminoacid_1_it) {
                auto &aminoacid = chain[aminoacid_1_it->first];
                aminoacid = false;
                // TODO: Test different level of nesting atom_1_it in *_2_it
                for (auto atom_1_it = aminoacid_1_it->second.begin(); atom_1_it != aminoacid_1_it->second.end(); ++atom_1_it) {
                  auto distances = DISTANCES.find(atom_1_it->first);
                  /*// NOTE: This should not happen thanks check during coordinates filling
                  if (distances == DISTANCES.end()) {
                    std::cerr << ITERATOR->getProteinName() << ": chemical element '" << atom_1_it->first << "' does not have specified van der Waals radius" << std::endl;
                    continue;
                  }*/
                  for (auto chain_2_it = coordinates.begin(); chain_2_it != coordinates.end(); ++chain_2_it) {
                    if (chain_1_it->first != chain_2_it->first) {
                      for (auto aminoacid_2_it = chain_2_it->second.begin(); aminoacid_2_it != chain_2_it->second.end(); ++aminoacid_2_it) {
                        // TODO: Test different level of nesting atom_1_it in *_2_it
                        for (auto atom_2_it = aminoacid_2_it->second.begin(); atom_2_it != aminoacid_2_it->second.end(); ++atom_2_it) {
                          auto distance = distances->second.find(atom_2_it->first);
                          /*// NOTE: This should not happen thanks check during coordinates filling if DISTANCES are symmetric
                          if (distance == distances->second.end()) {
                            std::cerr << ITERATOR->getProteinName() << ": chemical element '" << atom_2_it->first << "' does not have specified van der Waals radius" << std::endl;
                            continue;
                          }*/
                          if (std::pow(std::get<0>(atom_1_it->second) - std::get<0>(atom_2_it->second), 2) + std::pow(std::get<1>(atom_1_it->second) - std::get<1>(atom_2_it->second), 2)
                            + std::pow(std::get<2>(atom_1_it->second) - std::get<2>(atom_2_it->second), 2) <= distance->second) {
                            aminoacid = true;
                            goto found;
                          }
                        }
                      }
                    }
                  }
                }
found:;
              }
            }
          } while (ITERATOR->nextModel());
        }
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }
    };

    class InterfaceMediumFeature : public InterfaceFeature {
      private:
      std::map<std::string, double> DISTANCES_MAX;

      public:
      InterfaceMediumFeature(ProteinIterator* iterator, std::string distances, double distance) : InterfaceFeature(iterator, distances, distance) {
        double max = -1;
        for (auto elements_it = ELEMENTS.begin(); elements_it != ELEMENTS.end(); ++elements_it) {
          if (elements_it->second > max) {
            max = elements_it->second;
          }
        }
        for (auto elements_it = ELEMENTS.begin(); elements_it != ELEMENTS.end(); ++elements_it) {
          // NOTE: Small value added because of rounding as this is for a heuristic only and -1 vs -0.9999 does not matter in the case of test dm<0
          DISTANCES_MAX[elements_it->first] = max + distance + elements_it->second + 0.0001;
        }
      }

      void init(Protein* protein) override {
        INTERFACES.clear();

        std::set<std::string> warned;
        ITERATOR->init(protein);
        if (ITERATOR->resetModel()) {
          do {
            std::map<std::string, std::map<std::string, std::list<std::pair<std::string, Coordinate>>>> coordinates;
            // TODO: Compare effectivity of multimap with effectivity of map of vectors/lists/sets... and with effectivity of two-dimensional index (x;y)
            std::map<std::string, std::multimap<double, std::pair<std::string, Coordinate>>> counterparts;
            if (ITERATOR->resetChain()) {
              do {
                std::string chain_name = ITERATOR->getChainName();
                auto &chain = coordinates[chain_name];
                auto &counterpart = counterparts[chain_name];
                if (ITERATOR->resetAminoacid()) {
                  do {
                    auto &aminoacid = chain[ITERATOR->getAminoacidName()];
                    if (ITERATOR->resetAtom()) {
                      do {
                        if (DISTANCES.find(ITERATOR->element()) == DISTANCES.end()) {
                          if (warned.insert(ITERATOR->element()).second) {
                            std::cerr << ITERATOR->getProteinName() << ": chemical element '" << ITERATOR->element() << "' does not have specified van der Waals radius" << std::endl;
                          }
                        } else if (ITERATOR->computeCharacteristics()) {
                          auto position = ITERATOR->coordinates();
                          aminoacid.push_back({ITERATOR->element(), position});
                          counterpart.insert({std::get<0>(position),{ITERATOR->element(), position}});
                        } else {
                          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not contains a valid characteristic for an atom '" << ITERATOR->getAtomName() << "' in aminoacid '"
                            << ITERATOR->getAminoacidName() << "' in chain '" << ITERATOR->getChainName() << "' in model '" << ITERATOR->getModelName() << "'" << std::endl;
                        }
                      } while (ITERATOR->nextAtom());
                    }
                  } while (ITERATOR->nextAminoacid());
                }
              } while (ITERATOR->nextChain());
            }

            // TODO: Consider diagonal calculation, it could save up to 1/2 time, however it requires squared memory
            auto &model = INTERFACES[ITERATOR->getModelName()];
            for (auto chain_1_it = coordinates.begin(); chain_1_it != coordinates.end(); ++chain_1_it) {
              auto &chain = model[chain_1_it->first];
              for (auto aminoacid_1_it = chain_1_it->second.begin(); aminoacid_1_it != chain_1_it->second.end(); ++aminoacid_1_it) {
                auto &aminoacid = chain[aminoacid_1_it->first];
                aminoacid = false;
                // TODO: Test different level of nesting atom_1_it in *_2_it
                for (auto atom_1_it = aminoacid_1_it->second.begin(); atom_1_it != aminoacid_1_it->second.end(); ++atom_1_it) {
                  // NOTE: This is possible as unknown elements are skipped during filling of map
                  auto distances = DISTANCES[atom_1_it->first];
                  double max = DISTANCES_MAX[atom_1_it->first];
                  if (max < 0) {
                    continue;
                  }
                  for (auto chain_2_it = counterparts.begin(); chain_2_it != counterparts.end(); ++chain_2_it) {
                    if (chain_1_it->first != chain_2_it->first) {
                      for (auto atom_2_it = chain_2_it->second.lower_bound(std::get<0>(atom_1_it->second)-max);
                           atom_2_it != chain_2_it->second.upper_bound(std::get<0>(atom_1_it->second)+max); ++atom_2_it) {
                        if (std::pow(std::get<0>(atom_1_it->second) - std::get<0>(atom_2_it->second.second), 2) + std::pow(std::get<1>(atom_1_it->second) - std::get<1>(atom_2_it->second.second), 2)
                          + std::pow(std::get<2>(atom_1_it->second) - std::get<2>(atom_2_it->second.second), 2) <= distances[atom_2_it->second.first]) {
                          aminoacid = true;
                          goto found;
                        }
                      }
                    }
                  }
                }
found:;
              }
            }
          } while (ITERATOR->nextModel());
        }
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }
    };

    class InterfaceLargeFeature : public InterfaceFeature {
      private:
      double DISTANCE;

      public:
      InterfaceLargeFeature(ProteinIterator* iterator, std::string distances, double distance) : DISTANCE(distance), InterfaceFeature(iterator, distances, distance) { }

      void init(Protein* protein) override {
        INTERFACES.clear();

        std::set<std::string> warned;
        ITERATOR->init(protein);
        if (ITERATOR->resetModel()) {
          do {
            std::map<std::string, std::map<std::string, std::list<std::pair<std::string, Coordinate>>>> coordinates;
            if (ITERATOR->resetChain()) {
              do {
                auto &chain = coordinates[ITERATOR->getChainName()];
                if (ITERATOR->resetAminoacid()) {
                  do {
                    auto &aminoacid = chain[ITERATOR->getAminoacidName()];
                    if (ITERATOR->resetAtom()) {
                      do {
                        if (ELEMENTS.find(ITERATOR->element()) == ELEMENTS.end()) {
                          if (warned.insert(ITERATOR->element()).second) {
                            std::cerr << ITERATOR->getProteinName() << ": chemical element '" << ITERATOR->element() << "' does not have specified van der Waals radius" << std::endl;
                          }
                        } else if (ITERATOR->computeCharacteristics()) {
                          aminoacid.push_back({ITERATOR->element(), ITERATOR->coordinates()});
                        } else {
                          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not contains a valid characteristic for an atom '" << ITERATOR->getAtomName() << "' in aminoacid '"
                            << ITERATOR->getAminoacidName() << "' in chain '" << ITERATOR->getChainName() << "' in model '" << ITERATOR->getModelName() << "'" << std::endl;
                        }
                      } while (ITERATOR->nextAtom());
                    }
                  } while (ITERATOR->nextAminoacid());
                }
              } while (ITERATOR->nextChain());
            }

            std::map<std::string, Octree*> octrees;
            for (auto chain_it = coordinates.begin(); chain_it != coordinates.end(); ++chain_it) {
              std::vector<std::pair<Coordinate, double>> chain;
              for (auto residue_it = chain_it->second.begin(); residue_it != chain_it->second.end(); ++residue_it) {
                for (auto atom_it = residue_it->second.begin(); atom_it != residue_it->second.end(); ++atom_it) {
                  double distance = ELEMENTS[atom_it->first];
                  if (distance >= 0) {
                    chain.push_back(std::make_pair(atom_it->second, distance + DISTANCE));
                  }
                }
              }
              octrees.insert({{chain_it->first, new Octree(chain)}});
            }

            // TODO: Consider diagonal calculation, it could save up to 1/2 time, however it requires squared memory
            auto &model = INTERFACES[ITERATOR->getModelName()];
            for (auto chain_1_it = coordinates.begin(); chain_1_it != coordinates.end(); ++chain_1_it) {
              auto &chain = model[chain_1_it->first];
              for (auto aminoacid_it = chain_1_it->second.begin(); aminoacid_it != chain_1_it->second.end(); ++aminoacid_it) {
                auto &aminoacid = chain[aminoacid_it->first];
                aminoacid = false;
                // TODO: Test different level of nesting atom_1_it in *_2_it
                for (auto atom_it = aminoacid_it->second.begin(); atom_it != aminoacid_it->second.end(); ++atom_it) {
                  // NOTE: This should exist happen thanks check during coordinates filling
                  auto distance = ELEMENTS[atom_it->first];
                  if (distance >= 0) {
                    for (auto chain_2_it = octrees.begin(); chain_2_it != octrees.end(); ++chain_2_it) {
                      if (chain_1_it->first != chain_2_it->first && chain_2_it->second->contact(atom_it->second, distance)) {
                        aminoacid = true;
                        goto found;
                      }
                    }
                  }
                }
found:;
              }
            }
            for (auto octrees_it = octrees.begin(); octrees_it != octrees.end(); ++octrees_it) {
              delete octrees_it->second;
            }
          } while (ITERATOR->nextModel());
        }
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }
    };


    class TemperatureFeature : public Feature<std::string> {
      public:
      TemperatureFeature(ProteinIterator* iterator) : Feature(iterator) { }

      std::string title() override { return "temperature"; }

      void init(Protein* protein) override {
        ITERATOR->init(protein);
        ITERATOR->resetModel();
        ITERATOR->resetChain();
        ITERATOR->resetAminoacid();
      }

      std::string feature(const std::string &model, const std::string &chain, const std::string &aminoacid) override {
        // TODO: Compare effectivity of this heuristics; also consider a cacheing of names
        if (ITERATOR->getModelName() != model) {
          if (!((ITERATOR->nextModel() && ITERATOR->getModelName() == model) || ITERATOR->setModel(model))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetChain()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid chain for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAtom()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        if (ITERATOR->getChainName() != chain) {
          if (!((ITERATOR->nextChain() && ITERATOR->getChainName() == chain) || ITERATOR->setChain(chain))) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAminoacid()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid aminoacid for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
          if (!ITERATOR->resetAtom()) {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for chain '" << chain << "' in model '" << model << "'" << std::endl;
            return UNDEFINED;
          }
        }
        // NOTE: In this stage, it is expected that each query is on different aminoacid
        // TODO: Consider the first aminoacid
        if (!((ITERATOR->nextAminoacid() && ITERATOR->getAminoacidName() == aminoacid) || ITERATOR->setAminoacid(aminoacid))) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a aminoacid '" << aminoacid << "' in chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        if (!ITERATOR->resetAtom()) {
          std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not have a valid atom for aminoacid '" << aminoacid << "' chain '" << chain << "' in model '" << model << "'" << std::endl;
          return UNDEFINED;
        }
        // TODO: Well, there are several possibilities, how to determine temperature - C_alpha; average; average side chain; ignore hydrogens if available etc.
        float sum = 0;
        size_t count = 0;
        do {
          if (ITERATOR->computeCharacteristics()) {
            ++count;
            sum += ITERATOR->temperature();
          } else {
            std::cerr << "Protein '" << ITERATOR->getProteinName() << "' does not contains a valid characteristic for an atom '" << ITERATOR->getAtomName() << "' in aminoacid '"
                      << ITERATOR->getAminoacidName() << "' in chain '" << ITERATOR->getChainName() << "' in model '" << ITERATOR->getModelName() << "'" << std::endl;
          }
        } while (ITERATOR->nextAtom());
        return std::to_string(sum/count);
      }
    };

    class Features {
      private:
      BasicFilter* FILTER;
      Index INDEX;
      ProteinIterator* ITERATOR;
      std::map<std::string, std::string> FILES;

      // Add the given path together with its protein identifier in dictionary
      void index_pdb_file(const std::string &file) {
        try {
          auto ins = FILES.insert({ProteinParser::parse_id(file), file});
        } catch (...) { }
      }

      public:
      // index: path to a index file
      // iterator: how to interpret the index file
      Features(std::string index, ProteinIterator* iterator, BasicFilter* filter) : INDEX(index), ITERATOR(iterator), FILTER(filter) { }

      // Add all pdb files the given path together with their protein identifiers in dictionary
      void index_pdb(const std::string &item) {
        if (common::filesystem::is_directory(item)) {
          common::filesystem::RecursiveDirectoryFileIterator file_iterator(item);
          if (file_iterator.has_file()) {
            do {
              index_pdb_file(file_iterator.filename());
            } while (file_iterator.has_next());
          } else {
            std::cerr << "There is no file in the given directory." << std::endl;
          }
        } else {
          index_pdb_file(item);
        }
      }

      // Extract all features and store them in the file
      // file: path of the output file
      // features: feature extractors
      void extract_features(std::string &file, std::vector<Feature<std::string>*> &features) {
        if (!INDEX.reset()) {
          throw common::exception::TitledException("The index file is empty or it is not possible to read it.");
        }
        
        // Header line
        std::ofstream output(file);
        for (size_t i = 0; i < features.size(); i++) {
          if (i != 0) {
            output << '\t';
          }
          output << features[i]->title();
        }
        output << '\n';

        // Data
        size_t last = 0;
        // TODO: Is this necessary to prevent destruction until it is not problem?
        // NOTE: Expects that no protein has empty id_code and that the index file is valid (id_code in the first line)
        Protein protein;
        do {
          if (protein.ID_CODE != INDEX.protein()) {
            protein = ProteinParser::parse_protein(FILES[INDEX.protein()], FILTER);
            for (size_t i = 0; i < features.size(); i++) {
              features[i]->init(&protein);
            }
          }
          std::stringstream buffer;
          bool write = false;
          for (size_t i = 0; i < features.size(); i++) {
            if (i == 0) {
            } else {
              buffer << '\t';
            }
            std::string feature = features[i]->feature(INDEX.model(), INDEX.chain(), INDEX.aminoacid());
            if (feature != features[i]->UNDEFINED) {
              write = true;
              buffer << feature;
            }
          }
          if (write) {
            output << buffer.str();
            if (++last != INDEX.index()) {
              last = INDEX.index();
              output << '\t' << INDEX.index();
            }
            output << '\n';
          }
        } while (INDEX.next());

        output.flush();
        output.close();
      }
    };

    class FeatureReader {
      private:
      std::string NAME;
      std::ifstream INPUT;

      int COLUMN;
      int COLUMNS;

      size_t INDEX;
      std::string VALUE;

      public:
      FeatureReader(std::string file, std::string column) : NAME(file), INPUT(file), COLUMN(0), COLUMNS(0), INDEX(0), VALUE("") {
        std::string line;
        if (!std::getline(INPUT, line)) {
          throw common::exception::TitledException("The feature file '" + file + "' is empty");
        }
        std::stringstream parts(line);
        std::string part;
        while (std::getline(parts, part, '\t')) {
          if (part == column) {
            if (COLUMN < COLUMNS) {
              throw common::exception::TitledException("Features file '" + file + "' contains multiple columns with header '" + column + "'");
            }
          } else if (COLUMN == COLUMNS) {
            ++COLUMN;
          }
          ++COLUMNS;
        }
        if (COLUMN == COLUMNS) {
          throw common::exception::TitledException("Features file '" + file + "' contains multiple columns with header '" + column + "'");
        }
      }

      bool next(size_t index) {
        if (index < INDEX) {
          return false;
        }
        while (index > INDEX) {
          std::string line;
          if (INPUT.eof() || !std::getline(INPUT, line)) {
            INDEX = std::numeric_limits<int>::max();
            return false;
          }
          std::stringstream parts(line);
          std::string part;
          for (size_t i = 0; i <= COLUMN; i++) {
            if (!std::getline(parts, part, '\t')) {
              throw common::exception::TitledException("Unexpected format of features file '" + NAME + "': not enough elements in line '" + line + "'");
            }
          }
          VALUE = part;
          for (size_t i = COLUMN+1; i < COLUMNS; i++) {
            if (!std::getline(parts, part, '\t')) {
              throw common::exception::TitledException("Unexpected format of features file '" + NAME + "': not enough elements in line '" + line + "'");
            }
          }
          if (std::getline(parts, part, '\t')) {
            INDEX = std::stoi(part);
          } else {
            ++INDEX;
          }
          if (std::getline(parts, part, '\t')) {
            throw common::exception::TitledException("Unexpected format of features file '" + NAME + "': too many elements in line '" + line + "'");
          }
        }
        return index == INDEX && VALUE.size() > 0;
      }

      std::string value() {
        return VALUE;
      }
    };
    class FeaturesReader {
      private:
      std::string NAME;
      std::ifstream* INPUT;

      // NOTE: This is not optimal, however, header should occupy only a negligible piece of memory and thus this approach is prefered because of clarity
      std::vector<std::string> HEADERS;

      size_t INDEX;
      std::vector<std::string> LINE;

      public:
      FeaturesReader(std::string file) : NAME(file), INPUT(new std::ifstream(file)), INDEX(0) {
        std::string line;
        if (!std::getline(*INPUT, line)) {
          throw common::exception::TitledException("The feature file '" + file + "' is empty");
        }
        std::stringstream parts(line);
        std::string part;
        while (std::getline(parts, part, '\t')) {
          HEADERS.push_back(part);
        }
      }
			~FeaturesReader() {
				delete INPUT;
			}

      bool next_line(size_t index) {
        if (index < INDEX) {
          return false;
        }
        while (index > INDEX) {
          LINE.clear();
          std::string line;
          if (INPUT->eof() || !std::getline(*INPUT, line)) {
            INDEX = std::numeric_limits<int>::max();
            return false;
          }
          std::stringstream parts(line);
          std::string part;
          for (size_t i = 0; i < HEADERS.size(); i++) {
            if (!std::getline(parts, part, '\t')) {
              throw common::exception::TitledException("Unexpected format of features file '" + NAME + "': not enough elements in line '" + line + "'");
            }
            LINE.push_back(part);
          }
          if (std::getline(parts, part, '\t')) {
            INDEX = std::stoi(part);
          } else {
            ++INDEX;
          }
          if (std::getline(parts, part, '\t')) {
            throw common::exception::TitledException("Unexpected format of features file '" + NAME + "': too many elements in line '" + line + "'");
          }
        }
        return index == INDEX;
      }

      bool next_line() {
        while (INDEX != std::numeric_limits<int>::max()) {
          if (next_line(INDEX+1)) {
            return true;
          }
        }
        return false;
      }
      size_t index() {
        return INDEX;
      }

      size_t size() {
        return HEADERS.size();
      }

      std::string header(size_t i) {
        return HEADERS[i];
      }

      std::string value(size_t i) {
        return LINE[i];
      }
    };
  }
}