#pragma once

#include "index.h"
#include "features.h"
#include "subgraphs.h"
#include "../common/exception.h"
#include "../common/filesystem.h"
#include "../common/graph.h"
#include "../common/string.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <limits.h>

namespace inspire {
  namespace backend {
    class AbstractStream {
      public:
      AbstractStream() { }
      virtual ~AbstractStream() { }
      // NOTE: For performancy reasons, names of features are send accepted during initialization only and now only values are in <residue>.
      //       Thus, it is expected, feature values are in the same order as headers during initialization
      virtual void write(int id, std::vector<std::string> &residue, std::vector<bool> &fingerprint) = 0;
      virtual void finalize() = 0;
    };

    class QueryStream : public AbstractStream {
      // A file to store fingerprints.
      std::ofstream STREAM;
      // Headers of features to label feature values in the file.
      std::vector<std::string> HEADERS;
      public:
      QueryStream(std::string path, std::vector<std::string> &headers) : HEADERS(headers) {
        if (path.size() == 0 || path.back() == common::filesystem::directory_separator) {
          path += "query.fit";
        } else if (common::filesystem::is_directory(path)) {
          path += "/query.fit";
        } else {
          if (!common::string::ends_with(path, ".fit")) {
            path += ".fit";
          }
        }
        STREAM.open(path);
      }
      ~QueryStream() { }
      void write(int id, std::vector<std::string> &residue, std::vector<bool> &fingerprint) {
        STREAM << id << '\t';
        for (size_t i = 0; i < HEADERS.size(); i++) {
          STREAM << HEADERS[i] << ':' << residue[i] << '\t';
        }
        for (bool position : fingerprint) {
          STREAM << (position ? 1 : 0);
        }
        STREAM << '\n';
      }
      void finalize() {
        STREAM.flush();
        STREAM.close();
      }
    };

    // How is length of fingerprints encoded in file header
    // NOTE: The bigger the type is, the longer fingerprints can be used; while the smaller the type is, the less memory is used
    // NOTE: Actually, uint16_t should be enough
    typedef uint32_t size_f;

    class KnowledgebaseStream : public AbstractStream {
      private:
      std::string PATH;
      std::vector<std::string> HEADERS;
      /* List of files to store fingerprints */
      std::map<std::vector<std::string>, std::ofstream> STREAMS;
      /* Size of fingerprints in bytes */
      size_f LENGTH;

      std::map<std::vector<std::string>, std::ofstream>::iterator init(std::vector<std::string> &residues) {
        std::stringstream buffer;
        buffer << PATH;
        for (size_t i = 0; i < HEADERS.size() - 1; ++i) {
          std::string dir = HEADERS[i] + '_' + residues[i];
          if (!common::filesystem::is_portable_directory(dir)) {
            // NOTE: Previously created directories and/or files are not deleted.
            throw common::exception::TitledException("'" + dir + "' is not a valid portable directory name.");
          }
          buffer << dir << common::filesystem::directory_separator;
        }
        std::string filepath = buffer.str();
        if (!common::filesystem::exists(filepath) && !common::filesystem::create_directory_recursive(filepath)) {
          // NOTE: Previously created directories and/or files are not deleted.
          throw common::exception::TitledException("It is not possible to create a directory structure '" + buffer.str() + "'");
        }
        std::string file = HEADERS.back() + '_' + residues.back() + ".fin";
        if (!common::filesystem::is_portable_file(file)) {
          // NOTE: Previously created directories and/or files are not deleted.
          throw common::exception::TitledException("'" + file + "' is not a valid portable file name");
        }
        buffer << file;
        // CHECK: Does not this create a new vector?
        // CHECK: Does not it destroy stream on the end of this method?
        auto stream = STREAMS.insert({residues, std::ofstream(buffer.str(), std::ios::binary)}).first;
        stream->second.write(reinterpret_cast<const char*>(&LENGTH), sizeof(LENGTH));
        return stream;
      }

      public:
      KnowledgebaseStream(std::string path, std::vector<std::string> &headers, size_f length) : PATH(path), HEADERS(headers), LENGTH(length == 0 ? 0 : ((length-1)/CHAR_BIT+1)) {
        if (PATH.size() > 0) {
          if (common::filesystem::exists(PATH)) {
            if (!common::filesystem::is_directory(PATH)) {
              throw common::exception::TitledException("'" + PATH + "' exists but is not a directory.");
            }
          } else if (!common::filesystem::create_directory_recursive(PATH)) {
            throw common::exception::TitledException("It is not possible to create '" + PATH + "'");
          }
          if (PATH.back() != common::filesystem::directory_separator) {
            PATH.push_back(common::filesystem::directory_separator);
          }
        }
      }
      ~KnowledgebaseStream() { }
      void write(int id, std::vector<std::string> &residue, std::vector<bool> &fingerprint) {
        auto stream = STREAMS.find(residue);
        if (stream == STREAMS.end()) {
          stream = init(residue);
        }
        //std::ofstream &stream = STREAMS[residue];
        stream->second.write(reinterpret_cast<const char*>(&id), sizeof(id));
        // TODO: Move into special file to use in prediction for query fingerprints too.
        size_t bytes = 0;
        char bites = 0;
        char byte = 0;
        for (auto fingerprint_it = fingerprint.begin(); fingerprint_it != fingerprint.end(); ++fingerprint_it) {
          if (*fingerprint_it) {
            byte += (1 << bites);
          }
          if (++bites == CHAR_BIT) {
            stream->second.put(byte);
            bites = 0;
            if (++bytes == LENGTH) {
              break;
            }
            byte = 0;
          }
        }
        if (bites > 0) {
          stream->second.write(&byte, sizeof(byte));
          ++bytes;
        }
        if (bytes < LENGTH) {
          throw common::exception::TitledException("The fingerprint no. " + std::to_string(id) + " has an invalid length");
        }
        /*for (; bytes < LENGTH; ++bytes) {
          stream.put(0);
        }*/
      }
      void finalize() {
        for (auto streams_it = STREAMS.begin(); streams_it != STREAMS.end(); ++streams_it) {
          streams_it->second.flush();
          streams_it->second.close();
        }
        //? Is not necessary to destruct streams?
        STREAMS.clear();
      }
    };

    namespace fingerprint {

      // Print debug information about used indexes in fingerprint.
      //#define DEBUG_STDOUT

      namespace graph {

        class graph {
          public:
          typedef std::pair<int, int> edge;
          public:
          void add_undirected_edge(int left, int right) {
            edges.push_back(std::make_pair(left, right));
            this->vertices.insert(left);
            this->vertices.insert(right);
          }
          void add_vertex(int index) {
            this->vertices.insert(index);
          }
          public:
          std::pair<std::set<int>::const_iterator,
            std::set<int>::const_iterator> get_vertices_iterators() const {
            return std::make_pair(vertices.begin(), vertices.end());
          }
          std::pair<std::vector<edge>::const_iterator,
            std::vector<edge>::const_iterator> get_edges_iterators() const {
            return std::make_pair(edges.begin(), edges.end());
          }
          public:
          size_t get_vertices_count() {
            return vertices.size();
          }
          private:
          std::vector<edge> edges;
          std::set<int> vertices;
        };

        class graph_shortest_paths {
          public:
          static const int UNREACHABLE = -1;
          public:
          graph_shortest_paths() : vertices_distances(0) {

          }
          void compute(const graph& graph) {
            clear_data();
            floyd_warshall(graph);
          }
          int get_distance(int left, int right) const {
            common::graph::vertex left_vertex = name_to_vertex.at(left);
            common::graph::vertex right_vertex = name_to_vertex.at(right);
            int distance = vertices_distances[left_vertex][right_vertex];
            // Check for unreachable.
            if (distance == std::numeric_limits<int>::max()) {
              return UNREACHABLE;
            }
            return distance;
          }
          private: // boost::floyd_warshall_all_pairs_shortest_paths
          void floyd_warshall(const graph& graph) {
            // Create names for vertices.
            common::graph::undirected_graph g;
            auto vertices = graph.get_vertices_iterators();
            for (; vertices.first != vertices.second; ++vertices.first) {
              add_to_index(*vertices.first, common::graph::add_vertex(g));
            }
            // Add edges to the graph.				
            auto edges = graph.get_edges_iterators();
            for (; edges.first != edges.second; ++edges.first) {
              const auto& edge = *edges.first;
              auto first = get_index(edge.first);
              auto second = get_index(edge.second);
              common::graph::add_edge(first, second, g);
            }
            vertices_distances = common::graph::distance_matrix(name_to_vertex.size());
            common::graph::distance_matrix_map distances_map(vertices_distances, g);
            // Every edge has weight 1.
            common::graph::const_property_map const_one_property_map(1);
            common::graph::shortest_paths(g, distances_map, const_one_property_map);
          }
          private:
          void clear_data() {
            name_to_vertex.clear();
          }
          void add_to_index(int vertex_name, common::graph::vertex new_vertex) {
            //int new_index = name_to_vertex.size();
            name_to_vertex.insert(std::make_pair(vertex_name, new_vertex));
          }
          common::graph::vertex get_index(int vertex_name) {
            auto iter = name_to_vertex.find(vertex_name);
            if (iter == name_to_vertex.end()) {
              throw std::runtime_error("Missing vertex.");
            } else {
              return iter->second;
            }
          }
          private:
          std::map<int, common::graph::vertex> name_to_vertex;
          common::graph::distance_matrix vertices_distances;
        };

      }

      namespace properties {

        class verticies_pair_property {
          public:
          verticies_pair_property(size_t size) : result_size(size) {

          }
          size_t size() {
            return result_size;
          }
          virtual void initialize(const graph::graph& graph) = 0;
          virtual int64_t compute_descriptor(int left, int right) const = 0;
          protected:
          size_t result_size;
        };

        class vertex_property {
          public:
          vertex_property(size_t size, const std::string& property_name, std::map<std::string, int> &headers)
            : result_size(size) {
            initialize_property_index(property_name, headers);
          }
          size_t size() {
            return result_size;
          }
          virtual int64_t compute_descriptor(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) const = 0;
          private:
          void initialize_property_index(const std::string& property_name, std::map<std::string, int> &headers) {
            auto iter = headers.find(property_name);
            if (iter == headers.end()) {
              std::string message = "Missing HEADER: '" + property_name + "'";
              throw std::runtime_error(message.c_str());
            } else {
              property_index = headers[property_name];
            }
          }
          protected:
          std::string get_value(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) const {
            return residues[vertex].second[property_index];
          }
          protected:
          size_t result_size;
          int property_index;
        };

        class verticies_distance : public verticies_pair_property {
          public:
          using verticies_pair_property::verticies_pair_property;
          public:
          verticies_distance(size_t size) : verticies_pair_property(size) {
            unreachable_value =
#ifdef MAX
              pow(2, size) - 1
#else
              0
#endif
              ;
          }
          void initialize(const graph::graph& graph) {
            shortest_paths.compute(graph);
          }
          int64_t compute_descriptor(int left, int right) const {
            auto distance = shortest_paths.get_distance(left, right);
            if (distance == graph::graph_shortest_paths::UNREACHABLE) {
              return unreachable_value;
            } else {
              return distance;
            }
          }
          private:
          // Value used for disconnected vertices.
          int64_t unreachable_value;
          graph::graph_shortest_paths shortest_paths;
        };

        class value_property : public vertex_property {
          public:
          using vertex_property::vertex_property;
          public:
          int64_t compute_descriptor(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) const {
            return std::stoi(get_value(vertex, residues));
          }
        };

        class mapping_property : public vertex_property {
          public:
          using vertex_property::vertex_property;
          public:
          void add_mapping(std::string source_value, int64_t target_value) {
            value_map.insert(std::make_pair(source_value, target_value));
          }
          public:
          int64_t compute_descriptor(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) const {
            auto index = get_value(vertex, residues);
            return value_map.at(index);
          }
          private:
          std::map<std::string, int64_t> value_map;
        };

        class binning_property : public vertex_property {
          public:
          using vertex_property::vertex_property;
          public:
          void add_bin(float value) {
            bins.push_back(value);
          }
          public:
          int64_t compute_descriptor(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) const {
            float value = std::stof(get_value(vertex, residues));
            for (int index = 0; index < bins.size(); ++index) {
              if (value < bins[index]) {
                return index;
              }
            }
            return bins.size();
          }
          private:
          std::vector<float> bins;
        };

        verticies_pair_property* create_distance_property(common::graph::property_tree& settings) {
          size_t size = settings.get_child("size").get_value<size_t>();
          return new verticies_distance(size);
        }

        vertex_property* create_value_property(common::graph::property_tree& settings, std::map<std::string, int> &headers) {
          size_t size = settings.get_child("size").get_value<size_t>();
          std::string name = settings.get_child("property").get_value<std::string>();
          value_property* calculator = new value_property(size, name, headers);
          return calculator;
        }

        mapping_property* create_mapping_property(common::graph::property_tree& settings, std::map<std::string, int> &headers) {
          size_t size = settings.get_child("size").get_value<size_t>();
          std::string name = settings.get_child("property").get_value<std::string>();
          mapping_property* calculator = new mapping_property(size, name, headers);
          for (const auto& mapping : settings.get_child("map")) {
            calculator->add_mapping(mapping.first, mapping.second.get_value<int>());
          }
          return calculator;
        }

        vertex_property* create_binning_property(boost::property_tree::ptree& settings, std::map<std::string, int> &headers) {
          size_t size = settings.get_child("size").get_value<size_t>();
          std::string name = settings.get_child("property").get_value<std::string>();
          binning_property* calculator = new binning_property(size, name, headers);
          for (const auto& binning : settings.get_child("bins")) {
            calculator->add_bin(binning.second.get_value<float>());
          }
          return calculator;
        }

      }

      typedef std::vector<bool> fingerprint;
      typedef uint32_t size_f;

      struct configuration {
        public:
        size_f size;
        std::vector<std::shared_ptr<properties::verticies_pair_property> > verticies_pair_properties;
        std::vector<std::shared_ptr<properties::vertex_property> > vertex_properties;
        public:
        void load(const std::string& settings_path, std::map<std::string, int>& headers) {
          common::graph::property_tree root;
          common::graph::parse_json(settings_path, root);
          this->size = root.get_child("fingerprint.size").get_value<size_f>();
          for (auto& edge_settings : root.get_child("fingerprint.edge")) {
            add_verticies_pair_property(edge_settings.second);
          }
          for (auto& vertex_settings : root.get_child("fingerprint.vertex")) {
            add_vertex_property(vertex_settings.second, headers);
          }
          return;
        }
        size_t get_edge_size() const {
          size_t edge_size = 0;
          for (auto iter : verticies_pair_properties) {
            edge_size += iter->size();
          }
          return edge_size;
        }
        size_t get_vertex_size() const {
          size_t vertex_size = 0;
          for (auto iter : vertex_properties) {
            vertex_size += iter->size();
          }
          return vertex_size;
        }
        private:
        void add_verticies_pair_property(common::graph::property_tree& settings) {
          std::string type = settings.get_child("type").get_value<std::string>();
          if (type == "distance") {
            this->verticies_pair_properties.push_back(std::shared_ptr<properties::verticies_pair_property>(
              properties::create_distance_property(settings)));
          } else {
            std::string error = "Invalid calculator type: " + type;
            throw std::runtime_error(error.c_str());
          }
        }
        void add_vertex_property(common::graph::property_tree& settings, std::map<std::string, int> &headers) {
          std::string type = settings.get_child("type").get_value<std::string>();
          if (type == "property") {
            this->vertex_properties.push_back(std::shared_ptr<properties::vertex_property>(
              properties::create_value_property(settings, headers)));
          } else if (type == "mapping") {
            this->vertex_properties.push_back(std::shared_ptr<properties::vertex_property>(
              properties::create_mapping_property(settings, headers)));
          } else if (type == "binning") {
            this->vertex_properties.push_back(std::shared_ptr<properties::vertex_property>(
              properties::create_binning_property(settings, headers)));
          } else {
            std::string error = "Invalid calculator type: " + type;
            throw std::runtime_error(error.c_str());
          }
        }
      };

      class fingerprint_calculator {
        public:
        fingerprint_calculator(const configuration& calculator_configuration) : config(calculator_configuration) { }
        fingerprint compute(graph::graph& graph_object, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) {
          initialize_properties(graph_object);
#ifdef DEBUG_STDOUT
          std::cout << "fingeprint_calculator.compute" << std::endl;
#endif
          fingerprint compute_fingerprint;
          compute_fingerprint.resize(config.size, false);
          auto left_iterator = graph_object.get_vertices_iterators();
          for (; left_iterator.first != left_iterator.second; ++left_iterator.first) {
            auto right_iterator = graph_object.get_vertices_iterators();
            for (; right_iterator.first != right_iterator.second; ++right_iterator.first) {
              int left = *left_iterator.first;
              int right = *right_iterator.first;
#ifndef DIAGONAL
              if (left == right) {
                continue;
              }
#endif
              int64_t index = compute_index(graph_object, left, right, residues);
              compute_fingerprint[index] = true;
            }
          }
          return compute_fingerprint;
        }
        private:
        void initialize_properties(graph::graph& graph_object) {
          for (auto property : config.verticies_pair_properties) {
            property->initialize(graph_object);
          }
        }
        int64_t compute_index(graph::graph& graph_object, int left, int right, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) {
          int64_t left_code = compute_vertex_code(left, residues);
          int64_t edge_code = compute_vertices_pair_code(left, right);
          int64_t right_code = compute_vertex_code(right, residues);
          int64_t index = compose_index(left_code, edge_code, right_code);
#ifdef DEBUG_STDOUT
          std::cout << "\t" << left << " " << right << " " <<
            index << " : " << left_code << " " << edge_code << " " << right_code <<
            std::endl;
#endif
          index = index % config.size;
          return index;
        }
        int64_t compose_index(int64_t left_code, int64_t edge_code, int64_t right_code) {
          return
#ifdef ENN
          (edge_code << (2 * config.get_vertex_size())) +
            (left_code << config.get_vertex_size()) +
            right_code
#elif defined(NNE)
            (left_code << (config.get_edge_size() + config.get_vertex_size())) +
            (right_code << config.get_edge_size()) +
            edge_code
#else // NEN
            (left_code << (config.get_edge_size() + config.get_vertex_size())) +
            (edge_code << config.get_vertex_size()) +
            right_code
#endif	
            ;
        }
        int64_t compute_vertex_code(int vertex, std::map<int, std::pair<std::vector<int>, std::vector<std::string>>> &residues) {
          int64_t code = 0;
          size_t size = 0;
          for (auto property : config.vertex_properties) {
            code += property->compute_descriptor(vertex, residues) << size;
            size += property->size();
          }
          return code;
        }
        int64_t compute_vertices_pair_code(int left, int right) {
          int64_t code = 0;
          size_t size = 0;
          for (auto property : config.verticies_pair_properties) {
            code += property->compute_descriptor(left, right) << size;
            size += property->size();
          }
          return code;
        }
        private:
        const configuration& config;
      };

    }

    enum class FingerprintFormat { Binary, Text };

    class FingerprintWriter {
      private:
      Index INDEX;
      SubgraphReader SUBGRAPHS;
      AbstractStream* OUTPUT;
      // This is a definition of column names for the vector of features of the variable RESIDUES.
      // This splitted representation is 3x space saving and 15% quicker than single map<int, map<string, string> >.
      // Actually there is no space difference between map and unordered_map, but map is 10% quicker.
      // CHECK: Consider boost::bimap
      std::vector<std::string> HEADERS;
      std::map<std::string, int> HEADERS_MAP;

      std::vector<FeaturesReader> FEATURES;

      void process(std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > &model, fingerprint::fingerprint_calculator &calculator) {
        for (auto model_it = model.begin(); model_it != model.end(); ++model_it) {
          while (SUBGRAPHS.next(model_it->first)) {
            // Read verticies and build graph object.
            fingerprint::graph::graph graph_instance;
            const auto &subgraphs = SUBGRAPHS.value();
            for (auto subgraph_it = subgraphs.begin(); subgraph_it != subgraphs.end(); ++subgraph_it) {
              graph_instance.add_vertex(*subgraph_it);
            }

            for (auto source_it = subgraphs.begin(); source_it != subgraphs.end(); ++source_it) {
              for (auto target_it = model[*source_it].first.begin(); target_it != model[*source_it].first.end(); ++target_it) {
                if (*source_it <= *target_it) {
                  continue;
                }
                if (std::find(subgraphs.begin(), subgraphs.end(), *target_it) != subgraphs.end()) {
                  graph_instance.add_undirected_edge(*source_it, *target_it);
                }
              }
            }

            fingerprint::fingerprint fp = calculator.compute(graph_instance, model);
            OUTPUT->write(model_it->first, model_it->second.second, fp);
          }
        }
      }

      public:
      FingerprintWriter(std::string index, std::string subgraphs_file) : INDEX(index), SUBGRAPHS(subgraphs_file) {
        if (!INDEX.reset()) {
          throw common::exception::TitledException("The index file is empty or it is not possible to read it.");
        }
      }
      ~FingerprintWriter() {
        delete OUTPUT;
      }

      void add_features(std::string file)  {
        FEATURES.push_back(FeaturesReader(file));
        for (size_t i = 0; i < FEATURES.back().size(); i++) {
          auto ins = HEADERS_MAP.insert({FEATURES.back().header(i), HEADERS.size()});
          if (!ins.second) {
            throw common::exception::TitledException("Multiple columns contain the same header: '" + FEATURES.back().header(i));
          }
          HEADERS.push_back(FEATURES.back().header(i));
        }
      }

      void process(std::string settings_file, std::string edges_file, std::string output_path, FingerprintFormat format) {
        fingerprint::configuration calculator_configuration;
        calculator_configuration.load(settings_file, HEADERS_MAP);
        switch (format) {
          case FingerprintFormat::Binary:
            OUTPUT = new KnowledgebaseStream(output_path, HEADERS, calculator_configuration.size);
            break;
          case FingerprintFormat::Text:
            OUTPUT = new QueryStream(output_path, HEADERS);
            break;
          default:
            throw common::exception::TitledException("Unexpected output format.");
            break;
        }
        fingerprint::fingerprint_calculator calculator(calculator_configuration);

        EdgesReader edges(edges_file);

        std::string protein_id = INDEX.protein();
        std::string model_id = INDEX.model();

        // <id, <edges, features> >
        // NOTE: consider files are ordered and grouped by protein and biomolecule/model
        // NOTE: Definition of column names for the vector of strings is in the variable HEADERS and HEADERS_MAP.
        //       This splitted representation is 3x space saving and 15% quicker than single map<int, map<string, string> >.
        // NOTE: Actually there is no space different between map and unordered_map, but map is 10% quicker.
        std::map<int, std::pair<std::vector<int>, std::vector<std::string> > > model;

        do {
          if (protein_id != INDEX.protein()) {
            if (protein_id.size() > 0) {
              process(model, calculator);
              model.clear();
            }
            std::cout << protein_id << '\r';
            protein_id = INDEX.protein();
            model_id = INDEX.model();
          } else if (model_id != INDEX.model()) {
            process(model, calculator);
            model.clear();
            model_id = INDEX.model();
          }

          std::pair<std::map<int, std::pair<std::vector<int>, std::vector<std::string> > >::iterator, bool> node;
          if (edges.next(INDEX.index())) {
            node = model.insert({(int)INDEX.index(), {edges.value(), std::vector<std::string>()}});
          } else {
            std::cerr << "Residue no." << INDEX.index() << " is not presented in the edges file" << std::endl;
            node = model.insert({(int)INDEX.index(),{std::vector<int>(), std::vector<std::string>()}});
          }
          for (size_t i = 0; i < FEATURES.size(); i++) {
            if (FEATURES[i].next_line(INDEX.index())) {
              for (size_t j = 0; j < FEATURES[i].size(); j++) {
                node.first->second.second.push_back(FEATURES[i].value(j));
              }
            } else {
              std::cerr << "Residue no." << INDEX.index() << " is not presented in the features file no " << i+1 << std::endl;
              for (size_t j = 0; j < FEATURES[i].size(); j++) {
                node.first->second.second.push_back("");
              }
            }
          }
        } while (INDEX.next());
        process(model, calculator);
        OUTPUT->finalize();
      }
    };
  }
}