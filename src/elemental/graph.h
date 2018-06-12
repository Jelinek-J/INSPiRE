#pragma once

#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"


namespace elemental {
  namespace graph {
    typedef boost::property<boost::edge_weight_t, int> edge_weight_property;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, edge_weight_property> undirected_graph;
    typedef boost::exterior_vertex_property<undirected_graph, int> distance_property;
    typedef distance_property::matrix_type distance_matrix;
    typedef distance_property::matrix_map_type distance_matrix_map;
    typedef boost::constant_property_map<boost::graph_traits<undirected_graph>::edge_descriptor, int> const_property_map;
    typedef boost::adjacency_list<>::vertex_descriptor vertex;
    typedef boost::property_tree::ptree property_tree;

    vertex add_vertex(undirected_graph g) {
      return boost::add_vertex(g);
    }

    void add_edge(vertex v1, vertex v2, undirected_graph g) {
      boost::add_edge(v1, v2, g);
    }

    void shortest_paths(undirected_graph g, distance_matrix_map distances, const_property_map properties) {
      boost::floyd_warshall_all_pairs_shortest_paths(g, distances, boost::weight_map(properties));
    }

    void parse_json(const std::string &path, property_tree &root) {
      boost::property_tree::read_json(path, root);
    }
  }
}