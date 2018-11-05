#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>

namespace elemental {
  namespace xml {
    class Xml {
      private:
      boost::property_tree::ptree root;

      Xml(boost::property_tree::ptree root) : root(root) { }

      public:
      Xml(std::istream &stream) {
        boost::property_tree::read_xml(stream, root);
      }

      bool get_value(const std::string &path, std::string &value) {
        auto node = root.get_child_optional(path);
        if (node.is_initialized()) {
          value = node.get().data();
          return true;
        } else {
          return false;
        }
      }

      bool get_attribute(const std::string &key, std::string &value) {
        return get_value("<xmlattr>." + key, value);
      }

      bool get_attribute(const std::string &path, const std::string &key, std::string &value) {
        return get_value(path + ".<xmlattr>." + key, value);
      }

      bool get_nodes(const std::string &path, const std::string &key, std::vector<Xml> &nodes) {
        auto level = root.get_child_optional(path);
        if (level.is_initialized()) {
          nodes.clear();
          for (auto level_it = level->begin(); level_it != level->end(); ++level_it) {
            if (level_it->first == key) {
              nodes.push_back(level_it->second);
            }
          }
          return true;
        } else {
          return false;
        }
      }
    };

  }
}