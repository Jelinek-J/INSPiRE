#pragma once

#include <string>
#include <algorithm>

namespace elemental {
  namespace string {
    // Check, whether string 'str' starts with 'prefix'
    bool starts_with(const std::string &str, const std::string &prefix) {
      return !str.compare(0, prefix.size(), prefix);
    }

    // Check, whether string 'str' contains 'infix' at position 'at'
    bool contains_at(const std::string &str, const std::string &infix, const size_t &at) {
      if (str.size() < at + infix.size()) {
        return false;
      }
      return !str.compare(at, infix.size(), infix);
    }

    // Check, whether string 'str' ends with 'suffix'
    bool ends_with(const std::string &str, const std::string &suffix) {
      if (str.size() < suffix.size()) {
        return false;
      }
      return !str.compare(str.size()-suffix.size(), str.npos, suffix);
    }

    // Remove spaces from both begin and end
    std::string trim(const std::string &str) {
      if (str.empty()) {
        return "";
      }
      size_t first = str.find_first_not_of(' ');
      size_t last = str.find_last_not_of(' ');
      if (first == str.npos || last == str.npos || first > last) {
        return "";
      }
      return str.substr(first, last-first+1);
    }

    // Convert string to upper
    std::string to_upper(std::string &str) {
      transform(str.begin(), str.end(), str.begin(), toupper);
      return str;
    }
  }
}
