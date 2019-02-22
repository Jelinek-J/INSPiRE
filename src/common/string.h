#pragma once

#include <string>
#include <algorithm>

namespace common {
  namespace string {
    // Check, whether string 'str' starts with 'prefix'
    inline bool starts_with(const std::string &str, const std::string &prefix) {
      return !str.compare(0, prefix.size(), prefix);
    }

    inline bool contains(const std::string &str, char ch) {
      return str.find(ch) != str.npos;
    }

    // Check, whether string 'str' contains 'infix' at position 'at'
    inline bool contains_at(const std::string &str, const std::string &infix, const size_t &at) {
      if (str.size() < at + infix.size()) {
        return false;
      }
      return !str.compare(at, infix.size(), infix);
    }

    // Check, whether string 'str' ends with 'suffix'
    inline bool ends_with(const std::string &str, const std::string &suffix) {
      if (str.size() < suffix.size()) {
        return false;
      }
      return !str.compare(str.size()-suffix.size(), str.npos, suffix);
    }

    // Check, whether string 'str' ends with 'suffix'
    inline bool ends_with(const std::string &str, const char &suffix) {
      return !str.empty() && str.back() == suffix;
    }

    // Remove spaces from both begin and end
    inline std::string trim(const std::string &str) {
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
    inline std::string to_upper(std::string &str) {
      transform(str.begin(), str.end(), str.begin(), toupper);
      return str;
    }

    // Replace all occurances
    inline std::string replace_all(std::string str, const char find, const std::string replace) {
      for (size_t i = str.size(); i = str.rfind(find, i--) != str.npos; ) {
        str = str.replace(i, 1, replace);
      }
      return str;
    }
  }
}
