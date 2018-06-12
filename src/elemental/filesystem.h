#pragma once

#include <string>
#include "boost\filesystem.hpp"

// NOTE: If you do not like Boost or have problems with getting it to work, all usages of Boost are concentrated in namespace 'elemental'
// to make its replacement as easy as possible

namespace elemental {
  namespace filesystem {
    // Checks, whether the path exists.
    // boost::filesystem::exists, feel free to override it if you wants to aviod using Boost library
    inline bool exists(std::string path) {
      return boost::filesystem::exists(path);
    }

    // Checks, whether the path is a directory.
    // boost::filesystem::is_directory, feel free to override it if you wants to aviod using Boost library
    inline bool is_directory(std::string path) {
      return boost::filesystem::is_directory(path);
    }

    // Checks, whether the path is a regular file.
    // boost::filesystem::is_regular_file, feel free to override it if you wants to aviod using Boost library
    inline bool is_regular_file(std::string path) {
      return boost::filesystem::is_regular_file(boost::filesystem::path(path));
    }

    // Make a path relative to the given directory
    // boost::filesystem::relative, feel free to override it if you wants to aviod using Boost library
    inline std::string relative(std::string path, std::string directory) {
      return boost::filesystem::relative(boost::filesystem::path(path), boost::filesystem::path(directory)).string();
    }

    // Checks, whether the path exists.
    // boost::filesystem::exists, feel free to override it if you wants to aviod using Boost library
    char directory_separator = boost::filesystem::path::preferred_separator;

    // Create directory and all nonexisting ancestor directories
    inline bool create_directory_recursive(std::string path) {
      return boost::filesystem::create_directories(path);
    }

    // Check whether <path> is a portable directory name
    inline bool is_portable_directory(std::string name) {
      return boost::filesystem::portable_directory_name(name);
    }

    // Check whether <path> is a portable file name
    inline bool is_portable_file(std::string name) {
      return boost::filesystem::portable_file_name(name);
    }

    // Iterates through all files in a given directory and all its subdirectories
    class RecursiveDirectoryFileIterator {
      // End condition of the iterator
      const boost::filesystem::recursive_directory_iterator END;
      // Iterator of directory
      boost::filesystem::recursive_directory_iterator FILE;

      public:
      // Constructor accept a directory that should be iterated
      RecursiveDirectoryFileIterator(std::string directory) : FILE(directory) { }

      // Returns whether the iterator is set to a file
      bool has_file() {
        return FILE != END;
      }

      // Move to a next file if it exists and returns TRUE, else it returns false
      bool has_next() {
        while (++FILE != END) {
          if (boost::filesystem::is_regular_file(FILE->status())) {
            return true;
          }
        }
        return false;
      }

      // Returns a filename and path of the current file
      std::string filename() {
        return FILE->path().string();
      }

    };
  }
}