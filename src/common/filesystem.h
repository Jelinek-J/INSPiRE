#pragma once

#include <string>
#include "string.h"
#include "boost/filesystem.hpp"

// NOTE: If you do not like Boost or have problems with getting it to work, all usages of Boost are concentrated in namespace 'common'
// to make its replacement as easy as possible

namespace common {
  namespace filesystem {
    // Checks, whether the path exists.
    // boost::filesystem::exists, feel free to override it if you wants to aviod using Boost library
    inline bool exists(std::string path) {
      return boost::filesystem::exists(path);
    }

    // Checks, whether the path is a directory.
    // boost::filesystem::is_directory, feel free to override it if you wants to aviod using Boost library
    inline bool is_directory(std::string path) {
      return path.empty() || boost::filesystem::is_directory(path);
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

    // Generate unique file/directory name 
    inline std::string unique_name() {
      return boost::filesystem::unique_path().string();
    }

    // Check whether directory name ends with a directory separator and append it if not
    inline std::string enclose_directory_name(std::string name) {
      if (name.size() == 0 || name.back() == common::filesystem::directory_separator) {
        return name;
      }
      return name + directory_separator;
    }

    // Copy file from one location to another one
    inline void copy(std::string from, std::string to) {
      boost::filesystem::copy(from, to);
    }

    // Move file from one location to another one
    inline void move(std::string from, std::string to) {
      boost::filesystem::rename(from, to);
    }

    // Delete directory and all its content
    inline void remove_recursively(std::string directory) {
      boost::filesystem::remove_all(directory);
    }

    inline void remove_file(std::string file) {
      boost::filesystem::remove(file);
    }

    // Check, whether <path> ends with file or directory name (separator); in the second case, add the <filename>.
    // Then check, whether path ends with <extension> - if not, add the extension.
    // NOTE: To consider path as a directory, it must end with a directory separator, otherwise an extension is added only.
    // NOTE: Dot is not added prior extension implicitly.
    inline std::string complete(const std::string &path, const std::string &filename, const std::string &extension) {
      if (!extension.empty() && string::ends_with(filename, extension)) {
        std::string filename2 = filename.substr(0, filename.size()-extension.size());
        return complete(path, filename2, extension);
      }

      if (path.empty()) {
        return filename + extension;
      }
      if (string::ends_with(path, directory_separator)) {
        if (!is_directory(path)) {
          create_directory_recursive(path);
        }
        return path + filename + extension;
      }
      /* See the upper NOTE
      if (is_directory(path)) {
        return path + std::string(1, directory_separator) + filename + extension;
      }*/
      if (string::contains(path, directory_separator)) {
        std::string dir = path.substr(0, path.rfind(directory_separator));
        if (!is_directory(dir)) {
          create_directory_recursive(dir);
        }
        if (dir.size()+1 == path.size()) {
          return path + filename + extension;
        }
      }
      if (string::ends_with(path, extension)) {
        return path;
      }
      return path + extension;
    }

    // Identify basename in <pattern> and then call complete(<path>, basename, <extension>)
    inline std::string copy_filename(const std::string &path, const std::string &pattern, const std::string &extension) {
      size_t separator = pattern.rfind(directory_separator);
      size_t dot = pattern.rfind('.');
      if (separator == pattern.npos) {
        return complete(path, dot == pattern.npos || extension.empty() ? pattern : pattern.substr(0, dot), extension);
      } else {
        return complete(path, (dot == pattern.npos || dot < separator || extension.empty()) ? pattern.substr(separator+1) : pattern.substr(separator+1, dot-separator-1), extension);
      }
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
