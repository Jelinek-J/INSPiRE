#pragma once

#include "../common/filesystem.h"
#include "../common/exception.h"
#include <iostream>
#include <sstream>

namespace inspire {
  namespace backend {
    class Merger {
      private:
      static std::string basename(std::string &name) {
        size_t i = name.rfind(common::filesystem::directory_separator);
        std::string tmp = (i == name.npos ? name : name.substr(i+1));
        if (common::string::ends_with(tmp, ".fit")) {
          return tmp.substr(0, tmp.size()-4);
        } else {
          return tmp;
        }
      }

      static void merge_queries(std::string first, std::string second, std::string output) {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          output += basename(first);
          output += '-';
          output += basename(second);
          output += ".fit";
        } else if (!common::string::ends_with(output, ".fit")) {
          output += ".fit";
        }

        std::ifstream first_stream(first);
        std::ifstream second_stream(second);
        std::ofstream output_stream(output);
        std::string first_line;
        std::string second_line;
        while (std::getline(first_stream, first_line) && std::getline(second_stream, second_line)) {
          size_t ft = first_line.rfind('\t');
          size_t st = second_line.rfind('\t');
          if (ft != st || first_line.substr(0, ft) != second_line.substr(0, st)) {
            throw common::exception::TitledException("The fingerprint files have not the same keys");
          }
          output_stream << first_line << second_line.substr(st+1) << '\n';
        }
        output_stream.flush();
        output_stream.close();
        first_stream.close();
        second_stream.close();
      }

      static void merge_knowledgebase(std::string first, std::string second, std::string output) {
        if (!output.empty() && output.back() != common::filesystem::directory_separator) {
          output += common::filesystem::directory_separator;
        }
        std::string first_name = basename(first);
        std::string second_name = basename(second);
        if (first_name != second_name) {
          throw common::exception::TitledException("The files have a different names");
        }
        output += first_name;

        std::ifstream first_stream(first, std::ios::in | std::ios::binary);
        std::ifstream second_stream(second, std::ios::in | std::ios::binary);
        std::ofstream output_stream(output, std::ios::out | std::ios::binary);
        uint32_t first_length;
        first_stream.read(reinterpret_cast<char *>(&first_length), sizeof(first_length));
        uint32_t second_length;
        second_stream.read(reinterpret_cast<char *>(&second_length), sizeof(second_length));
        uint32_t length = first_length + second_length;
        output_stream.write(reinterpret_cast<const char*>(&length), sizeof(length));
        while (first_stream.peek() != EOF && second_stream.peek() != EOF) {
          uint32_t id;
          first_stream.read(reinterpret_cast<char *>(&id), sizeof(id));
          {
            uint32_t control_id;
            second_stream.read(reinterpret_cast<char *>(&control_id), sizeof(control_id));
            if (id != control_id) {
              throw common::exception::TitledException("The fingerprint files does not contains the same residues in the same order");
            }
          }
          output_stream.write(reinterpret_cast<const char*>(&id), sizeof(id));
          for (size_t i = 0; i < first_length; i++) {
            output_stream.put(first_stream.get());
          }
          for (size_t i = 0; i < second_length; i++) {
            output_stream.put(second_stream.get());
          }
        }
        output_stream.flush();
        output_stream.close();
        first_stream.close();
        second_stream.close();
      }

      static void merge_item(std::string first, std::string second, std::string output, std::string last) {
        if (common::string::ends_with(first, ".fit")) {
          merge_queries(first, second, output+last);
        } else if (common::string::ends_with(first, ".fin")) {
          merge_knowledgebase(first, second, output);
        } else {
          try {
            merge_queries(first, second, output+last);
          } catch (const std::exception&) {
            try {
              merge_knowledgebase(first, second, output);
            } catch (const std::exception&) {
              throw common::exception::TitledException("Unknown file format, it is not possible to merge files '" + first + "' and '" + second + "'");
            }
          }
        }
      }

      public:
      static void merge(std::string first, std::string second, std::string output) {
        if (!common::filesystem::exists(first)) {
          throw common::exception::TitledException(first + " does not exist");
        }
        if (!common::filesystem::exists(second)) {
          throw common::exception::TitledException(second + " does not exist");
        }
        if (common::filesystem::is_regular_file(first)) {
          if (common::filesystem::is_regular_file(second)) {
            merge_item(first, second, output, "");
          } else {
            throw common::exception::TitledException("The first file is a regular file while the second one not");
          }
        } else if (common::filesystem::is_directory(first)) {
          if (common::filesystem::is_directory(second)) {
            if (first.size() > 0 && first.back() != common::filesystem::directory_separator) {
              first += common::filesystem::directory_separator;
            }
            if (second.size() > 0 && second.back() != common::filesystem::directory_separator) {
              second += common::filesystem::directory_separator;
            }
            if (!common::filesystem::is_directory(output) && common::filesystem::exists(output)) {
              throw common::exception::TitledException(output + " is not a directory, but " + first + " and " + second + " are");
            }
            if (output.size() > 0 && output.back() != common::filesystem::directory_separator) {
              output += common::filesystem::directory_separator;
            }

            common::filesystem::RecursiveDirectoryFileIterator file_iterator(first);
            if (file_iterator.has_file()) {
              do {
                std::string relative = common::filesystem::relative(file_iterator.filename(), first);
                std::string counterpart = second + relative;
                if (!common::filesystem::exists(counterpart)) {
                  std::cerr << file_iterator.filename() << " does not have a counterpart in " << second << std::endl;
                } else {
                  std::string path = relative;
                  size_t last = path.rfind(common::filesystem::directory_separator);
                  if (last != path.npos) {
                    path = path.substr(0, last+1);
                    path = output + path;
                  } else {
                    path = output;
                  }
                  common::filesystem::create_directory_recursive(path);
                  merge_item(file_iterator.filename(), counterpart, path, relative.substr(last+1));
                }
              } while (file_iterator.has_next());
            }
          } else {
            throw common::exception::TitledException("The first file is a directory while the second one not");
          }
        } else {
          throw common::exception::TitledException("The first file is not regular file nor a directory");
        }
      }
    };
  }
}