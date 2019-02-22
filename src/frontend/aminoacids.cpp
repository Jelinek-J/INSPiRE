// aminoacids.cpp : Defines the entry point for the console application.
//

#include "../common/string.h"
#include "../common/filesystem.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>

//#define TESTING

namespace inspire {
  namespace frontend {
    class Cif {
      protected:
      inline static std::string validate(const std::string &value) {
        std::string ret = common::string::trim(value);
        if (ret == "?") {
          return "";
        }
        if (ret.size() > 2 && ret[0] == '"' && ret[ret.size()-1] == '"') {
          return ret.substr(1, ret.size()-2);
        }
        return ret;
      }
    };

    class Parents : public Cif {
      private:
      inline static void insert(std::map<std::string, char> &dictionary,
                                const std::string &code, const std::string &name, const std::string &synonyms, const std::string &parent, const char &one) {
        if (code.size() > 0 && parent.size() > 0) {
          auto it = dictionary.find(parent);
          if (it == dictionary.end()) {
            std::cerr << "Component '" << code << "' has ";
            if (parent.find(',') == parent.npos) {
              std::cerr << "unknown parent aminoacid code";
            } else {
              std::cerr << "ambiguous parent aminoacid";
            }
            std::cerr << ": '" << parent << "'" << std::endl;
          } else {
            dictionary[code] = it->second;
          }
        }
      }

      public:
      static void parse_parents(std::string input_path, std::string output_path) {
        std::ifstream input(input_path);
        if (common::filesystem::is_directory(output_path) || output_path.size() == 0) {
          if (output_path.size() > 0 && output_path[output_path.size()-1] != common::filesystem::directory_separator) {
            output_path.push_back(common::filesystem::directory_separator);
          }
          output_path += "aminoacid.moc";
        } else if (!common::string::ends_with(output_path, ".moc")) {
          output_path += ".moc";
        }
        std::ofstream output(output_path);

        std::map<std::string, char> dictionary;
        dictionary["ALA"] = 'A';
        dictionary["ASX"] = 'B';
        dictionary["CYS"] = 'C';
        dictionary["ASP"] = 'D';
        dictionary["GLU"] = 'E';
        dictionary["PHE"] = 'F';
        dictionary["GLY"] = 'G';
        dictionary["HIS"] = 'H';
        dictionary["ILE"] = 'I';
        dictionary["XLE"] = 'J';
        dictionary["LYS"] = 'K';
        dictionary["LEU"] = 'L';
        dictionary["MET"] = 'M';
        dictionary["ASN"] = 'N';
        dictionary["PYL"] = 'O';
        dictionary["PRO"] = 'P';
        dictionary["GLN"] = 'Q';
        dictionary["ARG"] = 'R';
        dictionary["SER"] = 'S';
        dictionary["THR"] = 'T';
        dictionary["SEC"] = 'U';
        dictionary["VAL"] = 'V';
        dictionary["TRP"] = 'W';
        dictionary["UNK"] = 'X';
        dictionary["TYR"] = 'Y';
        dictionary["GLX"] = 'Z';
        dictionary["DA"] = '0';
        dictionary["DC"] = '1';
        dictionary["DG"] = '2';
        dictionary["DT"] = '3';
        dictionary["DU"] = '4';
        dictionary["A"] = '5';
        dictionary["C"] = '6';
        dictionary["G"] = '7';
        dictionary["T"] = '8';
        dictionary["U"] = '9';


        std::string code;
        std::string name;
        std::string synonyms;
        std::string parent;
        char one;

        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          if (common::string::starts_with(line, "data_")) {
            insert(dictionary, code, name, synonyms, parent, one);
            if (line.size() <= 8) {
              code = line.substr(5);
            } else {
              code = "";
              std::cerr << "Unexpected length of a item code: '" << line << std::endl;
            }
            name = "";
            synonyms = "";
            parent = "";
            one = ' ';
          } else if (common::string::starts_with(line, "_chem_comp.name")) {
            name = validate(line.substr(15));
          } else if (common::string::starts_with(line, "_chem_comp.mon_nstd_parent_comp_id")) {
            parent = validate(line.substr(34));
            parent = common::string::to_upper(parent);
          } else if (common::string::starts_with(line, "_chem_comp.pdbx_synonyms")) {
            synonyms = validate(line.substr(24));
          } else if (common::string::starts_with(line, "_chem_comp.one_letter_code")) {
            std::string tmp = common::string::trim(line.substr(26));
            if (tmp.size() != 1) {
              std::cerr << "Ambiguous one letter code: '" << tmp << "'" << std::endl;
            } else if (tmp == "?") {
              one = ' ';
            } else {
              one = tmp[0];
            }
          }
        }
        insert(dictionary, code, name, synonyms, parent, one);
        input.close();

        for (auto it = dictionary.begin(); it != dictionary.end(); ++it) {
          output << it->first << '\t' << it->second << '\n';
        }
        output.flush();
        output.close();
      }
    };

    class Compositions : public Cif {
      public:
      static void parse_composition(std::string input_path, std::string output_path) {
        std::ifstream input(input_path);
        if (common::filesystem::is_directory(output_path) || output_path.size() == 0) {
          if (output_path.size() > 0 && output_path[output_path.size()-1] != common::filesystem::directory_separator) {
            output_path.push_back(common::filesystem::directory_separator);
          }
          output_path += "composition.cit";
        } else if (!common::string::ends_with(output_path, ".cit")) {
          output_path += ".cit";
        }
        std::ofstream output(output_path);

        std::map<std::string, std::vector<std::string>> dictionary;
        std::string code;
        std::vector<std::string> compositions;
        std::string line;
        while (!input.eof() && std::getline(input, line)) {
          if (common::string::starts_with(line, "data_")) {
            if (line.size() <= 8) {
              code = line.substr(5);
            } else {
              code = "";
              std::cerr << "Unexpected length of a item code: '" << line << "'" << std::endl;
            }
            compositions.clear();
            compositions.push_back("");
          } else if (common::string::starts_with(line, "loop_")) {
            if (!input.eof() && std::getline(input, line)) {
              if (common::string::starts_with(line, "_chem_comp_atom.")) {
                size_t id = -1;
                size_t alt = -1;
                size_t type = -1;
                size_t leave = -1;
                size_t i = 0;
                do {
                  if (common::string::contains_at(line, "atom_id", 16)) {
                    id = i;
                  } else if (common::string::contains_at(line, "alt_atom_id", 16)) {
                    alt = i;
                  } else if (common::string::contains_at(line, "type_symbol", 16)) {
                    type = i;
                  } else if (common::string::contains_at(line, "pdbx_leaving_atom_flag", 16)) {
                    leave = i;
                  }
                  ++i;
                } while (!input.eof() && std::getline(input, line) && common::string::starts_with(line, "_"));
                if (id == -1) {
                  std::cerr << "Component '" << code << "' does not contains 'atom_id' within a '_chem_comp_atom' section" << std::endl;
                  return;
                }
                if (alt == -1) {
                  std::cerr << "Component '" << code << "' does not contains 'alt_atom_id' within a '_chem_comp_atom' section" << std::endl;
                  return;
                }
                if (type == -1) {
                  std::cerr << "Component '" << code << "' does not contains 'type_symbol' within a '_chem_comp_atom' section" << std::endl;
                  return;
                }
                if (leave == -1) {
                  std::cerr << "Component '" << code << "' does not contains 'pdbx_leaving_atom_flag' within a '_chem_comp_atom' section" << std::endl;
                  return;
                }
                // TODO: Not optimal
                size_t max = id;
                max = std::max(max, alt);
                max = std::max(max, type);
                max = std::max(max, leave);
                ++max;
                std::vector<std::string> first;
                std::vector<std::string> second;
                do {
                  std::vector<std::string> items;
                  std::stringstream parts(line);
                  for (size_t i = 0; i < max; i++) {
                    std::string part;
                    while (std::getline(parts, part, ' ') && part.empty()) { }
                    items.push_back(validate(part));
                  }
                  if (items[type] != "H" && items[leave] != "Y") {
                    first.push_back(items[id]);
                    second.push_back(items[alt]);
                  }
                } while (!input.eof() && std::getline(input, line) && !common::string::starts_with(line, "#"));
                std::sort(first.begin(), first.end());
                std::sort(second.begin(), second.end());
                output << code << '\t';
                for (auto first_it = first.begin(); first_it != first.end(); ++first_it) {
                  if (first_it != first.begin()) {
                    output << ' ';
                  }
                  output << *first_it;
                }
                if (first != second) {
                  for (auto second_it = second.begin(); second_it != second.end(); ++second_it) {
                    output << (second_it == second.begin() ? '\t' : ' ');
                    output << *second_it;
                  }
                }
                output << '\n';
              }
            } else {
              std::cerr << "Unexpected end of file within a 'loop_' section";
              return;
            }
          }
        }
        output.flush();
        output.close();
      }
    };
  }
}

void help() {
  std::cout << "Help\n\n";

  std::cout << "Extract informations about residues from IUPAC's components file.\n\n";

  std::cout << "Usage:\t<CIF-FILE> ( -c<OUTPUT-PATH> | -p<OUTPUT-PATH> | [-a]<OUTPUT-PATH> )+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<CIF-FILE>         \tPath to the components file\n";
  std::cout << "        \t-h                 \tShow informations about the program\n";
  std::cout << "    Extracted informations:\n";
  std::cout << "        \t-c<OUTPUT-PATH>    \tExtracts chemical composition for each residue (hydrogens, deuteriums and atoms with leaving flag are skipped).\n";
  std::cout << "        \t                   \tIf <OUTPUT-PATH> is a directory or ends with a directory separator, 'composition.cit' is used as the file name;\n";
  std::cout << "        \t                   \tif the path is not a directory but does not but not ends with '.cit', the extension is appended.\n";
  std::cout << "        \t-p<OUTPUT-PATH>    \tExtracts parent aminacid single-letter codes where available.\n";
  std::cout << "        \t                   \tIf <OUTPUT-PATH> is a directory or ends with a directory separator, 'aminoacid.moc' is used as the file name;\n";
  std::cout << "        \t                   \tif the path is not a directory but does not but  not ends  with  '.moc', the extension is appended.\n";
  std::cout << "        \t[-a]<OUTPUT-PATH>  \tShortcut for '-c<OUTPUT-PATH> -p<OUTPUT-PATH>'.  -a is mandatory, if <OUTPUT-PATH> starts with a hyphen-minus sign.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 3;
  const char* arg[] = {argv[0], "C:\\Users\\jan.jelinek\\Desktop\\components_iupac.cif", "-cC:\\Inspire\\"};
  argv = arg;
#endif // TESTING

  if (argc < 3) {
    if (argc != 0) {
      std::cerr << "Invalid number of arguments.\n\n";
    }
    help();
    return 0;
  }

  for (size_t i = 2; i < argc; i++) {
    switch (strlen(argv[i])) {
      case 0:
        inspire::frontend::Parents::parse_parents(argv[1], "");
        inspire::frontend::Compositions::parse_composition(argv[1], "");
        break;
      case 1:
        inspire::frontend::Parents::parse_parents(argv[1], argv[2]);
        inspire::frontend::Compositions::parse_composition(argv[1], argv[2]);
        break;
      default:
        if (argv[i][0] == '-') {
          switch (argv[i][1]) {
            case 'c':
              inspire::frontend::Compositions::parse_composition(argv[1], std::string(argv[2]).substr(2));
              break;
            case 'p':
              inspire::frontend::Parents::parse_parents(argv[1], std::string(argv[2]).substr(2));
              break;
            case 'a':
              inspire::frontend::Parents::parse_parents(argv[1], std::string(argv[2]).substr(2));
              inspire::frontend::Compositions::parse_composition(argv[1], std::string(argv[2]).substr(2));
              break;
            default:
              std::cerr << "Unrecognized switcher '" << argv[i] << "'" << std::endl;
              break;
          }
        } else {
          inspire::frontend::Parents::parse_parents(argv[1], argv[2]);
          inspire::frontend::Compositions::parse_composition(argv[1], argv[2]);
        }
        break;
    }
  }
  
  return 0;
}

