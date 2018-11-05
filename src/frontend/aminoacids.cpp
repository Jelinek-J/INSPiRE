// aminoacids.cpp : Defines the entry point for the console application.
//
//#define TESTING

#include "../elemental/string.h"
#include "../elemental/filesystem.h"
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include <set>

namespace inspire {
  namespace frontend {
    class Aminoacids {
      public:
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

      inline static std::string validate(const std::string &value) {
        std::string ret = elemental::string::trim(value);
        if (ret == "?") {
          return "";
        }
        if (ret.size() > 2 && ret[0] == '"' && ret[ret.size()-1] == '"') {
          return ret.substr(1, ret.size()-2);
        }
        return ret;
      }

      static void parse(std::string input_path, std::string output_path) {
        std::ifstream input(input_path);
        if (elemental::filesystem::is_directory(output_path) || output_path.size() == 0) {
          if (output_path.size() > 0 && output_path[output_path.size()-1] != elemental::filesystem::directory_separator) {
            output_path.push_back(elemental::filesystem::directory_separator);
          }
          output_path += "aminoacid.moc";
        } else if (!elemental::string::ends_with(output_path, ".moc")) {
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
          if (elemental::string::starts_with(line, "data_")) {
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
          } else if (elemental::string::starts_with(line, "_chem_comp.name")) {
            name = validate(line.substr(15));
          } else if (elemental::string::starts_with(line, "_chem_comp.mon_nstd_parent_comp_id")) {
            parent = validate(line.substr(34));
            parent = elemental::string::to_upper(parent);
          } else if (elemental::string::starts_with(line, "_chem_comp.pdbx_synonyms")) {
            synonyms = validate(line.substr(24));
          } else if (elemental::string::starts_with(line, "_chem_comp.one_letter_code")) {
            std::string tmp = elemental::string::trim(line.substr(26));
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
  }
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 3;
  const char* arg[] = {argv[0], "C:\\Users\\jan.jelinek\\Desktop\\components.cif", "C:\\Inspire\\"};
  argv = arg;
#endif // TESTING

  if (argc != 3) {
    if (argc != 0) {
      std::cerr << "Invalid number of arguments.\n\n";
    }
    std::cout << "Help\n\n";

    std::cout << "Extract parent aminoacid codes from PDB's cif file and store them in a file.\n";
    std::cout << "<cif_file> <output_path\t<cif_file> path to the input file and <output_path> set where to store output file\n";
    std::cout << "                       \tIf <output_path> is a directory or ends with a directory separator, 'aminoacid.moc' is used as the file name;\n";
    std::cout << "                       \tif the path is not a directory but does not but not ends with '.moc', the extension is appended.\n\n";
    return 0;
  }

  inspire::frontend::Aminoacids::parse(argv[1], argv[2]);
  
  return 0;
}

