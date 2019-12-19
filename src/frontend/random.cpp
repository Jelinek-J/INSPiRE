// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../common/string.h"
#include "../common/random.h"
#include "../backend/index.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Extracts required features from proteins indexed in a given index file.\n\n";

  std::cout << "Usage:\t[-b|-c|-bc|-w] <INDEX-FILE> <OUTPUT-PATH> [-s] (-a<TRANSFORMATION-FILE>|-c|-e|-i<RADII-FILE>[<DISTANCE>]|-r<RADII-FILE>;<COMPOSITION-FILE>;<MAX-SASA-FILE>|-t)+ (<PROTEINS-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INDEX-FILE>     \tPath to a index file\n";
  std::cout << "        \t-s               \tSave each feature in a separate file (otherwise all features are stored together in one file)\n";
  std::cout << "        \t<OUTPUT-PATH>    \tWhere to store output file.\n";
  std::cout << "        \t                 \tIf '-s' is typed, <OUTPUT-PATH> must be a directory; and each feature is stored in a separated file named according to the corresponding feature with '.tur' as an extension.\n";
  std::cout << "        \t                 \tOtherwise(if '-s' is not typed), all features are stored in the same file.\n";
  std::cout << "        \t                 \t  If <OUTPUT-PATH> is a directory or ends with a directory separator, 'features.tur' is used as the file name.\n";
  std::cout << "        \t                 \t  If <OUTPUT-PATH> does not end with '.tur' extension, the extension is appended.\n";
  std::cout << "        \t<PROTEINS-PATH>  \tPath to a protein or a directory with proteins\n";
  std::cout << "        \t-h, --help       \tShow informations about the program\n";
  std::cout << "    Iterators: (if no iterator is specified, only the first biomolecule from the first model with the first crystallographic transformation is used)\n";
  std::cout << "        \t-b   \tAll biomolecules and models, but only the first crystallographic transformation are used\n";
  std::cout << "        \t-c   \tAll crystallographic transformations, but only the first biomolecule and model are used\n";
  std::cout << "        \t-bc  \tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "        \t-w   \tIgnore both biomolecules and crystallographic transformation, use all chains as they are\n";
  std::cout << "    Features:\n";
  std::cout << "        \t-a<TRANSFORMATION-FILE>\n";
  std::cout << "        \t    \tAminoacid type three-letter code is transformed with transformation defined in <TRANSFORMATION-FILE>, that should be in format 'key\\tvalue'\n";
  std::cout << "        \t-c  \tCoordinates of carbon_alpha of an aminoacid\n";
  std::cout << "        \t-e  \tAtomic composition of a residue (helium and deuterium are skipped)\n";
  std::cout << "        \t-i<RADII-FILE>[<DISTANCE>]\n";
  std::cout << "        \t    \tWhether a residue is an interfacial residue with <RADII-FILE> defining radiuses of chemical elements and\n";
  std::cout << "        \t    \t<DISTANCE> sets the maximal allowed distance of two van der Waals radiuses (0.5A is a default value).\n";
  std::cout << "        \t    \t<DISTANCE> must be separated by a space from <RADII-FILE>\n";
  std::cout << "        \t-r<RADII-FILE>;<COMPOSITION-FILE>;<MAX-SASA-FILE>\n";
  std::cout << "        \t    \tRelative solvent accessible surface area with residues' composition defined in <COMPOSITION-FILE>, atomic radiuses defined in <RADII-FILE> and\n";
  std::cout << "        \t    \treference solvent accessible surface areas defined in <MAX-SASA-FILE>.\n";
  std::cout << "        \t-t  \tTemperature factor of an aminoacid\n\n";

  std::cout << "Feature Files Format:\n";
  std::cout << "\tHeader line:\tThe first line of each feature file; names of columns (features) are separated by a tabulator.\n";
  std::cout << "\tData line:  \tEach feature is separated by a tabulator, if it is not possible to extract an feature, the corresponding cell is empty.\n";
  std::cout << "\t            \tDefaultly, index of each line is an index of the previous line plus one (the first line has index '1').\n";
  std::cout << "\t            \tIf a line has not consecutive index, the index is specified in an extra(last) column.\n\n";



  std::cout << "Help\n\n";

  std::cout << "Check what complexes are valid to be in a knowledge-base\n";
  std::cout << "<index_file> <index_log> <aminoacid_file> <aminoacid_log> <coordinates_file> <coordinates_log> (<pdb_directory>|<pdb_file>)+\t.\n";
  std::cout << "\tCheck, whether:\n";
  std::cout << "\t\t1) File contains REMARKS 350 (or equivalent in mmcif/ xml);\n";
  std::cout << "\t\t2) Is not a monomer;\n";
  std::cout << "\t\t3) Does not contain DNA, RNA, unknown amino acid (X, B, J, Z), modified aminoacid nor unknown residue;\n";
  std::cout << "\t\t4) Contains no HETATM except H2O;\n";
  std::cout << "\t\t5) No chain contains less than 20 aminoacids;\n";
  std::cout << "\t\t6) Contains less than 10 incomplete aminoacids;\n";
  std::cout << "\t\t7) No chain contains more than 1% of incomplete amino acids;\n";
  std::cout << "Each identifier (chain, complex) also applies to the following lines until it is redefined.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  argc = 3;
  const char* arg[] = { argv[0], "C:\\Inspire\\random\\residues.ind", "C:\\Inspire\\random\\mined.rty"
  };
  argv = arg;
#endif // TESTING
  if (argc < 3) {
    if (argc != 0) {
      std::cerr << "Unexpected number of arguments:";
      for (size_t i = 1; i < argc; i++) {
        std::cerr << "\n[" << i << "]:\t" << argv[i];
      }
      std::cerr << std::endl;
    }
    help();
    return 0;
  }

  try {
    std::map<std::string, std::set<std::string>> collisions;
    std::ifstream similarities(argv[2]);
    std::string line;
    while (std::getline(similarities, line)) {
      std::set<std::string> &id = collisions[line.substr(0, 4)];
      while (std::getline(similarities, line) && !line.empty()) {
        if (!common::string::starts_with(line, "0.0")) {
          size_t first = line.find('\t');
          if (first == line.npos) {
            throw common::exception::TitledException("Unexpected format of similarity file: '" + line + "'");
            return 1;
          }
          size_t second = line.rfind('\t');
          if (first == second) {
            throw common::exception::TitledException("Unexpected format of similarity file: '" + line + "'");
            return 2;
          }
          id.emplace(line.substr(first+1, second-first-1));
        }
      }
    }

    std::vector<std::string> order;
    std::map<std::string, char> proteins;
    inspire::backend::Index index(argv[1]);
    if (index.reset()) {
      std::string prev_protein;
      std::string prev_chain;
      std::vector<char> chains;
      bool add = false;
      do {
        std::string new_protein = index.protein();
        if (prev_protein != new_protein) {
          if (!chains.empty()) {
            proteins[prev_protein] = chains[common::random::random_size_t(chains.size())];
            order.push_back(prev_protein);
            chains.clear();
          }
          prev_protein = new_protein;
          prev_chain = "";
          add = collisions.find(new_protein) != collisions.end();
        }
        if (add && index.model().empty()) {
          std::string new_chain = index.chain();
          if (new_chain.size() == 1 && prev_chain != new_chain) {
            chains.push_back(new_chain[0]);
            prev_chain = new_chain;
          }
        }
      } while (index.next());
      if (!chains.empty()) {
        proteins[prev_protein] = chains[common::random::random_size_t(chains.size())];
        order.push_back(prev_protein);
      }
    }
    std::random_shuffle(order.begin(), order.end());

    std::set<std::string> selected;
    std::set<std::string> forbidden;
    while (selected.size() < 100 && selected.size() + forbidden.size() < order.size()) {
      std::string next = order[common::random::random_size_t(order.size())];
      if (selected.find(next) == selected.end() && forbidden.find(next) == forbidden.end()) {
        std::set<std::string> &collision = collisions[next];
        bool add = true;
        for (auto collision_it = collision.begin(); collision_it != collision.end(); ++collision_it) {
          if (selected.find(*collision_it) != selected.end()) {
            add = false;
            break;
          }
        }
        if (add) {
          selected.emplace(next);
          forbidden.insert(collision.begin(), collision.end());
        } else {
          forbidden.emplace(next);
        }
      }
    }
    if (selected.size() < 100) {
      std::cerr << "Not enough unrelated structures, only " << selected.size() << " structures taken" << std::endl;
    }
    for (auto selected_id = selected.begin(); selected_id != selected.end(); ++selected_id) {
      if (selected_id != selected.begin()) {
        std::cout << ", ";
      }
      std::cout << *selected_id << '.' << proteins[*selected_id];
    }
    std::cout << std::endl;
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
    return 2;
  } catch (...) {
    std::cerr << "UNKNOWN ERROR" << std::endl;
    help();
    return 3;
  }

  return 0;
}

