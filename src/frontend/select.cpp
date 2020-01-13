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

  std::cout << "Extracts required chains from a given index file. Note that only chains from the first model are preserved.\n\n";

  std::cout << "Usage:\t<INPUT-INDEX-FILE> <OUTPUT-INDEX-FILE> (<PROTEIN-ID>(.<CHAIN>)*)+";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INPUT-INDEX-FILE>     \tPath to a index file to be filtered.\n";
  std::cout << "        \t<OUTPUT-INDEX-FILE>    \tPath where to store filtered input.\n";
  std::cout << "        \t<PROTEINS-ID>.<CHAIN>  \tWhat chains should be preserved. There can be specified multiple chains separated by a dot '.'. If no chain is specified, all chains are preserved.\n";
  std::cout << "        \t-h, --help             \tShow informations about the program.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  argc = 103;
  const char* arg[] = { argv[0], "C:\\Inspire\\random\\residues.ind", "C:\\Inspire\\random\\benchmark.ind",
                        "5MT7.A", "5MT9.D", "5MUQ.H", "5MYL.A", "5MZZ.B", "5NDC.A", "5NDT.B", "5NJA.E", "5NSM.D", "5NWX.B", "5NYJ.J", "5NZP.B", "5NZT.A", "5O4H.A", "5O95.B", "5OJI.A", "5OL0.A", "5OUQ.A", "5OY3.A", "5UEE.B", "5UH8.I", "5UKY.B", "5USE.B", "5UU3.A", "5UZ0.B", "5UZI.B", "5V0Z.C", "5V3G.F", "5V3U.B", "5V4Y.A", "5V6X.A", "5VAX.E", "5VKY.A", "5VS0.D", "5VSV.B", "5VY9.F", "5VZ8.A", "5W4C.A", "5W5P.A", "5WAV.B", "5WCD.L", "5WEU.A", "5WKH.A", "5WWE.B", "5X06.F", "5X54.A", "5X8F.C", "5X9R.B", "5XB1.B", "5XBS.B", "5XHC.C", "5XI7.A", "5XNC.A", "5XOM.B", "5XTB.Q", "5XTO.A", "5XVC.L", "5Y0R.A", "5Y2U.B", "5Y6U.B", "5Y96.B", "5YAT.B", "5YBF.A", "5YEW.B", "5YHW.H", "5YQG.C", "5YS6.A", "5YVN.A", "5YY9.C", "5YZ0.D", "5YZL.B", "5YZN.A", "6AMA.N", "6AQT.B", "6AU5.A", "6AVZ.A", "6AZS.A", "6B21.A", "6B5F.A", "6B60.B", "6BB1.D", "6BB2.A", "6BDK.A", "6BHG.A", "6BIY.C", "6BN7.B", "6BNQ.H", "6BQ6.B", "6BVS.A", "6BWL.A", "6BWO.A", "6EN0.A", "6EP3.A", "6ET9.C", "6EVQ.A", "6EY5.B", "6EZG.A", "6F5N.B", "6FAN.E", "6FE4.E"
  };
  argv = arg;
#endif // TESTING
  if (argc < 4) {
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
    std::map<std::string, std::set<std::string>> select;
    for (size_t i = 3; i < argc; i++) {
      std::stringstream chains(argv[i]);
      std::string protein;
      if (!std::getline(chains, protein, '.')) {
        throw common::exception::TitledException("Unexpected error during parsing of chain(s) identifier: '" + std::string(argv[i]) +"'");
      }
      std::string chain;
      auto& subselect = select[protein];
      while (std::getline(chains, chain, '.')) {
        subselect.emplace(chain);
      }
    }

    inspire::backend::Index index(argv[1]);
    std::ofstream output(argv[2]);
    if (index.reset()) {
      std::string prev_protein;
      std::string prev_chain;
      auto select_it = select.end();
      bool add_chain = false;
      std::string suffix;
      do {
        std::string new_protein = index.protein();
        if (prev_protein != new_protein) {
          select_it = select.find(new_protein);
          bool add_protein = select_it != select.end();
          bool add_chain = false;
          prev_protein = new_protein;
          prev_chain = "";
          if (select_it != select.end()) {
            suffix = "\t\t" + new_protein;
          }
        }
        if (select_it != select.end()) {
          if (index.model().empty()) {
            std::string new_chain = index.chain();
            if (prev_chain != new_chain) {
              if (select_it->second.size() > 0 && select_it->second.find(new_chain) == select_it->second.end()) {
                add_chain = false;
                prev_chain = "";
              } else {
                add_chain = true;
                prev_chain = new_chain;
                suffix = "\t" + new_chain + suffix;
              }
            }
            if (add_chain) {
              output << index.aminoacid() << suffix << '\n';
              suffix = "";
            }
          } else {
            add_chain = false;
            prev_chain = "";
          }
        }
      } while (index.next());
    }
    output.flush();
    output.close();
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

