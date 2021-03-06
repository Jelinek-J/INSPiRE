// mine.cpp : Defines the entry point for the console application.
//

#include "../common/filesystem.h"
#include "../common/string.h"
#include "../common/exception.h"
#include "../backend/mine.h"
#include <iostream>
#include <fstream>
#include <sstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "For each query fingerprint finds indices of k-most similar fingerprints with the same features of the central residue in the knowledge-base.\n\n";

  std::cout << "Usage:\t([-t <THREADS>] [-n <COUNT>] [-C] [-c <CENTRAL-FEATURES>] -k <KNOWLEDGE-BASE> ([-S] [-s <SIBLINGS-FILE>] [-p] <QUERY-FILE> <OUTPUT-PATH>)+)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t-t <THREADS>           \tNumber of threads that should be used for data mining. Default value is 1.\n";
  std::cout << "        \t-n <COUNT>             \tNumber of the most similar elements that will be returned.\n";
  std::cout << "        \t                       \tIf multiple fingerprints have the same similarity score, all fingerprints with the score equal to the score of the <COUNT>-th most similar element are returned too.\n";
  std::cout << "        \t                       \tDefault value is 1.\n";
  std::cout << "        \t-c <CENTRAL-FEATURES>  \tWhat features of central residues will be used for prefiltering of knowledge-base.\n";
  std::cout << "        \t                       \tMultiple features must be separated by a directory separator.\n";
  std::cout << "        \t-C                     \tClean filtering previously set with '-c' switch.\n";
  std::cout << "        \t-k <KNOWLEDGE-BASE>    \tPath to the root directory of a knowledge-base.\n";
  std::cout << "        \t-s <SIBLINGS-FILE>     \tPath to a file with defining, what knowledge-base's fingerprints should be skipped when mining most similar fingerprints for individual query fingerprints.\n";
  std::cout << "        \t                       \tThis is usefull for benchmarking to exclude fingerprints from the same protein/ benchmark instead of construction of new knowledge-base for each benchmark.\n";
  std::cout << "        \t-S                     \tClean siblings previously set with '-s' switch.\n";
  std::cout << "        \t<QUERY-FILE>           \tA path to a file with query fingerprints for which the most similar elements should be find.\n";
  std::cout << "        \t-p                     \tA switcher saying that the next one argument will be a <QUERY-FILE>.\n";
  std::cout << "        \t                       \tThis switcher is mandatory if a <QUERY-FILE> starts with a hyphen-minus sign.\n";
  std::cout << "        \t<OUTPUT-PATH>          \tWhere to store output file.\n";
  std::cout << "        \t                       \tIf <OUTPUT-PATH> is empty or ends with a directory separator, <QUERY-FILE>'s basename is used as the file name with '.med' as an extension.\n";
  std::cout << "        \t                       \tIf <OUTPUT-PATH> does not end with '.med' extension, the extension is appended.\n";
  std::cout << "        \t-h                     \tShow informations about the program\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 13;
  const char* args[] = {argv[0],
    "-t", "1",
    "-n", "1",
    "-c", "aminoacid",
    "-k", "C:\\Inspire\\precompiled\\fingerprints\\",
    "-s", "C:\\Inspire\\gvin\\related.pan",
    "C:\\Inspire\\test\\inspire\\fingerprints.fit", "C:\\Inspire\\test\\mint"
  };
  argv = args;
#endif // TESTING
  // NOTE: In the case of default (cs-CZ?) it writes '�' instead of (non-breakable?) space.
//  std::cout.imbue(std::locale("en_UK"));
  // Special case, the later initialization is useless
  if (argc <= 1 || common::string::starts_with(argv[1], "-h")) {
    help();
    return 0;
  }

  inspire::backend::Mine* mine = nullptr;
  try {
    std::set<std::string> filters;
    int threads = 1;
    int limit = 1;

    for (int i = 1; i < argc; i++) {
      if (std::strlen(argv[i]) == 2 && argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'C':
            filters.clear();
            break;
          case 'c':
            if (++i >= argc) {
              std::cerr << "Error: Considered features of central nodes are not specified";
              help();
              return 4564;
            }
            {
              filters.clear();
              std::stringstream arg(argv[i]);
              std::string part;
              while (std::getline(arg, part, common::filesystem::directory_separator)) {
                filters.insert(part);
              }
            }
            break;
          case 'h':
            help();
            break;
          case 'k':
            if (++i >= argc) {
              std::cerr << "Error: Knowledge-base file is not specified";
              help();
              return 7968;
            }
            if (mine != nullptr) {
              delete mine;
            }
            mine = new inspire::backend::Mine(argv[i], filters, threads, limit);
            // NOTE: For the case that the multiple predictions will be specified the next prediction will be without filtering.
            break;
          case 'n':
            if (++i >= argc) {
              std::cerr << "Error: Number of threads is not specified";
              help();
              return 3215;
            }
            limit = std::stoi(argv[i]);
            if (mine != nullptr) {
              mine->limit(limit);
            }
            break;
          case 'p': // Just for the case that the filename start as any switch
            if (++i >= argc) {
              std::cerr << "Error: Input file is not specified";
              help();
              return 4568;
            }
            if (++i >= argc) {
              std::cerr << "Error: Output file identifier is not specified";
              help();
              return 1542;
            }
            if (mine == nullptr) {
              std::cerr << "Error: Knowledge base is not specified yet";
              help();
              return 7912;
            }
            mine->select(argv[i-1], argv[i]);
            break;
          case 'S':
            if (mine != nullptr) {
              mine->clear_excludes();
            }
            break;
          case 's':
            if (++i >= argc) {
              std::cerr << "Error: File with definition of excludes is not specified";
              help();
              return 1553;
            }
            if (mine == nullptr) {
              std::cerr << "Error: Knowledge base is not specified yet";
              help();
              return 7913;
            }
            mine->load_excludes(argv[i]);
            break;
          case 't':
            if (++i >= argc) {
              std::cerr << "Error: Number of threads is not specified";
              help();
              return 8975;
            }
            threads = std::stoi(argv[i]);
            if (mine != nullptr) {
              mine->threads(threads);
            }
            break;
          default: // Filename arbitrary starts with "-"
            if (++i >= argc) {
              std::cerr << "Error: Output file identifier is not specified";
              help();
              return 1542;
            }
            if (mine == nullptr) {
              std::cerr << "Error: Knowledge base is not specified yet";
              help();
              return 7914;
            }
            mine->select(argv[i - 1], argv[i]);
            break;
        }
      } else {
        if (++i >= argc) {
          std::cerr << "Error: Output file identifier is not specified";
          help();
          return 1542;
        }
        if (mine == nullptr) {
          std::cerr << "Error: Knowledge base is not specified yet";
          help();
          return 7915;
        }
        mine->select(argv[i - 1], argv[i]);
      }
    }
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
#ifdef TESTING
    log << "ERROR: " << e.what() << std::endl;
#endif // TESTING
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
    return 2;
#ifdef TESTING
    log << "ERROR: " << e.what() << std::endl;
#endif // TESTING
  } catch (...) {
    std::cerr << "UNKNOWN ERROR" << std::endl;
    help();
    return 3;
#ifdef TESTING
    log << "UNKNOWN ERROR" << std::endl;
#endif // TESTING
  }
  if (mine != nullptr) {
    delete mine;
  }

  return 0;
}

