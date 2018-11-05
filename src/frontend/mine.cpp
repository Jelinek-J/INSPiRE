// mine.cpp : Defines the entry point for the console application.
//

#include "../elemental/filesystem.h"
#include "../elemental/string.h"
#include "../elemental/exception.h"
#include "../backend/mine.h"
#include <iostream>
#include <fstream>
#include <sstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "([-t <threads>] [-n <nearest elements>] [-C] [-c <central features>] -k <knowledge-base> ([-S] [-s <siblings>] [-p] <input> <output>)+)+\n";
  std::cout << "-t <threads>\tNumber of threads that should be used for data mining. Default value is 1.\n";
  std::cout << "-n <nearest elements>\tNumber of the most similar elements that will be returned. Default value is 1.\n";
  std::cout << "-CtClean filtering previously set with '-c' switch.\n";
  std::cout << "-c <central features>\tFeatures of central nodes used for filtering of knowledge-base separated by a directory separator.\n";
  // TODO: Path to file with binning.
  std::cout << "-k <knowledge-base>\tPath to the root directory of knowledge-base.\n";
  std::cout << "-S\tClean siblings previously set with '-s' switch.\n";
  std::cout << "-s <siblings>\tPath to the file with a specification of siblings that should be skipped.\n";
  std::cout << "[-p] <input>\tA path to A file with elements that should be predicted\n            \t-p switcher is necessary in the case that the filename starts with '-'.\n\n";
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
    "-k", "C:\\Inspire\\test\\fingerprints",
    "-s", "C:\\Inspire\\test\\siblings.exc",
    "C:\\Inspire\\test\\fingerprints.fit", "C:\\Inspire\\test\\mined"
  };
  argv = args;
#endif // TESTING
  // NOTE: In the case of default (cs-CZ?) it writes 'á' instead of (non-breakable?) space.
//  std::cout.imbue(std::locale("en_UK"));
  // Special case, the later initialization is useless
  if (argc <= 1 || elemental::string::starts_with(argv[1], "-h")) {
    help();
    return 0;
  }

  try {
    inspire::backend::Mine* mine = nullptr;
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
              while (std::getline(arg, part, elemental::filesystem::directory_separator)) {
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
  } catch (const elemental::exception::TitledException& e) {
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

  return 0;
}

