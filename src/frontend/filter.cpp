// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/filter.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Filters feature(s) file(s) base on the corresponding index files. It is usefull e.g. for benchmarking or reusing of index-independent features.\n\n";

  std::cout << "Usage:\t<ORIGINAL-INDEX> <FILTERED-INDEX> (<ORIGINAL-FEATURES> <FILTERED-FEATURES>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<ORIGINAL-INDEX>     \tIndex file for feature file(s) to be filtered\n";
  std::cout << "        \t<FILTERED-INDEX>     \tIndex file by which to filter\n";
  std::cout << "        \t<ORIGINAL-FEATURES>  \tFeatures file that will be filtered\n";
  std::cout << "        \t<FILTERED-FEATURES>  \tWhere to store output file with filtered features.\n";
  std::cout << "        \t                     \tIf <FILTERED-FEATURES> is empty or ends with a directory separator, <ORIGINAL-FEATURES>'s basename is used as the file name;\n";
  std::cout << "        \t                     \tif <FILTERED-FEATURES> is not a directory but does not but not ends with '.tur' extension, the extension is appended.\n";
  std::cout << "        \t-h                   \tShow informations about the program\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 5;
  const char* args[] = {argv[0],
    "C:\\Inspire\\filter\\real\\residues.ind",
    "C:\\Inspire\\filter\\query\\residues.ind",
    "C:\\Inspire\\filter\\real\\interfaces.tur",
    "C:\\Inspire\\filter\\query\\interfaces.tur"
  };
  argv = args;
#endif // TESTING

  if (argc < 3 || argc %2 != 1) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  inspire::backend::Filter filter(argv[1], argv[2]);

  try {
    for (size_t i = 3; i < argc; i+=2) {
      filter.filter(argv[i], argv[i+1]);
    }
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
#ifdef TESTING
    log << "ERROR: " << e.what() << std::endl;
#endif // TESTING
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
#ifdef TESTING
    log << "ERROR: " << e.what() << std::endl;
#endif // TESTING
    return 2;
  } catch (...) {
    std::cerr << "UNKNOWN ERROR" << std::endl;
    help();
#ifdef TESTING
    log << "UNKNOWN ERROR" << std::endl;
#endif // TESTING
    return 3;
  }

  return 0;
}
