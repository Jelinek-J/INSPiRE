// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/filter.h"
#include "../elemental/exception.h"
#include "../elemental/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "<orifinal_index> <filtered_index> (<original_features> <filtered_features>)+\tFilters feature(s) file(s) base on the corresponding index files.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 5;
  const char* args[] = {argv[0],
    "C:\\Inspire\\basic2\\gvin\\construction\\residues.ind",
    "C:\\Inspire\\basic2\\gvin\\residues.ind",
    "C:\\Inspire\\basic2\\gvin\\construction\\interfaces.tur",
    "C:\\Inspire\\basic2\\gvin\\interfaces2.tur"
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
  } catch (const elemental::exception::TitledException& e) {
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
