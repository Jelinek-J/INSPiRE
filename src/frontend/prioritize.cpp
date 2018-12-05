// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/prioritize.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "((-r[<weight>]) (-a[<weight>]) (-d[<weight>]) ([-] <mine> <output>)+)+\tPrioritize residues in <mine> file based on their statistics.\n";
  std::cout << "                                                                      \tWeights of individual attributes can be set.\n";
  std::cout << "                                                                      \t-r\tratio of positive and negative examples is used;\n";
  std::cout << "                                                                      \t-a\tabsolute number of positive examples is used;\n";
  std::cout << "                                                                      \t-d\tdistance of examples from the query is used.\n";
  std::cout << "                                                                      \tIf <weight> is not specified, it is set to '1'.\n";
  std::cout << "                                                                      \t-\tis required if a <mine> starts with '-'.";
  std::cout << "NOTE: Principially, <interfaces> could be an arbitrary features file.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-prioritize.log");
  argc = 6;
  const char* args[] = {argv[0],
    "-a10", "-r1", "-d0.1",
    "C:\\Inspire\\gvin\\ratios.sas",
    "C:\\Inspire\\gvin\\prioritized-a_r_d"
  };
  argv = args;
#endif // TESTING

  if (argc < 4) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  try {
    size_t argv_index = 1;
    inspire::backend::Prioritizer* prioritizer = nullptr;
    double ratio = 0;
    double absolute = 0;
    double distance = 0;
    while (argv_index < argc) {
      for (; argv_index < argc && strlen(argv[argv_index]) >= 2 && argv[argv_index][0] == '-'; ++argv_index) {
        switch (argv[argv_index][1]) {
          case 'r':
            if (prioritizer != nullptr) {
              delete prioritizer;
              absolute = 0;
              distance = 0;
            }
            if (strlen(argv[argv_index]) == 2) {
              ratio = 1;
            } else {
              ratio = std::stod(std::string(argv[argv_index]).substr(2));
            }
            break;
          case 'a':
            if (prioritizer != nullptr) {
              delete prioritizer;
              ratio = 0;
              distance = 0;
            }
            if (strlen(argv[argv_index]) == 2) {
              absolute = 1;
            } else {
              absolute = std::stod(std::string(argv[argv_index]).substr(2));
            }
            break;
          case 'd':
            if (prioritizer != nullptr) {
              delete prioritizer;
              ratio = 0;
              absolute = 0;
            }
            if (strlen(argv[argv_index]) == 2) {
              distance = 1;
            } else {
              distance = std::stod(std::string(argv[argv_index]).substr(2));
            }
            break;
          default:
            std::cerr << "Unknown specifier '" << argv[argv_index] << "'." << std::endl;
            help();
            return 6;
        }
      }
      prioritizer = new inspire::backend::LinearPrioritizer(ratio, absolute, distance);
      for (; argv_index < argc && (strlen(argv[argv_index]) < 1 || argv[argv_index][0] != '-' || (strlen(argv[argv_index]) == 1 && ++argv_index < argc)); ++argv_index) {
        std::string input(argv[argv_index++]);
        std::string output(argv[argv_index++]);
        prioritizer->prioritize(input, output);
      }
    }

    delete prioritizer;
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

  return 0;
}
