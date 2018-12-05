// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/optimize.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "(-f<precision (-l<interfaces> ([-] <stats> [-o<output>])+)+)+\tSelect the optimal parameters for the given predictor based on the given input.\n";
  std::cout << "                                                             \t-f<precision>\tlinear selector will be optimised iwth the given precission;\n";
  std::cout << "                                                             \t<thresholds> must contains number of 2-combination of <x> thresholds,\n";
  std::cout << "                                                             \twhere <x> is number of labels in <mine> and they are ordered from the earlier to later\n";
  std::cout << "                                                             \t-\tis required if a <mine> starts with '-'.\n";
  std::cout << "                                                             \t<stats>\tset used for optimization.\n";
  std::cout << "                                                             \t-l<interface>\treal labels.\n";
  std::cout << "                                                             \t-o<output>\tlocation to store output.\n";
  std::cout << "NOTE: Principially, <interfaces> could be an arbitrary features file.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-optimize.log");
  argc = 4;
  const char* args[] = {argv[0],
    "-lC:\\Inspire\\basic2\\gvin\\interfaces.tur",
    "C:\\Inspire\\basic2\\gvin2\\ratios.sas",
    "-oC:\\Inspire\\basic2\\gvin2\\dependence"
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
    std::string function = "f1000";
    inspire::backend::Optimizer* optimizer = nullptr;
    while (argv_index < argc && strlen(argv[argv_index]) >= 2 && argv[argv_index][0] == '-') {
      switch (argv[argv_index][1]) {
        case 'f':
          if (optimizer != nullptr) {
            delete optimizer;
          }
          function = std::string(argv[argv_index++]).substr(1);
          break;
        case 'l':
          if (optimizer != nullptr) {
            delete optimizer;
          }
          if (common::string::starts_with(function, "f")) {
            optimizer = new inspire::backend::FractionalOptimizer(std::string(argv[argv_index]).substr(2), std::stoi(function.substr(1)));
          }
          
          while (++argv_index < argc && (strlen(argv[argv_index]) < 2 || argv[argv_index][0] != '-')) {
            if (argv[argv_index] == std::string("-")) {
              ++argv_index;
            }
            if (argv_index+1 >= argc) {
              std::cerr << "Missing output specifier for input '" << argv[argv_index] << "'." << std::endl;
              help();
              return 6;
            }
            std::string input(argv[argv_index++]);
            if (common::string::starts_with(argv[argv_index], "-o")) {
              std::string output(argv[argv_index++]);
              optimizer->optimize(input, output.substr(2));
            } else {
              optimizer->optimize(input);
            }
          }
          break;
        default:
          std::cerr << "Missing selector at position #" << argv_index << " '" << argv[argv_index] << "'." << std::endl;
          help();
          return 5;
          break;
      }
    }

    delete optimizer;
    if (argv_index < argc && strlen(argv[1]) < 2 || argv[1][0] != '-') {
      std::cerr << "Missing selector at position #" << argv_index << " '" << argv[argv_index] << "'." << std::endl;
      help();
      return 4;
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

  return 0;
}
