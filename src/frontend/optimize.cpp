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

  std::cout << "Optimize parameters of a classifier based on the selected feature and mined statistics.\n\n";

  std::cout << "Usage:\t(-f<PRECISION> (-l<FEATURES-FILE> ([-] <STATISTICS-FILE> -o<OUTPUT-PATH>)+)+)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t-f<PRECISION>      \tFractional binary classificator's threshold will be optimised with the given precission,\n";
  std::cout << "        \t                   \ti.e. (i+0.5)/<PRECISION> for i from 0 to (<PRECISION>-1) is tested as threshold for ratio of positive and negative cases.\n";
  std::cout << "        \t                   \tThreshold is optimized using Matthews correlation coefficient as a prediction quality metric.\n";
  std::cout << "        \t-l<FEATURES-FILE>  \tReal labels of residues their classification is optimized.\n";
  std::cout << "        \t<STATISTICS-FILE>  \tSet of mined statistics used for optimization.\n";
  std::cout << "        \t<OUTPUT-PATH>      \tWhere to store output file.\n";
  std::cout << "        \t                   \tIf OUTPUT-PATH is empty or ends with a directory separator, STATISTICS-FILE's basename is used as the file name with '.sed' as an extension.\n";
  std::cout << "        \t                   \tIf OUTPUT-PATH does not end with '.sed' extension, the extension is appended.\n";
  std::cout << "        \t-h                 \tShow informations about the program\n\n";
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

  inspire::backend::Optimizer* optimizer = nullptr;
  try {
    size_t argv_index = 1;
    std::string function = "f1000";
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
          } else {
            std::cerr << "Unknown optimize '" << function << "'." << std::endl;
            help();
            return 7;
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
  if (optimizer != nullptr) {
    delete optimizer;
  }

  return 0;
}
