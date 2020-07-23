// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/predict.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Selects a winning label for each query residue according the given metric.\n\n";

  std::cout << "Usage:\t((-f<THRESHOLDS> | -w<WEIGHTS>) ([-] <STATISTICS-FILE> <OUTPUT-PATH>)+)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t-f<THRESHOLDS>     \tFractional binary tree classificator with the given thresholds will be used for the prediction,\n";
  std::cout << "        \t                   \ti.e. it is tested whether number of hits with label 0 divided by a number of hits with labels 0 or 1 is at least threshold_0_1 and\n";
  std::cout << "        \t                   \tso on with the winning label and label 2...\n";
  std::cout << "        \t                   \t<THRESHOLDS> must contains number of 2-combination of x space-separated decimal numbers greater than or equal to 0 and less than or equal to 1,\n";
  std::cout << "        \t                   \twhere x is number of labels in <STATISTICS-FILE> and order of labels in STATISTICS-FILE's header must be preserved.\n";
  std::cout << "        \t-w<WEIGHTS>        \tWeighted classificator with the given weights will be used for the prediction.\n";
  std::cout << "        \t                   \tI.e. it is selected a label whose number of hits times weight of the label is the highest.\n";
  std::cout << "        \t                   \t<WEIGHTS> must contains x-1 space-separated decimal numbers greater than or equal to 0,\n";
  std::cout << "        \t                   \twhere x is a number of labels in <STATISTICS-FILE>.\n";
  std::cout << "        \t                   \tOrder weights must be the same as the order of labels in STATISTICS-FILE's;\n";
  std::cout << "        \t                   \tA weight  of the last label is calculated as 1-SUM(<WEIGHTS>), thus SUM(<WEIGHTS>) must be less than or equal to one.\n";
  std::cout << "        \t<STATISTICS-FILE>  \tSet of mined statistics used for prediction.\n";
  std::cout << "        \t                   \tPreceding '-' is mandatory, if a <STATISTICS-FILE> starts with a hyphen-minus sign.\n";
  std::cout << "        \t<OUTPUT-PATH>      \tWhere to store output file.\n";
  std::cout << "        \t                   \tIf <OUTPUT-PATH> is empty or ends with a directory separator, <STATISTICS-FILE>'s basename is used as the file name with '.pec' as an extension.\n";
  std::cout << "        \t                   \tIf <OUTPUT-PATH> does not end with '.pec' extension, the extension is appended.\n";
  std::cout << "        \t-h                 \tShow informations about the program\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-predict.log");
  argc = 4;
  const char* args[] = {argv[0],
    "-f 0.33333",
    "C:\\Inspire\\test\\mined.sas",
    "C:\\Inspire\\test\\mined"
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

  inspire::backend::Predictor* predictor = nullptr;
  try {
    size_t argv_index = 1;
    while (argv_index < argc && strlen(argv[argv_index]) >= 2 && argv[argv_index][0] == '-') {
      switch (argv[argv_index][1]) {
        case 'f':
          if (predictor != nullptr) {
            delete predictor;
          }
          predictor = new inspire::backend::FractionalPredictor(common::string::trim(std::string(argv[argv_index]).substr(2)));
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
            std::string output(argv[argv_index++]);
            predictor->predict(input, output);
          }
          break;
        case 'w':
          if (predictor != nullptr) {
            delete predictor;
          }
          predictor = new inspire::backend::WeightedPredictor(common::string::trim(std::string(argv[argv_index]).substr(2)));
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
            std::string output(argv[argv_index++]);
            predictor->predict(input, output);
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
  if (predictor != nullptr) {
    delete predictor;
  }

  return 0;
}
