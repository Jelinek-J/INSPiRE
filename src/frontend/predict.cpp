// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/predict.h"
#include "../elemental/exception.h"
#include "../elemental/string.h"
#include <iostream>
#include <fstream>

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "((-f <tresholds>) ([-] <mine> <output>)+)+\tSelect the winning label according the given metric.\n";
  std::cout << "                                      \t-f <thresholds>\tlinear selector will be used;";
  std::cout << "                                      \t               \t<thresholds> must contains number of 2-combination of <x> thresholds,";
  std::cout << "                                      \t               \twhere <x> is number of labels in <mine> and they are ordered from the earlier to later\n";
  std::cout << "                                      \t-\tis required if a <mine> starts with '-'.";
  std::cout << "NOTE: Principially, <interfaces> could be an arbitrary features file.\n\n";
}

int main(int argc, const char** argv) {
  if (argc < 4) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  try {
    size_t argv_index = 1;
    inspire::backend::Predictor* predictor = nullptr;
    while (argv_index < argc && strlen(argv[argv_index]) >= 2 && argv[argv_index][0] == '-') {
      switch (argv[argv_index][1]) {
        case 'f':
          if (predictor != nullptr) {
            delete predictor;
          }
          predictor = new inspire::backend::FractionalPredictor(elemental::string::trim(std::string(argv[argv_index]).substr(2)));
          while (++argv_index < argc && (strlen(argv[argv_index]) < 2 || argv[argv_index][0] != '-')) {
            if (argv[argv_index] == "-") {
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

    delete predictor;
    if (argv_index < argc && strlen(argv[1]) < 2 || argv[1][0] != '-') {
      std::cerr << "Missing selector at position #" << argv_index << " '" << argv[argv_index] << "'." << std::endl;
      help();
      return 4;
    }
  } catch (const elemental::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
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
