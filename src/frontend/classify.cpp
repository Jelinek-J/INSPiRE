// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/classify.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Parse a file with mined residues, classify them according to given labels file and write summary statistics.\n\n";

  std::cout << "Usage:\t<FEATURES-FILE> [-f <FEATURE-NAMES>] (<MINED-FILE> <OUTPUT-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<MINED-FILE>        \tPath to a file with mined residues\n";
  std::cout << "        \t<FEATURES-FILE>     \tPath to a features file that will by used for classification\n";
  std::cout << "        \t-f <FEATURE-NAMES>  \tSpecify, what features should be used for clasification; multiple features can be separated by a directory separator\n";
  std::cout << "        \t<OUTPUT-PATH>       \tWhere to store output file.\n";
  std::cout << "        \t                    \tIf <OUTPUT-PATH> is a directory or ends with a directory separator, <PREDICTION-FILE>'s basename is used as the file name;\n";
  std::cout << "        \t                    \tif the path is not a directory but does not but not ends with '.sas' extension, the extension is appended.\n";
  std::cout << "        \t-h                  \tShow informations about the program\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 6;
  const char* args[] = {argv[0],
    "C:\\Inspire\\repair\\interface.tur",
    "-f", "interface",
    "C:\\Inspire\\repair\\mined.med",
    "C:\\Inspire\\repair\\stats"
  };
  argv = args;
#endif // TESTING

  if ((argc < 6 && (argc < 4 || common::string::starts_with(argv[2], "-f"))) || ((argc % 2) != 0)) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  inspire::backend::Classifier* classifier = nullptr;
  try {
    size_t start;
    if (common::string::starts_with(argv[2], "-f")) {
      std::vector<std::string> features;
      std::stringstream parts(argv[3]);
      std::string part;
      while (std::getline(parts,part,common::filesystem::directory_separator)) {
        features.push_back(part);
      }
      classifier = new inspire::backend::Classifier(argv[1], features);
      start = 4;
    } else {
      classifier = new inspire::backend::Classifier(argv[1]);
      start = 2;
    }
    for (size_t i = start; i < argc; i+=2) {
      classifier->classify(argv[i], argv[i+1]);
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
  if (classifier == nullptr) {
    delete classifier;
  }

  return 0;
}
