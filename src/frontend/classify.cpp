// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/classify.h"
#include "../elemental/exception.h"
#include "../elemental/string.h"
#include <iostream>
#include <fstream>

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "<interfaces> [-f <feature>] (<mine> <output>)+\tParse file with mined residues and clasify them according to labels in <interfaces>\n";
  std::cout << "                                              \t-f\tuse only <feature> for labeling; multiple features can be separated by a directory separator\n";
  std::cout << "NOTE: Principially, <interfaces> could be an arbitrary features file.\n\n";
}

int main(int argc, const char** argv) {
  if ((argc < 6 && (argc < 4 || elemental::string::starts_with(argv[2], "-f"))) || ((argc % 2) != 0)) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  try {
    inspire::backend::Classifier* classifier;
    size_t start;
    if (elemental::string::starts_with(argv[2], "-f")) {
      std::vector<std::string> features;
      std::stringstream parts(argv[3]);
      std::string part;
      while (std::getline(parts,part,elemental::filesystem::directory_separator)) {
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
    delete classifier;
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
