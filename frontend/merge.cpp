// mine.cpp : Defines the entry point for the console application.
//

#include "../common/filesystem.h"
#include "../common/string.h"
#include "../common/exception.h"
#include "../backend/merge.h"
#include <iostream>
#include <fstream>
#include <sstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Pairs fingerprints with the same index in different sets of fingerprints and concatenate them.\n\n";

  std::cout << "Usage:\t(<FINGERPRINTS> <FINGERPRINTS> <OUTPUT-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<FINGERPRINTS>  \tSet of fingerprints.\n";
  std::cout << "        \t                \tBoth sets must be of the same type (knowledge-base, or query) and they must have the same directory structure (in the case of knowledge-base file format).\n";
  std::cout << "        \t<OUTPUT-PATH>   \tWhere to store output fingerprints.\n";
  std::cout << "        \t                \tIf <FINGERPRINTS>s are in text format and <OUTPUT-PATH> is empty or ends with a directory separator,\n";
  std::cout << "        \t                \t\t<FINGERPRINTS>'s basenames are concatenated with hyphen-minus sign is used as a file name.\n";
  std::cout << "        \t                \t\tIf filename does not end with '.fit' extension, the extension is appended.\n";
  std::cout << "        \t                \tIf <FINGERPRINTS>s are in binary format and <OUTPUT-PATH> is considered as a directory. Directory structure of <FINGERPRINTS>s is preserved.\n";
  std::cout << "        \t-h                 \tShow informations about the program\n\n";

  std::cout << "Notes:\tEach index number must be unique in each set of fingerprints.\n";
  std::cout << "      \tFingerprints must be sorted by index numbers due to performancy reasons.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 7;
  const char* args[] = {argv[0],
    "C:\\Inspire\\basic2\\gvin2\\fingerprints.fit",
    "C:\\Inspire\\basic2\\gvin3\\fingerprints.fit",
    "C:\\Inspire\\basic2\\gvin4\\fingerprints.fit",
    "C:\\Inspire\\basic2\\fingerprints\\",
    "C:\\Inspire\\basic2\\fingerprints2\\",
    "C:\\Inspire\\basic2\\fingerprints3\\"
  };
  argv = args;
#endif // TESTING
  // NOTE: In the case of default (cs-CZ?) it writes 'á' instead of (non-breakable?) space.
//  std::cout.imbue(std::locale("en_UK"));
  // Special case, the later initialization is useless
  if (argc < 4 || argc % 3 != 1) {
    help();
    return 0;
  }


  try {
    for (size_t i = 1; i < argc; i+=3) {
      inspire::backend::Merger::merge(argv[i], argv[i+1], argv[i+2]);
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

