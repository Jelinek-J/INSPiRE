// mine.cpp : Defines the entry point for the console application.
//

#include "../elemental/filesystem.h"
#include "../elemental/string.h"
#include "../elemental/exception.h"
#include "../backend/merge.h"
#include <iostream>
#include <fstream>
#include <sstream>

#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "(<first> <second> <output>)+\tMerges fingerprints in <first> and <second> and save it into <output>\n\n";
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
  } catch (const elemental::exception::TitledException& e) {
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

