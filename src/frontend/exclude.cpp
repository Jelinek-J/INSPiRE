// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/exclude.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "<knowledge_base> <query> <output> (<protein_id>)+\tGenerate exclude files to define what fingerprints should be skipped in knowledge-base.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 62;
  const char* args[] = {argv[0],
    "C:\\Inspire\\basic2\\construction\\residues.ind",
    "C:\\Inspire\\basic2\\gvin\\residues.ind",
    "C:\\Inspire\\basic2\\gvin\\",
    "1A22", "1AOK", "1APY", "1AUT", "1AVO", "1B9X", "1BCC", "1BFO", "1BH8", "1BVN", "1BWV", "1DEV", "1E3A", "1F60", "1FLC", "1FPN",
    "1FS0", "1G2C", "1G72", "1G8J", "1G9M", "1GKA", "1H8T", "1HDM", "1HFE", "1HWM", "1I2M", "1I7Q", "1I8I", "1II8", "1JIW", "1JKJ",
    "1JSD", "1KDQ", "1KIU", "1KMI", "1MCO", "1MDA", "1MG2", "1NEX", "1NPE", "1OCC", "1OF5", "1ORY", "1P5U", "1Q90", "1QGK", "1R8O",
    "1SB2", "1SKY", "1SQB", "1T0F", "1T0Q", "1ULI", "1W1W", "1Y14", "1ZBA", "2PCD"
  };
  argv = args;
#endif // TESTING

  if (argc < 5) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  try {
    inspire::backend::Excluder excluder;
    for (size_t i = 4; i < argc; i++) {
      excluder.add(argv[i]);
    }
    excluder.generate(argv[1], argv[2], argv[3]);
  } catch (const common::exception::TitledException& e) {
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
