// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/similarity.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "(-x|-n) (-p|-c) <knowledge_base_index> (<query_index> <mine> <output>)+\tParse file with mined residues and identify the most similar:\n";
  std::cout << "                                                                       \t\t-p\tproteins; or\n";
  std::cout << "                                                                       \t\t-c\tchains. And use:\n";
  std::cout << "                                                                       \t\t-x\tthe higher of query/ template scores; or\n";
  std::cout << "                                                                       \t\t-c\tthe lower of query/ template scores.\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 7;
  const char* args[] = {argv[0],
    "-x",
    "-p",
    "C:\\Inspire\\pdb-new-all\\construction\\residues.ind",
    "C:\\Inspire\\pdb-new-all\\prediction\\residues.ind",
    "C:\\Inspire\\pdb-new-all\\prediction\\mined.med",
    "C:\\Inspire\\pdb-new-all\\prediction\\similarity"
  };
  argv = args;
#endif // TESTING

  if (argc < 7 || argc %3 != 1) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }
  if (strlen(argv[1]) < 2 || argv[1][0] != '-') {
    std::cerr << "The first argument must identify the score selector (ma<x>imal or mi<n>imal).";
    help();
    return 0;
  }
  if (strlen(argv[2]) < 2 || argv[2][0] != '-') {
    std::cerr << "The second argument must identify the type of comparator (<p>roteins or <c>hains).";
    help();
    return 0;
  }

  inspire::backend::Combinator* combinator = nullptr;
  inspire::backend::Similariter* similariter = nullptr;
  try {
    switch (argv[1][1]) {
      case 'x':
        combinator = new inspire::backend::MaximalCombinator();
        break;
      case 'n':
        combinator = new inspire::backend::MinimalCombinator();
        break;
      default:
        break;
    }
    switch (argv[2][1]) {
      case 'p':
        similariter = new inspire::backend::ProteinSimilariter(combinator, argv[3]);
        break;
      case 'c':
        similariter = new inspire::backend::ChainSimilariter(combinator, argv[3]);
        break;
      default:
        break;
    }
    for (size_t i = 4; i < argc; i+=3) {
      similariter->analyze(argv[i], argv[i+1], argv[i+2]);
    }
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
  if (similariter != nullptr) {
    delete similariter;
  }
  if (combinator != nullptr) {
    delete combinator;
  }

  return 0;
}
