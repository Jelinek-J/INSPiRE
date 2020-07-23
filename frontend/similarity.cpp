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

  std::cout << "Parse file with mined residues and identify the most similar proteins/chains.\n\n";

  std::cout << "Usage:\t-x|-n -p|-c <KNOWLEDGE-BASE-INDEX> (<QUERY-INDEX> <MINED-FILE> <OUTPUT-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<KNOWLEDGE-BASE-INDEX>  \tPath to the index file of the knowledge-base used for mining\n";
  std::cout << "        \t<QUERY-INDEX>           \tPath to the index file of a query\n";
  std::cout << "        \t<MINED-FILE>            \tPath to a file with mined residues\n";
  std::cout << "        \t<OUTPUT-PATH>           \tWhere to store output file.\n";
  std::cout << "        \t                        \tIf <OUTPUT-PATH> is empty or ends with a directory separator, <MINED-FILE>'s basename is used as the file name with '.rty' as an extension.\n";
  std::cout << "        \t                        \tIf OUTPUT-PATH does not end with '.rty' extension, the extension is appended.\n";
  std::cout << "        \t-h                      \tShow informations about the program\n";
  std::cout << "    Entities: identify the most similar\n";
  std::cout << "        \t-p                      \tproteins;\n";
  std::cout << "        \t-c                      \tchains.\n";
  std::cout << "    Aggregator: (how should be combined scores of query protein/ chain and template protein/ chain)\n";
  std::cout << "        \t-x                      \tthe higher of query/ template scores will be taken;\n";
  std::cout << "        \t-n                      \tthe lower of query/ template scores will be taken.\n\n";
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
