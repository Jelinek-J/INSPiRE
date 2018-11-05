// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/assign.h"
#include "../elemental/exception.h"
#include "../elemental/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "-h\tPrint this message.\n\n";
  std::cout << "((-x|-c<delimiter>|-l) <index> <prediction> <output>)+\tParse prediction file, label them according to protein, model, chain and residue, and store as:\n";
  std::cout << "                                                      \t\t-x\tXML file;\n";
  std::cout << "                                                      \t\t-c\t<delimiter>-separated value file; or\n";
  std::cout << "                                                      \t\t-l\tstructured list.\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 5;
  const char* args[] = {argv[0],
    "-x",
    "C:\\Inspire\\repair\\residues.ind",
    "C:\\Inspire\\repair\\gvin-511-aa-c6.pec",
    "C:\\Inspire\\repair\\"
  };
  argv = args;
#endif // TESTING

  if (argc < 5 || argc %4 != 1) {
    if (argc > 1) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }
  try {
    for (size_t i = 1; i < argc; i+=4) {
      if (strlen(argv[i]) < 2 || argv[i][0] != '-') {
        std::cerr << "The first argument must identify the score selector (ma<x>imal or mi<n>imal).";
        help();
        return 4;
      }
      std::string output(argv[i+3]);
      switch (argv[i][1]) {
        case 'x':
          inspire::backend::Assignator::Xml(argv[i+1], argv[i+2], output);
          break;
        case 'c':
          if (strlen(argv[i]) < 3) {
            std::cerr << "Missing delimiter";
            help();
            return 5;
          }
          {
            char ch = argv[i][2];
            if (strlen(argv[i]) > 3 && argv[i][2] == '\\') {
              switch (argv[i][3]) {
                case 't':
                  ch = '\t';
                  break;
                case 'b':
                  ch = '\b';
                  break;
                case 'n':
                  ch = '\n';
                  break;
                case 'r':
                  ch = '\r';
                  break;
                case 'f':
                  ch = '\f';
                  break;
                case 'a':
                  ch = '\a';
                  break;
                case '\\':
                  ch = '\\';
                  break;
                default:
                  throw elemental::exception::TitledException("Unknown escape character");
                  help();
                  return 6;
              }
            }
            inspire::backend::Assignator::Csv(ch, argv[i+1], argv[i+2], output);
          }
          break;
        case 'l':
          inspire::backend::Assignator::List(argv[i+1], argv[i+2], output);
          break;
        default:
          break;
      }
    }
  } catch (const elemental::exception::TitledException& e) {
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
