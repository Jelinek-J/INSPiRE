// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/assign.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Parse prediction file, label them according to protein, model, chain and residue, and store them in a file.\n\n";

  std::cout << "Usage:\t(-x|-c<DELIMITER>|-l <INDEX-FILE> <PREDICTION-FILE> <OUTPUT-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<PREDICTION-FILE> \tPath to a file with prediction that should be transformed to a human-readable format\n";
  std::cout << "        \t<INDEX-FILE>      \tPath to a index file with labels of residues\n";
  std::cout << "        \t<OUTPUT-PATH>     \tWhere to store output file.\n";
  std::cout << "        \t                  \tIf <OUTPUT-PATH> is a directory or ends with a directory separator, <PREDICTION-FILE>'s basename is used as the file name;\n";
  std::cout << "        \t                  \tif the path is not a directory but does not but not ends with an extension corresponding to a choosen file format, the extension is appended.\n";
  std::cout << "        \t-h                \tShow informations about the program\n";
  std::cout << "    Output file format switchers:\n";
  std::cout << "        \t-x            \tXML file format. This file format has '.xml' as the extension.\n";
  std::cout << "        \t-c<DELIMITER> \tDelimiter-separated-value file format with <DELIMITER> as the delimiter.\n";
  std::cout << "        \t              \tThis file format has '.csv' as the extension.\n";
  std::cout << "        \t-l            \tSpace efficient variation of tab-separated-value file, where repeated values are ommited.\n";
  std::cout << "        \t              \tThis file format has '.pes' as the extension.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 5;
  const char* args[] = {argv[0],
    "-l",
    "C:\\Inspire\\assign\\residues.ind",
    "C:\\Inspire\\assign\\prediction.pec",
    "C:\\Inspire\\assign\\"
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
        std::cerr << "The first argument must identify the file format (<x>ml, <c>sv, or tab-separated va<l>ues with ommited repeated values).";
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
            std::cerr << "Missing delimiter" << std::endl;
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
                  throw common::exception::TitledException("Unknown escape character");
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
