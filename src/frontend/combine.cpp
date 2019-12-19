// mine.cpp : Defines the entry point for the console application.
//

#include "../backend/combine.h"
#include "../common/exception.h"
#include "../common/string.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Combine multiple optimization files into one file\n\n";

  std::cout << "Usage:\t[-m] [-o] <OUTPUT-FILE> -t (<HEADER> <INPUT-FILE>)+\n";
  std::cout << "      \t[-m] [-o] <OUTPUT-FILE> [-p<PREFIX>] [-s<SUFFIX>] [-i] <INPUT-FILE>+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t-o <OUTPUT-FILE>   \tWhere to store the output file. '-o' is mandatory only in the case <OUTPUT-FILE> starts with minus sign.\n";
  std::cout << "        \t-i <INPUT-FILE>    \tInput file that should be added to the output file. '-i' is mandatory only in the case the first <INPUT-FILE> starts with minus sign and implicit headers.\n";
  std::cout << "        \t-m                 \tUse median to combine values from input files. (Defaultly, the mean value is used.)\n";
  std::cout << "        \t-t                 \tA flag indicating that each input file will be preceded by a header used in the output file.\n";
  std::cout << "        \t<HEADER>           \tHeader that shoul be used in the output file to represen the following input file.\n";
  std::cout << "        \t-p<PREFIX>         \tPrefix <PREFIX> should be removed from input file names to get corresponding headers for input files.\n";
  std::cout << "        \t-s<SUFFIX>         \tSuffix <SUFFIX> should be removed from input file names to get corresponding headers for input files.\n";
  std::cout << "        \t-h                 \tShow informations about the program\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-optimize.log");
  argc = 7;
  const char* args[] = {argv[0],
    "C:\\Inspire\\combine\\pure-sasapercentil8.den",
    "-t",
    "test",
    "C:\\Inspire\\combine\\2LTD.B-pure-sasapercentil8.den",
    "pokus",
    "C:\\Inspire\\combine\\3VNI.D-pure-sasapercentil8.den"
  };
  argv = args;
#endif // TESTING

  if (argc < 4) {
    if (argc > 1 && strcmp(argv[1], "-h")) {
      std::cerr << "Not enough arguments" << std::endl;
    }
    help();
    return 0;
  }

  inspire::backend::MutuallyOptimizer* optimizer;
  try {
    size_t argv_index = 1;
    bool median = false;
    if (strcmp("-m", argv[argv_index]) == 0) {
      median = true;
      ++argv_index;
    }
    if (strcmp("-o", argv[argv_index]) == 0) {
      ++argv_index;
    }
    if (median) {
      optimizer = new inspire::backend::MedianMutuallyOptimize(argv[argv_index]);
    } else {
      optimizer = new inspire::backend::AverageMutuallyOptimize(argv[argv_index]);
    }

    if (strcmp("-t", argv[++argv_index])) {
      std::string prefix;
      if (common::string::starts_with(argv[argv_index], "-p")) {
        prefix = std::string(argv[argv_index]).substr(2);
        ++argv_index;
      }
      std::string suffix;
      if (common::string::starts_with(argv[argv_index], "-s")) {
        suffix = std::string(argv[argv_index]).substr(2);
        ++argv_index;
      }
      if (strcmp("-i", argv[argv_index]) == 0) {
        ++argv_index;
      }
      for ( ; argv_index < argc; ++argv_index) {
        std::string header = argv[argv_index];
        if (common::string::starts_with(header, prefix)) {
          header = header.substr(prefix.size());
        } else {
          throw common::exception::TitledException("Filename '" + std::string(argv[argv_index]) + "' does not start with the given prefix '" + prefix + "'.");
        }
        if (common::string::ends_with(header, suffix)) {
          header = header.substr(0, header.size() - suffix.size());
        } else {
          throw common::exception::TitledException("Filename '" + std::string(argv[argv_index]) + "' does not end with the given suffix '" + suffix + "'.");
        }
        optimizer->add_input(header, argv[argv_index]);
      }
    } else {
      ++argv_index;
      if ((argc-argv_index)%2) {
        throw common::exception::TitledException("Invalid pairing of headers and input files: even number of corresponding arguments expected.");
      }
      while (argv_index < argc) {
        optimizer->add_input(argv[argv_index], argv[argv_index+1]);
        argv_index += 2;
      }
    }
    optimizer->combine();
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
  if (optimizer != nullptr) {
    delete optimizer;
  }

  return 0;
}
