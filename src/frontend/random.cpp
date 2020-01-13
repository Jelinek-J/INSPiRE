// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../backend/random.h"
#include "../common/exception.h"

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Randomly selects 100 protein chains (at most one chain from each protein) with mutual similarity less than a given threshold.\n";
  std::cout << "Selected chains will be printed in format <PROTEIN-ID>.<CHAIN-ID>.\n\n";

  std::cout << "Usage:\t<OUTPUT-FILE> [-c<COUNT>] [-l<LIMIT>] <INDEX-FILE> <SIMILARITY-FILE>\n";
  std::cout << "      \t<OUTPUT-FILE> ([-c<COUNT>] [-l<LIMIT>] [-] <HEADER> <INDEX-FILE> <SIMILARITY-FILE>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<OUTPUT-FILE>    \tWhere to store identifiers of selected chains.\n";
  std::cout << "        \t-c<COUNT>        \tHow many of dissimilar chains should be selected. The default value is 100.\n";
  std::cout << "        \t-l<LIMIT>        \tLimit to consider chains as similar. The default value is 0.1.\n";
  std::cout << "        \t- <HEADER>       \tHeader to use in the output file for a following pair of input files. The hyphen/minus sign is mandatory only if <HEADER> starts with a hyphen/minus sign.\n";
  std::cout << "        \t<INDEX-FILE>     \tPath to a index file containing chains that can be selected.\n";
  std::cout << "        \t<SIMILARITY-FILE>\tPath to a similarity file defining similarity of chains.\n";
  std::cout << "        \t-h, --help       \tShow informations about the program.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  argc = 3;
  const char* arg[] = { argv[0], "C:\\Inspire\\random\\residues.ind", "C:\\Inspire\\random\\mined.rty"
  };
  argv = arg;
#endif // TESTING
  if (argc < 4) {
    if (argc != 2 || argv[1] != "-h") {
      std::cerr << "Unexpected number of arguments:";
      for (size_t i = 1; i < argc; i++) {
        std::cerr << "\n[" << i << "]:\t" << argv[i];
      }
      std::cerr << std::endl;
    }
    help();
    return 0;
  }

  size_t count = 100;
  size_t limit = 0.1;

  try {
    bool first = true;
    inspire::backend::Random random(argv[1]);
    for (size_t argi = 2; argi < argc; ++argi) {
      if (!(first && argi+2 == argc) && strlen(argv[argi]) > 1 && argv[argi][0] == '-') {
        switch (argv[argi][1]) {
          case 'c':
            if (strlen(argv[argi]) == 2) {
              count = std::atoi(argv[++argi]);
            } else {
              count = std::stoi(std::string(argv[argi]).substr(2));
            }
            break;
          case 'l':
            if (strlen(argv[argi]) == 2) {
              limit = std::stod(std::string(argv[++argi]));
            } else {
              limit = std::stod(std::string(argv[argi]).substr(2));
            }
            break;
          default:
            std::cerr << "Unknown argument no. " << argi << ": '" << argv[argi] << "'\n" << std::endl;
            help();
            return 5;
        }
      } else {
        if (argv[argi] == "-") {
          ++argi;
        }
        if (first && argi+2 == argc) {
          random.select(count, limit, argv[argi], argv[argi+1]);
        } else {
          random.select(count, limit, argv[argi], argv[argi+1], argv[argi+2]);
          ++argi;
          first = false;
        }
        ++argi;
      }
    }
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
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

