// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../backend/random.h"
#include "../common/exception.h"

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Randomly selects 100 protein chains (at most one chain from each protein) with mutual similarity less than a given threshold.\n";
  std::cout << "Selected chains will be printed in format <PROTEIN-ID>.<CHAIN-ID>.\n\n";

  std::cout << "Usage:\t<OUTPUT-FILE> [-n<COUNT>] [-l<LIMIT>] [-f|-b|-c|-bc|-w] <INDEX-FILE> <SIMILARITY-FILE>\n";
  std::cout << "      \t<OUTPUT-FILE> ([-n<COUNT>] [-l<LIMIT>] [-f|-b|-c|-bc|-w] [-] <HEADER> <INDEX-FILE> <SIMILARITY-FILE>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<OUTPUT-FILE>    \tWhere to store identifiers of selected chains.\n";
  std::cout << "        \t-n<COUNT>        \tHow many of dissimilar chains should be selected. The default value is 100.\n";
  std::cout << "        \t-l<LIMIT>        \tLimit to consider chains as similar. The default value is 0.1.\n";
  std::cout << "        \t- <HEADER>       \tHeader to use in the output file for a following pair of input files. The hyphen/minus sign is mandatory only if <HEADER> starts with a hyphen/minus sign.\n";
  std::cout << "        \t<INDEX-FILE>     \tPath to a index file containing chains that can be selected.\n";
  std::cout << "        \t<SIMILARITY-FILE>\tPath to a similarity file defining similarity of chains.\n";
  std::cout << "        \t-h, --help       \tShow informations about the program.\n";
  std::cout << "\tIterators: What iterator was used to construction of the index file:\n";
  std::cout << "\t\t-b \tAll biomolecules and models, but only the first crystallographic transformation were used;\n";
  std::cout << "\t\t-c \tAll crystallographic transformations, but only the first biomolecule and model were used;\n";
  std::cout << "\t\t-bc\tAll biomolecules, models and crystallographic transformations were used;\n";
  std::cout << "\t\t-f \tOnly the first biomolecule from the first model with the first crystallographic transformation was used (this is the default option);\n";
  std::cout << "\t\t-w \tBoth biomolecules and crystallographic transformation were ignored, all chains were used as they were.\n\n";
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
  double limit = 0.1;

  try {
    bool first = true;
    inspire::backend::Random random(common::filesystem::complete(argv[1], "selected", ".let"));
    inspire::backend::ProteinIterator* it = nullptr;
    for (size_t argi = 2; argi < argc; ++argi) {
      if (!(first && argi+2 == argc) && strlen(argv[argi]) > 1 && argv[argi][0] == '-') {
        switch (argv[argi][1]) {
          case 'n':
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
          case 'b':
            if (strlen(argv[argi]) > 2 && argv[argi][2] == 'c') {
              it = new inspire::backend::AllExceptAltLocIterator();
            } else {
              it = new inspire::backend::BiomoleculesIterator();
            }
            break;
          case 'c':
            it = new inspire::backend::FirstModelCrystallographicIterator();
            break;
          case 'f':
            it = new inspire::backend::FirstModelIterator();
            break;
          case 'w':
            it = new inspire::backend::ExplicitIterator();
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
        if (it == nullptr) {
          it = new inspire::backend::FirstModelIterator();
        }
        if (first && argi+2 == argc) {
          random.select(count, limit, it, argv[argi], argv[argi+1]);
        } else {
          random.select(count, limit, argv[argi], it, argv[argi+1], argv[argi+2]);
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

