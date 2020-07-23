// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../backend/subgraphs.h"

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Extract subgraphs for each indexed residue (called as central residue) with defined position.\n\n";

  std::cout << "Usage:\t<INDEX-FILE> <POSITIONS-PATH> <OUTPUT-PATH> ( -c <LIMIT>+ | -d <DISTANCE>+ | -e <DISTANCE> <LIMIT>+ )+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INDEX-FILE>       \tPath to a index file\n";
  std::cout << "        \t<POSITIONS-PATH>   \tFeatures file with coordinates of residues\n";
  std::cout << "        \t<OUTPUT-PATH>      \tWhere to store output file(s).\n";
  std::cout << "        \t                   \tTo the <OUTPUT-PATH> is appended code of the used subgraph type(see that section) and '.sup' extension.\n";
  std::cout << "        \t-h                 \tShow informations about the program\n";
  std::cout << "    Subgraphs: in each subgraph is a central residue and\n";
  std::cout << "        \t-c <LIMIT>         \t<LIMIT>-nearest neighbours (if several neighbours are within the same distance as the <LIMIT>th-nearest neighbour, they all are taken).\n";
  std::cout << "        \t                   \tFilename suffix is 'c<LIMIT>.sup'.\n";
  std::cout << "        \t-d DISTANCE        \tAll residues that are at most <DISTANCE> Angstroems away.\n";
  std::cout << "        \t                   \tFilename suffix is 'd<DISTANCE>.sup'.\n";
  std::cout << "        \t-e DISTANCE LIMIT  \tAll residues that are at most <LIMIT> edges away, and\n";
  std::cout << "        \t                   \ttwo residues are considered to be connected by an edge if they are at most <DISTANCE> Angstroems distant.\n";
  std::cout << "        \t                   \tFilename suffix is 'e<DISTANCE>_<LIMIT>.sup' (<DISTANCE> and <LIMIT> are separated by underscore).\n\n";

  std::cout << "Note: Currently, DISTANCE is parsed and then converted back to string to create filenames, so formatting is not preserved.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 29;
  const char* arg[] = {argv[0], "C:\\Inspire\\test\\residue.ind", "C:\\Inspire\\test\\features.tur", "C:\\Inspire\\test\\s-",
    "-c", "0", "6", "12", "18",
    "-d", "0", "6", "12", "18",
    "-e", "5", "1", "2", "3",
    "-e", "6", "1", "2", "3",
    "-e", "7", "1", "2", "3"
  };
  argv = arg;
#endif // TESTING
  if (argc < 5) {
    help();
    return 0;
  }

  try {
    inspire::backend::ChainSubgraphs subgraphs(argv[1], argv[2], argv[3]);

    for (size_t i = 4; i < argc; i++) {
      if (strlen(argv[i]) == 2 && argv[i][0] == '-') {
        switch (argv[i][1]) {
          case 'c':
            while (++i < argc && (strlen(argv[i]) == 0 || argv[i][0] != '-')) {
              subgraphs.k_nearest(atoi(argv[i]));
            }
            --i;
            break;
          case 'd':
            while (++i < argc && (strlen(argv[i]) == 0 || argv[i][0] != '-')) {
              subgraphs.distance_limit(atof(argv[i]));
            }
            --i;
            break;
          case 'e':
            if (++i < argc && (strlen(argv[i]) == 0 || argv[i][0] != '-')) {
              double length = atof(argv[i]);
              while (++i < argc && (strlen(argv[i]) == 0 || argv[i][0] != '-')) {
                subgraphs.edge_limit(length, atoi(argv[i]));
              }
            }
            --i;
            break;
          default:
            std::cerr << "Unknown argument " << argv[i];
            break;
        }
      } else {
        std::cerr << "Unknown argument " << argv[i];
      }
    }

    if (subgraphs.empty()) {
      std::cerr << "No valid subgraph type selected.";
      help();
      return 11;
    }

    subgraphs.extract_subgraphs();
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
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

