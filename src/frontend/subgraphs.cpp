// subgraphs.cpp : Defines the entry point for the console application.
//

//#define TESTING

#include "../backend/subgraphs.h"

void help() {
  std::cout << "Help\n\n";
  std::cout << "Extract subgraphs for residues in the knowledge-base.\n\n";
  std::cout << "<index> <positions> <path> (-c <limit>|-d <distance>|-e <distance> <limit>)+\n";
  std::cout << "\tExtract subgraphs for each residue in the index file that has defined position in the file <positions>.\n";
  std::cout << "\tGraphs are stored in a file '<path>{type}[<distance>_]<limit>', where {type} is 'c', 'd' or 'e'.\n\n";
  std::cout << "\tIn a subgraph is a central residue and:\n";
  std::cout << "\t-c:\t<limit>-nearest neighbours (if several neighbours are within the same distance, they all are taken);\n";
  std::cout << "\t-d:\tall residues that are at most <distance>Angstroems away (according to carbon alpha);\n";
  std::cout << "\t-e:\tall residues that are at most <limit> edges away and\n";
  std::cout << "\t\ttwo residues are considered to be connected by an edge if they are at most <distance> distant.\n";
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

