// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../backend/fingerprints.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Constructs fingerprints defined in settings for subgraphs specified in '.sup' files with features of nodes defined in '.tur' files.\n\n";

  std::cout << "Usage:\t(k|q) <OUTPUT-PATH> <SETTINGS-FILE> <INDEX-FILE> <NODES-FILE> <EDGES-FILE> (<FEATURES-PATH>)*\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INDEX-FILE>     \tPath to a index file\n";
  std::cout << "        \t<SETTINGS-FILE>  \tPath to a file defining how should be fingerprints constructed\n";
  std::cout << "        \t<FEATURES-PATH>  \tFeatures file or directory with features files required to construction of fingerprints or that could be used latter durring prediction\n";
  std::cout << "        \t<OUTPUT-PATH>    \tWhere to store output file(s).\n";
  std::cout << "        \t                 \tWith 'q' switcher:\tIf <OUTPUT-PATH> is empty or ends with a directory separator, 'query.fit' is used as the file name.\n";
  std::cout << "        \t                 \t                  \tIf <OUTPUT-PATH> does not end with '.fit' extension, the extension is appended.\n";
  std::cout << "        \t                 \tWith 'k' switcher:\t<OUTPUT-PATH> is interpreted as a directory and fingerprints are separated in directories and\n";
  std::cout << "        \t                 \t                  \tfiles named after features of central residues with '.fin' as a file extension.\n";
  std::cout << "        \t                 \t                  \tOrder of directories hieararchy is the same as the order of features in FEATURES-PATHs.\n";
  std::cout << "        \t-h               \tShow informations about the program\n";
  std::cout << "    Subgraphs Definition:\n";
  std::cout << "        \t<NODES-FILE>     \tDefines for what subgraphs should be fingerprints constructed\n";
  std::cout << "        \t<EDGES-FILE>     \tDefines what nodes are connected by an edge\n";
  std::cout << "    File Formats:\n";
  std::cout << "        \tk                \tFingerprints are stored in binary files splitted in directories based on features of central residues\n";
  std::cout << "        \tq                \tAll fingerprints are stored in a single text file\n\n";

  std::cout << "Notes:\tDue to performancy reasons (saving RAM consumption), there are following expectations:\n";
  std::cout << "      \t\tRecords in <INDEX-FILE> are grouped by 'protein_name', then by 'model_name' and finally by 'chain_name' columns;\n";
  std::cout << "      \t\tLines in <NODES-FILE>, <EDGES-FILE> and files in <FEATURES-PATH>s are sorted by index; and\n";
  std::cout << "      \t\tAll indices within a single row in both <NODES-FILE> and <EDGES-FILE> are within the same model, i.e. you cannot combine two pdb files nor two different models.\n";
  std::cout << "      \t<NODES-FILE>, <EDGES-FILE> have the same file format, but their interpretation is different:\n";
  std::cout << "      \t\tin the first case the whole row identify nodes in the graph with first node as a central residue;\n";
  std::cout << "      \t\tin the second case the first node is connected with all other nodes by an edge.\n";
  std::cout << "      \tAlthough <EDGES-FILE>'s file format allows interpretation as directed edges, they are used as undirected edges.\n\n";  
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\test\\error-fingerprints.log");
  argc = 10;
  const char* arg[] = { argv[0], "k", "C:\\Inspire\\test\\fingerprints\\expand",
    "C:\\Inspire\\test\\fingerprints\\mass\\settings.json",
    "C:\\Inspire\\test\\construction\\residues.ind",
  	/**/ "C:\\Inspire\\test\\construction\\nodes.sup" /*/ "D:\\fingerprints\\old\\benchmark\\gvin-c_12.sup" /**/,
    "C:\\Inspire\\test\\construction\\edges.sup",
  	//"C:\\Inspire\\repair\\affinity.tur",
    //"C:\\Inspire\\repair\\interface.tur",
    //"C:\\Inspire\\repair\\aminoacid.tur",
    //"C:\\Inspire\\repair\\rasa10e.tur",
  	//"C:\\Inspire\\repair\\temperature10.tur",
    "C:\\Inspire\\test\\construction\\integer.tur",
    "C:\\Inspire\\test\\construction\\double.tur",
    "C:\\Inspire\\test\\construction\\word.tur"
  };
  argv = arg;
#endif // TESTING
  if (argc < 7) {
    if (argc > 1) {
      std::cerr << "Unexpected number of arguments:";
      for (size_t i = 1; i < argc; i++) {
        std::cerr << "\n[" << i << "]:\t" << argv[i];
      }
      std::cerr << std::endl;
    }
    help();
    return 0;
  }

  try {
    inspire::backend::FingerprintWriter fingerprints(argv[4], argv[5]);
    for (size_t i = 7; i < argc; i++) {
      fingerprints.add_features(argv[i]);
    }
    if (strlen(argv[1]) != 1) {
      std::cerr << "Unexpected length of type switche: '" << argv[1] << "'" << std::endl;
      help();
    }
    switch (argv[1][0]) {
      case 'k':
        fingerprints.process(argv[3], argv[6], argv[2], inspire::backend::FingerprintFormat::Binary);
        break;
      case 'q':
        fingerprints.process(argv[3], argv[6], argv[2], inspire::backend::FingerprintFormat::Text);
        break;
      default:
        break;
    }
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

