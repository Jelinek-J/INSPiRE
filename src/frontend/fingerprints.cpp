// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../backend/fingerprints.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";
  std::cout << "Construct fingerprints for specified residues.\n\n";
  std::cout << "<type> <output> <settings> <index> <nodes> <edges> (<tur_directory>|<tur_file>)*\n";
  std::cout << "\tConstruct fingerprints according to <settings> from every graph in <nodes> with edges from <edges>.\n";
  std::cout << "\tFeatures of individual residues will be taken from <tur_file> and/or all *.tur files in <tur_dicrectory>.\n";
  std::cout << "\tThe order of features in <tur_file> and/or <tur_dicrectory> is used for order of directories hierarchy.\n";
  std::cout << "\tFingerprints will be stored in <output>, it should be an empty directory or does not exist.\n";
  std::cout << "\t<type> can be 'k' for knowledge-base, or 'q' for query.\n";
  std::cout << "NOTE: Due to performancy reasons (for large knowledge bases it is not possible to load everything into RAM), there are following expectations:\n";
  std::cout << "\tRecords in <index> are grouped by protein_name and then by model_name columns;\n";
  std::cout << "\tLines in files <nodes>, <edges>, <tur_directory>s and <tur_file>s are sorted by index; and\n";
  std::cout << "\tAll indices within a single row in both <nodes> and <edges> are within the same model, i.e. you cannot combine two pdb files nor two different models.";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\error-subgraphs.log");
  argc = 9;
  const char* arg[] = { argv[0], "q", "C:\\Inspire\\repair\\fingerprints",
    "C:\\Inspire\\repair\\settings-511-t.json",
    "C:\\Inspire\\repair\\residue.ind",
  	/**/ "C:\\Inspire\\repair\\c6.sup" /*/ "D:\\fingerprints\\old\\benchmark\\gvin-c_12.sup" /**/,
  	"C:\\Inspire\\repair\\d6.sup",
  	//"C:\\Inspire\\repair\\affinity.tur",
    //"C:\\Inspire\\repair\\interface.tur",
    "C:\\Inspire\\repair\\aminoacid.tur",
    //"C:\\Inspire\\repair\\rasa10e.tur",
  	"C:\\Inspire\\repair\\temperature10.tur"
  };
  argv = arg;
#endif // TESTING
  if (argc < 7) {
    if (argc != 0) {
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

