// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../backend/validate.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Check what complexes are valid to be in a knowledge-base.\n\n";
  std::cout << "<index_file> <index_log> <aminoacid_file> <aminoacid_log> <coordinates_file> <coordinates_log> (<pdb_directory>|<pdb_file>)+\n";
  std::cout << "-h\n\n";
  std::cout << "\t<index_file>      \tPath to a index file;\n";
  std::cout << "\t<index_log>       \tPath to a log file created during <index_file> construction;\n";
  std::cout << "\t<aminoacid_file>  \tPath to a feature file with residue type three-letters codes;\n";
  std::cout << "\t<aminoacid_log>   \tPath to a file with transformation - it is used to define valid residue codes;\n";
  std::cout << "\t<coordinates_file>\tPath to a feature file with measured composition of residues;\n";
  std::cout << "\t<coordinates_log> \tPath to a file saying what are valid compositions of residue types;\n";
  std::cout << "\t<pdb_path>        \tPath to a protein or a directory with proteins that should be validated;\n";
  std::cout << "\t-h                \tShow this help.\n";
  std::cout << "\tIterators: What iterator was used to construction of the index file\n";
  std::cout << "\t           (if no iterator is specified, only the first biomolecule from the first model with the first crystallographic transformation was used):\n";
  std::cout << "\t\t-b \tAll biomolecules and models, but only the first crystallographic transformation were used;\n";
  std::cout << "\t\t-c \tAll crystallographic transformations, but only the first biomolecule and model were used;\n";
  std::cout << "\t\t-bc\tAll biomolecules, models and crystallographic transformations were used;\n";
  std::cout << "\t\t-w \tBoth biomolecules and crystallographic transformation were ignored, all chains were used as they were.\n\n";
  std::cout << "Program check, whether:\n";
  std::cout << "\t1) File contains REMARKS 350 (or equivalent in mmcif/ xml);\n";
  std::cout << "\t2) Is not a monomer;\n";
  std::cout << "\t3) Does not contain DNA, RNA, unknown amino acid (X, B, J, Z), modified aminoacid nor unknown residue;\n";
  std::cout << "\t4) Contains no HETATM except H2O;\n";
  std::cout << "\t5) No chain contains less than 20 aminoacids;\n";
  std::cout << "\t6) Contains less than 10 incomplete aminoacids;\n";
  std::cout << "\t7) No chain contains more than 1% of incomplete amino acids;\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\test\\error-validate.log");
  argc = 10;
  const char* arg[] = { argv[0],
    "C:\\Inspire\\validate\\residue.ind", "C:\\Inspire\\validate\\index.log",
    "C:\\Inspire\\validate\\aminoacid.tur", "C:\\Inspire\\validate\\pure.moc",
    "C:\\Inspire\\validate\\composition.tur", "C:\\Inspire\\validate\\composition.cit",
    "-b", "C:\\Inspire\\validate\\2010\\"
  };
  argv = arg;
#endif // TESTING
  if (argc < 8) {
    if (argc != 1 && argv[1] != "-h") {
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
    size_t files = 7;
    inspire::backend::ProteinIterator* it;
    if (std::strlen(argv[7]) > 1 && argv[7][0] == '-') {
      if (argv[7] == std::string("-c")) {
        it = new inspire::backend::FirstModelCrystallographicIterator();
      } else if (argv[7] == std::string("-b")) {
        it = new inspire::backend::BiomoleculesIterator();
      } else if (argv[7] == std::string("-bc")) {
        it = new inspire::backend::AllExceptAltLocIterator();
      } else if (argv[7] == std::string("-w")) {
        it = new inspire::backend::ExplicitIterator();
      } else {
        std::cerr << "ERROR: Modifier '" << argv[7] << "' is not currently supported." << std::endl;
        help();
        return 1;
      }
      ++files;
    } else {
      it = new inspire::backend::FirstModelIterator();
    }
    inspire::backend::Validate validate;
    for (; files < argc; ++files) {
      validate.add_files(argv[files]);
    }
    validate.validate(argv[1], it, argv[2], argv[3], argv[4], argv[5], argv[6]);
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    help();
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

