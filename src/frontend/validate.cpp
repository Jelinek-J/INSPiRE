// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../backend/validate.h"
#include <iostream>
#include <fstream>

//#define TESTING

void help() {
  std::cout << "Help\n\n";

  std::cout << "Check what complexes are valid to be in a knowledge-base\n";
  std::cout << "<index_file> <index_log> <aminoacid_file> <aminoacid_log> <coordinates_file> <coordinates_log> (<pdb_directory>|<pdb_file>)+\t.\n";
  std::cout << "\tCheck, whether:\n";
  std::cout << "\t\t1) File contains REMARKS 350 (or equivalent in mmcif/ xml);\n";
  std::cout << "\t\t2) Is not a monomer;\n";
  std::cout << "\t\t3) Does not contain DNA, RNA, unknown amino acid (X, B, J, Z), modified aminoacid nor unknown residue;\n";
  std::cout << "\t\t4) Contains no HETATM except H2O;\n";
  std::cout << "\t\t5) No chain contains less than 20 aminoacids;\n";
  std::cout << "\t\t6) Contains less than 10 incomplete aminoacids;\n";
  std::cout << "\t\t7) No chain contains more than 1% of incomplete amino acids;\n";
  std::cout << "Each identifier (chain, complex) also applies to the following lines until it is redefined.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\test\\error-validate.log");
  argc = 10;
  const char* arg[] = { argv[0], "C:\\Inspire\\validate\\selected\\",
    "C:\\Inspire\\validate\\residue.ind", "C:\\Inspire\\validate\\index.log",
    "C:\\Inspire\\validate\\aminoacid.tur", "C:\\Inspire\\validate\\pure.moc",
    "C:\\Inspire\\validate\\composition.tur", "C:\\Inspire\\validate\\composition.cit",
    "-b", "C:\\Inspire\\validate\\2010\\"
  };
  argv = arg;
#endif // TESTING
  if (argc < 9) {
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
    size_t files = 8;
    inspire::backend::ProteinIterator* it;
    if (std::strlen(argv[8]) > 1 && argv[8][0] == '-') {
      if (argv[8] == std::string("-c")) {
        it = new inspire::backend::FirstModelCrystallographicIterator();
      } else if (argv[8] == std::string("-b")) {
        it = new inspire::backend::BiomoleculesIterator();
      } else if (argv[8] == std::string("-bc")) {
        it = new inspire::backend::AllExceptAltLocIterator();
      } else if (argv[8] == std::string("-w")) {
        it = new inspire::backend::ExplicitIterator();
      } else {
        std::cerr << "ERROR: Modifier '" << argv[8] << "' is not currently supported." << std::endl;
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
    validate.validate(argv[1], argv[2], it, argv[3], argv[4], argv[5], argv[6], argv[7]);
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

