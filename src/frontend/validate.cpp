// subgraphs.cpp : Defines the entry point for the console application.
//

#include "../common/exception.h"
#include "../backend/validate.h"
#include <iostream>
#include <fstream>

//#define TESTING

static double DEFAULT_INCOMPLETE_REL = 0.05;
static size_t DEFAULT_INCOMPLETE_ABS = 50;
static size_t DEFAULT_MIN_LENGTH = 20;
static size_t DEFAULT_MIN_CHAINS = 2;

void help() {
  std::cout << "Help\n\n";

  std::cout << "Check what complexes are valid to be in a knowledge-base.\n\n";
  std::cout << "[-a[<incomplete_abs>]] [-b] [-c[<min_chains>]] [-d] [-l[<min_length>]] [-n] [-r[<incomplete_rel>]] [-u] <index_file> <index_log> <residues_file> <residues_valid> <compositions_file> <compositions_valid> [-ib|-ic|-ibc|-iw] <pdb_path>+\n";
  std::cout << "-h\n\n";
  std::cout << "\t\t<index_file>        \tPath to a index file;\n";
  std::cout << "\t\t<index_log>         \tPath to a log file created during <index_file> construction;\n";
  std::cout << "\t\t<residues_file>     \tPath to a feature file with residue type three-letters codes;\n";
  std::cout << "\t\t<residues_valid>    \tPath to a file with transformation - it is used to define valid residue codes;\n";
  std::cout << "\t\t<compositions_file> \tPath to a feature file with measured composition of residues;\n";
  std::cout << "\t\t<compositions_valid>\tPath to a file saying what are valid compositions of residue types;\n";
  std::cout << "\t\t<pdb_path>          \tPath to a protein or a directory with proteins that should be validated;\n";
  std::cout << "\t\t-h                  \tShow this help.\n";
  std::cout << "\tValidators: What should be checked and what thresholds should be used.\n";
  std::cout << "\t            Not all of the previously mentioned files are mandatory for all validations, if a file is not used in any test,\n";
  std::cout << "\t            a dummy value can be used for the corresponding argument (but the argument cannot be leaved out to have input less prone to errors).\n";
  std::cout << "\t\t-a<incomplete_abs>\tThere can be at most <incomplete_abs> incomplete (or missing) residues within the whole protein.\n";
  std::cout << "\t\t                  \tThe default value of <incomplete_abs> (if missing or empty) is 50.\n";
  std::cout << "\t\t                  \tMandatory are <index_file>, <index_log>, <residues_file>, <compositions_file>, <compositions_valid> and an used iterator.\n";
  std::cout << "\t\t-b                \tFilter out proteins without Remark 350 in PDB file (or equivalent in PDBML/PDBc file).\n";
  std::cout << "\t\t                  \tMandatory is <index_log>.\n";
  std::cout << "\t\t-c<min_chains>    \tThere must be at least <min_chains> chains within all biomolecules.\n";
  std::cout << "\t\t                  \tThe default value of <min_chains> (if missing or empty) is 2.\n";
  std::cout << "\t\t                  \tMandatory is <index_file>.\n";
  std::cout << "\t\t-d                \tAll biomolecules must have the same number of chains.\n";
  std::cout << "\t\t                  \tMandatory is <index_file>.\n";
  std::cout << "\t\t-l<min_length>    \tEach chain must have at least <min_length> residues.\n";
  std::cout << "\t\t                  \tThe default value of <min_length> (if missing or empty) is 20.\n";
  std::cout << "\t\t                  \tMandatory are <index_file>, <index_log> and an used iterator.\n";
  std::cout << "\t\t-n                \tThere cannot be any HETATM record.\n";
  std::cout << "\t\t                  \tMandatory are <pdb_path>s.\n";
  std::cout << "\t\t-r<incomplete_rel>\tThere can be at most 100*<incomplete_rel> % of residues incomplete (or missing) within each chain.\n";
  std::cout << "\t\t                  \tThe default value of <incomplete_rel> (if missing or empty) is 0.05.\n";
  std::cout << "\t\t                  \tMandatory are <index_file>, <index_log>, <residues_file>, <compositions_file>, <compositions_valid> and an used iterator.\n";
  std::cout << "\t\t-u                \tThere cannot be any unknown residue.\n";
  std::cout << "\t\t                  \tMandatory are <index_file>, <residues_file> and <residues_valid>.\n";
  std::cout << "\tIterators: What iterator was used to construction of the index file\n";
  std::cout << "\t           (if no iterator is specified, only the first biomolecule from the first model with the first crystallographic transformation was used):\n";
  std::cout << "\t\t-ib \tAll biomolecules and models, but only the first crystallographic transformation were used;\n";
  std::cout << "\t\t-ic \tAll crystallographic transformations, but only the first biomolecule and model were used;\n";
  std::cout << "\t\t-ibc\tAll biomolecules, models and crystallographic transformations were used;\n";
  std::cout << "\t\t-iw \tBoth biomolecules and crystallographic transformation were ignored, all chains were used as they were.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log("C:\\Inspire\\test\\error-validate.log");
  argc = 10;
  const char* arg[] = { argv[0],
    "C:\\Inspire\\validate\\residue.ind", "C:\\Inspire\\validate\\index.log",
    "C:\\Inspire\\validate\\aminoacid.tur", "C:\\Inspire\\validate\\pure.nor",
    "C:\\Inspire\\validate\\composition.tur", "C:\\Inspire\\validate\\composition.cit",
    "-b", "C:\\Inspire\\validate\\2010\\"
  };
  argv = arg;
#endif // TESTING

  if (argc < 7) {
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
    double incomplete_rel = -1.0;
    size_t incomplete_abs = -1;
    bool heteroatoms = false;
    bool unknown_residues = false;
    size_t min_length = 0;
    bool ambiguity = false;
    size_t min_chains = 0;
    bool remark_350 = false;

    size_t argi = 0;
    while ((++argi)+6 < argc && strlen(argv[argi]) > 1 && argv[argi][0] == '-') {
      switch (argv[argi][1]) {
        case 'a':
          if (strlen(argv[argi]) > 2) {
            incomplete_abs = std::stoll(std::string(argv[argi]).substr(2));
          } else {
            if (strlen(argv[argi+1]) == 0) {
              incomplete_abs = DEFAULT_INCOMPLETE_ABS;
              ++argi;
            } else if (argv[argi+1][0] == '-') {
              incomplete_abs = DEFAULT_INCOMPLETE_ABS;
            } else {
              incomplete_abs = std::atoll(argv[++argi]);
            }
          }
          break;
        case 'b':
          remark_350 = true;
          break;
        case 'c':
          if (strlen(argv[argi]) > 2) {
            min_chains = std::stoll(std::string(argv[argi]).substr(2));
          } else {
            if (strlen(argv[argi+1]) == 0) {
              min_chains = DEFAULT_MIN_CHAINS;
              ++argi;
            } else if (argv[argi+1][0] == '-') {
              min_chains = DEFAULT_MIN_CHAINS;
            } else {
              min_chains = std::atoll(argv[++argi]);
            }
          }
          break;
        case 'd':
          ambiguity = true;
          break;
        case 'n':
          heteroatoms = true;
          break;
        case 'l':
          if (strlen(argv[argi]) > 2) {
            min_length = std::stoll(std::string(argv[argi]).substr(2));
          } else {
            if (strlen(argv[argi+1]) == 0) {
              min_length = DEFAULT_MIN_LENGTH;
              ++argi;
            } else if (argv[argi+1][0] == '-') {
              min_length = DEFAULT_MIN_LENGTH;
            } else {
              min_length = std::atoll(argv[++argi]);
            }
          }
          break;
        case 'r':
          if (strlen(argv[argi]) > 2) {
            incomplete_rel = std::stod(std::string(argv[argi]).substr(2));
          } else {
            if (strlen(argv[argi+1]) == 0) {
              incomplete_rel = DEFAULT_INCOMPLETE_REL;
              ++argi;
            } else if (argv[argi+1][0] == '-') {
              incomplete_rel = DEFAULT_INCOMPLETE_REL;
            } else {
              incomplete_rel = std::stod(std::string(argv[++argi]));
            }
          }
          break;
        case 'u':
          unknown_residues = true;
          break;
        default:
          break;
      }
    }

    inspire::backend::Validate validator(argv[argi]);
    inspire::backend::ProteinIterator* it;
    std::set<std::string> paths;
    size_t argj = argi+6;
    if (std::strlen(argv[argj]) > 2 && argv[argj][0] == '-') {
      if (argv[argj] == std::string("-ic")) {
        it = new inspire::backend::FirstModelCrystallographicIterator();
      } else if (argv[argj] == std::string("-ib")) {
        it = new inspire::backend::BiomoleculesIterator();
      } else if (argv[argj] == std::string("-ibc")) {
        it = new inspire::backend::AllExceptAltLocIterator();
      } else if (argv[argj] == std::string("-iw")) {
        it = new inspire::backend::ExplicitIterator();
      } else {
        std::cerr << "ERROR: Modifier '" << argv[argj] << "' is not currently supported." << std::endl;
        help();
        return 1;
      }
    } else {
      --argj;
      it = new inspire::backend::FirstModelIterator();
    }
    while (++argj < argc) {
      paths.emplace(argv[argj]);
    }

    if (incomplete_rel != -1.0) {
      validator.validate_incomplete_rel(incomplete_rel, argv[argi], it, argv[argi+1], argv[argi+2], argv[argi+4], argv[argi+5]);
    }
    if (incomplete_abs != -1) {
      validator.validate_incomplete_abs(incomplete_abs, argv[argi], it, argv[argi+1], argv[argi+2], argv[argi+4], argv[argi+5]);
    }
    if (heteroatoms) {
      validator.validate_heteroatoms(paths);
    }
    if (unknown_residues) {
      validator.validate_residues(argv[argi], argv[argi+2], argv[argi+3]);
    }
    if (min_length != 0) {
      validator.validate_length(min_length, argv[argi], it, argv[argi+1]);
    }
    if (ambiguity) {
      validator.validate_ambiguity(argv[argi]);
    }
    if (min_chains != 0) {
      validator.validate_chains_count(min_chains, argv[argi]);
    }
    if (remark_350) {
      validator.validate_remark_350(argv[argi+1]);
    }
    validator.print_proteins();

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

