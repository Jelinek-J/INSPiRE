// inspire.cpp : Defines the entry point for the console application.
//

#include "../backend/pdb.h"
#include "../backend/protein.h"
#include "../backend/index.h"
#include "../elemental/filesystem.h"
#include "../elemental/string.h"

#define TESTING

namespace inspire {
  namespace frontend {
    // Class for parsing pdb files and index them
    class Indexer {
#ifdef TESTING
      // Errors log for testing reasons
      std::ofstream LOG;
#endif // TESTING
      // Object that ensures indexing proteins
      backend::Indexer INDEX;
      // Filters pdb files to skip uninportant informations
      // NOTE: currently, no information is skipped
      backend::BasicFilter FILTER;

      public:
      // Initialize indexer with output index location in 'file' and protein iterator 'it'
      Indexer(std::string file, backend::ProteinIterator* it)
#ifdef TESTING
        : LOG("C:\\Inspire\\error.log"),
#endif // TESTING
        INDEX(file, it)
      {

      }

      ~Indexer() {
#ifdef TESTING
        LOG.flush();
        LOG.close();
#endif // TESTING
      }

      // Load a protein from pdb 'file' and index it
      void process(std::string file) {
        if (elemental::string::ends_with(file, ".pdb")) {
          try {
            std::ifstream input;
            input.open(file);
            std::cout << file << "    ";
            backend::Protein protein = backend::Pdb::parse_pdb(input, FILTER);
            std::cout << "parsed    ";
            INDEX.index(&protein);
            std::cout << "indexed\r";
          } catch (const elemental::exception::TitledException& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
#ifdef TESTING
            LOG << file << "    ERROR: " << e.what() << std::endl;
#endif // TESTING
          } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
#ifdef TESTING
            LOG << file << "    ERROR: " << e.what() << std::endl;
#endif // TESTING
          } catch (...) {
            std::cerr << "UNKNOWN ERROR" << std::endl;
#ifdef TESTING
            LOG << file << "    UNKNOWN ERROR" << std::endl;
#endif // TESTING
          }
        }
      }

    };
  }
}

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Create index of residues in the knowledge-base.\n";
  std::cout << "<index_name> [-(b|c|bc)] (<pdb_directory>|<pdb_file>)+\tIndex all *.pdb files from the every <pdb_directory> and every <pdb_file>.\n";
  std::cout << "                                          \tIndex is stored in a file '<index_name>.ind', resp. '<index_name>/residue.ind' if <index_name> is a directory.\n\n";

  std::cout << "-cb\tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "-b\tAll biomolecules and models, but only the first crystallographic transformation are used.\n";
  std::cout << "-c\tAll crystallographic transformations, but only the first biomolecule and model are used.\n";
  std::cout << "\tIf not switch is typed, only the first biomolecule, model and crystallographic transformations are used.\n\n";

  std::cout << "Format of lines in index files is '<residue_id>[\t<chain_id>[\t<pdb_file_basename>]]'.\n";
  std::cout << "Each residue is on separate line and the number of residue is the number of line where it is. The first index is 1.\n";
  std::cout << "Each identifier (chain, complex) also applies to the following lines until it is redefined.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 3;
  const char* arg[] = { argv[0], "C:\\Inspire\\test\\", "C:\\pdb\\test\\" };
  argv = arg;
#endif // TESTING

  if (argc < 3 && (argc < 4 || std::strlen(argv[2]) == 0 || argv[2][0] == '-')) {
    help();
    return 0;
  }

  size_t start;
  inspire::backend::ProteinIterator* it;
  if (std::strlen(argv[2]) > 1 && argv[2][0] == '-') {
    if (argv[2] == "-c") {
      it = new inspire::backend::FirstModelCrystallographicIterator();
    } else if (argv[2] == "-b") {
      it = new inspire::backend::BiomoleculesIterator();
    } else if (argv[2] == "-bc") {
      it = new inspire::backend::AllExceptAltLocIterator();
    } else {
      std::cerr << "ERROR: Modifier '" << argv[2] << "' is not currently supported." << std::endl;
      help();
      return 1;
    }
    start = 3;
  } else {
    it = new inspire::backend::FirstModelIterator();
    start = 2;
  }
  inspire::frontend::Indexer indexer(argv[1], it);

  for (size_t i = start; i < argc; i++) {
    if (elemental::filesystem::is_directory(argv[i])) {
      elemental::filesystem::RecursiveDirectoryFileIterator file_iterator(argv[i]);
      if (file_iterator.has_file()) {
        do {
          indexer.process(file_iterator.filename());
        } while (file_iterator.has_next());
      } else {
        std::cerr << "There is no file in the given directory." << std::endl;
      }
    } else {
      indexer.process(argv[i]);
    }
  }
  
  return 0;
}

