// inspire.cpp : Defines the entry point for the console application.
//

#include "../backend/index.h"
#include "../backend/pdb.h"
#include "../backend/iterators.h"
#include "../elemental/filesystem.h"
#include "../elemental/string.h"
#include <iostream>

//#define TESTING

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Create index of residues in the knowledge-base.\n";
  std::cout << "<index_name> [-(b|c|bc)] (<pdb_directory>|<pdb_file>)+\tIndex all PDB files from the every <pdb_directory> and every <pdb_file>.\n";
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
  const char* arg[] = { argv[0], "C:\\Inspire\\diff\\", "C:\\pdb\\" };
  argv = arg;
#endif // TESTING

  if (argc < 3 && (argc < 4 || std::strlen(argv[2]) == 0 || argv[2][0] == '-')) {
    help();
    return 0;
  }

  size_t start;
  inspire::backend::ProteinIterator* it;
  if (std::strlen(argv[2]) > 1 && argv[2][0] == '-') {
    if (argv[2] == std::string("-c")) {
      it = new inspire::backend::FirstModelCrystallographicIterator();
    } else if (argv[2] == std::string("-b")) {
      it = new inspire::backend::BiomoleculesIterator();
    } else if (argv[2] == std::string("-bc")) {
      it = new inspire::backend::AllExceptAltLocIterator();
    } else if (argv[1] == std::string("-w")) {
      it = new inspire::backend::ExplicitIterator();
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
  inspire::backend::BasicFilter* filter = new inspire::backend::BasicFilter();

  inspire::backend::Indexer indexer(argv[1], it, filter);

  for (size_t i = start; i < argc; i++) {
    indexer.process(argv[i]);
  }

  delete filter;
  delete it;

  return 0;
}

