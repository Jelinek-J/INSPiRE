// inspire.cpp : Defines the entry point for the console application.
//

#include "../backend/index.h"
#include "../backend/pdb.h"
#include "../backend/iterators.h"
#include "../common/filesystem.h"
#include "../common/string.h"
#include <iostream>

//#define TESTING

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Create index of all residues in selected protein files.\n\n";

  std::cout << "Usage:\t<OUTPUT-PATH> [-b|-c|-bc|-w] (<PROTEINS-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<PROTEINS-PATH>     \tPath to a protein or a directory with proteins that should be indexed.\n";
  std::cout << "        \t                    \tSupported file formats are PDB('*.pdb'), PDBx/mmCIF('*.cif') and ('*.xml').\n";
  std::cout << "        \t<OUTPUT-PATH> [-s]  \tWhere to store output file.\n";
  std::cout << "        \t                    \tIf <OUTPUT-PATH> is empty or ends with a directory separator, 'residue.ind' is used as the file name.\n";
  std::cout << "        \t                    \tIf <OUTPUT-PATH> does not end with '.ind' extension, the extension is appended.\n";
  std::cout << "        \t-h                  \tShow informations about the program\n";
  std::cout << "    Iterators: If no iterator is specified, only the first biomolecule from the first model with the first crystallographic transformation is used\n";
  std::cout << "        \t-b                  \tAll biomolecules and models, but only the first crystallographic transformation are used\n";
  std::cout << "        \t-c                  \tAll crystallographic transformations, but only the first biomolecule and model are used\n";
  std::cout << "        \t-bc                 \tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "        \t-w                  \tIgnore both biomolecules and crystallographic transformation, use all chains as they are\n\n";

  std::cout << "Index File Format:\n";
  std::cout << "\tOnly the most probable alternative location is used for each atom. (It can be changed in some future iterators.)\n";
  std::cout << "\tEach residue is on separate line and the number of the residue is the number of its line.\n";
  std::cout << "\tThe index of the first line is 1.\n";
  std::cout << "\tFor each residue, index file contains residue identifier, chain identifier, model identifier and protein identifier in this order and separated by tab.\n";
  std::cout << "\tIf any identifier except residue identifier is missing, it is the same as in the previous line and trailing tabs are trimmed.\n\n";

  std::cout << "Notes:\tOnly the most probable alternative location is used for each atom. (It can be changed in some future iterators.)\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 4;
  const char* arg[] = { argv[0], "C:\\Inspire\\test\\", "-b", "C:\\pdb\\pdb\\multimer\\selected\\" };
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
    } else if (argv[2] == std::string("-w")) {
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

