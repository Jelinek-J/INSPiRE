// features.cpp : Defines the entry point for the console application.
//

//#define FREESASA

#include "../backend/features.h"
#include "../common/exception.h"
#ifdef FREESASA
#include "../common/sasa.h"
#endif // FREESASA

//#define TESTING

namespace inspire {
  namespace frontend {
    class Features {
#ifdef TESTING
      // Errors log for testing reasons
      std::ofstream LOG;
#endif // TESTING
      backend::Features FEATURES;

      public:
      Features(std::string index, backend::ProteinIterator* it, backend::BasicFilter* filter) :
#ifdef TESTING
        LOG("C:\\Inspire\\error.log"),
#endif // TESTING
        FEATURES(index, it, filter) {
      }

      ~Features() {
#ifdef TESTING
        LOG.flush();
        LOG.close();
#endif // TESTING
      }

      void index(std::string file) {
#ifdef TESTING
        std::cout << "Indexing: " << file << '\r';
#endif // DEBUG
        try {
          FEATURES.index_pdb(file);
        } catch (const common::exception::TitledException& e) {
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
#ifdef TESTING
        std::cout << "Indexed: " << file << std::endl;
#endif // DEBUG
      }

      void extract(std::string &file, std::vector<backend::Feature<std::string>*> features) {
#ifdef TESTING
        std::cout << "Extracting: " << file << '\r';
#endif // DEBUG
        if (common::string::ends_with(file, ".tur")) {
          try {
            FEATURES.extract_features(file, features);
          } catch (const common::exception::TitledException& e) {
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
#ifdef TESTING
        std::cout << "Extracted: " << file << std::endl;
#endif // DEBUG
      }
    };
  }
}

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Extracts required features from proteins indexed in a given index file.\n\n";

  std::cout << "Usage:\t[-b|-c|-bc|-w] <INDEX-FILE> <OUTPUT-PATH> [-] (-a<TRANSFORMATION-FILE>|-c|-e|-i<RADII-FILE>[<DISTANCE>]";
#ifdef FREESASA
  std::cout << "|-r<RADII-FILE>;<COMPOSITION-FILE>;<MAX-SASA-FILE>";
#endif // FREESASA
  std::cout << "|-t)+ (<PROTEINS-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INDEX-FILE>     \tPath to a index file\n";
  std::cout << "        \t-               \tSave each feature in a separate file (otherwise all features are stored together in one file)\n";
  std::cout << "        \t<OUTPUT-PATH>    \tWhere to store output file.\n";
  std::cout << "        \t                 \tIf '-s' is typed, <OUTPUT-PATH> must be a directory; and each feature is stored in a separated file named according to the corresponding feature with '.tur' as an extension.\n";
  std::cout << "        \t                 \tOtherwise(if '-s' is not typed), all features are stored in the same file.\n";
  std::cout << "        \t                 \t  If <OUTPUT-PATH> is a directory or ends with a directory separator, 'features.tur' is used as the file name.\n";
  std::cout << "        \t                 \t  If <OUTPUT-PATH> does not end with '.tur' extension, the extension is appended.\n";
  std::cout << "        \t<PROTEINS-PATH>  \tPath to a protein or a directory with proteins\n";
  std::cout << "        \t-h, --help       \tShow informations about the program\n";
  std::cout << "    Iterators: (if no iterator is specified, only the first biomolecule from the first model with the first crystallographic transformation is used)\n";
  std::cout << "        \t-b   \tAll biomolecules and models, but only the first crystallographic transformation are used\n";
  std::cout << "        \t-c   \tAll crystallographic transformations, but only the first biomolecule and model are used\n";
  std::cout << "        \t-bc  \tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "        \t-w   \tIgnore both biomolecules and crystallographic transformation, use all chains as they are\n";
  std::cout << "    Features:\n";
  std::cout << "        \t-a<TRANSFORMATION-FILE>\n";
  std::cout << "        \t    \tAminoacid type three-letter code is transformed with transformation defined in <TRANSFORMATION-FILE>, that should be in format 'key\\tvalue'\n";
  std::cout << "        \t-c  \tCoordinates of carbon_alpha of an aminoacid\n";
  std::cout << "        \t-e  \tAtomic composition of a residue (helium and deuterium are skipped)\n";
  std::cout << "        \t-i<RADII-FILE>[<DISTANCE>]\n";
  std::cout << "        \t    \tWhether a residue is an interfacial residue with <RADII-FILE> defining radiuses of chemical elements and\n";
  std::cout << "        \t    \t<DISTANCE> sets the maximal allowed distance of two van der Waals radiuses (0.5A is a default value).\n";
  std::cout << "        \t    \t<DISTANCE> must be separated by a space from <RADII-FILE>\n";
#ifdef FREESASA
  std::cout << "        \t-r<RADII-FILE>;<COMPOSITION-FILE>;<MAX-SASA-FILE>\n";
  std::cout << "        \t    \tRelative solvent accessible surface area with residues' composition defined in <COMPOSITION-FILE>, atomic radiuses defined in <RADII-FILE> and\n";
  std::cout << "        \t    \treference solvent accessible surface areas defined in <MAX-SASA-FILE>.\n";
  std::cout << "        \t-s<RADII-FILE>;<COMPOSITION-FILE>\n";
  std::cout << "        \t    \tRelative solvent accessible surface area with residues' composition defined in <COMPOSITION-FILE>, atomic radiuses defined in <RADII-FILE>.\n";
#endif // FREESASA
  std::cout << "        \t-t  \tTemperature factor of an aminoacid\n\n";

  std::cout << "Feature Files Format:\n";
  std::cout << "\tHeader line:\tThe first line of each feature file; names of columns (features) are separated by a tabulator.\n";
  std::cout << "\tData line:  \tEach feature is separated by a tabulator, if it is not possible to extract an feature, the corresponding cell is empty.\n";
  std::cout << "\t            \tDefaultly, index of each line is an index of the previous line plus one (the first line has index '1').\n";
  std::cout << "\t            \tIf a line has not consecutive index, the index is specified in an extra(last) column.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 7;
  const char* arg[] = {argv[0], "-b", "C:\\Inspire\\test\\2760\\example\\construction\\residues.ind", "C:\\Inspire\\test\\2760\\example\\", "-s", "-rC:\\Inspire\\test\\2760\\example\\radiuses.rus;C:\\Inspire\\test\\2760\\example\\composition.cop", "C:\\Inspire\\pdb\\5j7v.pdb"};
  argv = arg;
#endif // TESTING

  if (argc < 5 && (argc < 6 || std::strlen(argv[2]) == 0 || argv[2][0] == '-')) {
    help();
    return 0;
  }

  size_t start;
  inspire::backend::ProteinIterator* it;
  if (std::strlen(argv[1]) > 1 && argv[1][0] == '-') {
    if (argv[1] == std::string("-c")) {
      it = new inspire::backend::FirstModelCrystallographicIterator();
    } else if (argv[1] == std::string("-b")) {
      it = new inspire::backend::BiomoleculesIterator();
    } else if (argv[1] == std::string("-bc")) {
      it = new inspire::backend::AllExceptAltLocIterator();
    } else if (argv[1] == std::string("-w")) {
      it = new inspire::backend::ExplicitIterator();
    } else {
      std::cerr << "ERROR: Modifier '" << argv[1] << "' is not currently supported." << std::endl;
      help();
      return 1;
    }
    start = 3;
  } else {
    it = new inspire::backend::FirstModelIterator();
    start = 2;
  }
  inspire::backend::BasicFilter* filter = new inspire::backend::BasicFilter();

  inspire::frontend::Features extractor(argv[start-1], it, filter);

  size_t stop = start+2;
  while (stop < argc && (strlen(argv[stop]) == 0 || argv[stop][0] == '-')) {
    ++stop;
  }
  if (stop == argc) {
    std::cerr << "Missing <pdb_directory> or <pdb_file> argument." << std::endl;
    help();
    return 0;
  }
  for (size_t i = stop; i < argc; i++) {
    extractor.index(argv[i]);
  }

  bool separate;
  std::string path(argv[start]);
  if (argv[start+1] == std::string("-")) {
    if (!common::filesystem::is_directory(argv[start])) {
      std::cerr << "'" << argv[start] << "' is not a directory or does not exists." << std::endl;
      return 5;
    }
    if (strlen(argv[start+1]) == 0) {
      std::cerr << "Argument #" << start+1 << " is empty." << std::endl;
      return 8;
    }
    if (argv[start+1][0] != '-') {
      std::cerr << "No feature is selected for extraction." << std::endl;
      return 6;
    }

    if (path.size() > 0 && path[path.size()-1] != common::filesystem::directory_separator) {
      path.push_back(common::filesystem::directory_separator);
    }
    separate = true;
    start += 2;
  } else {
    if (path.empty() || path.back() == common::filesystem::directory_separator) {
      path += "features.tur";
    } else if (!common::string::ends_with(path, ".tur")) {
      path += ".tur";
    }
    separate = false;
    ++start;
  }

  std::vector<inspire::backend::Feature<std::string>*> features;
  std::vector<inspire::backend::IFeature*> inner_features;
  for (size_t i = start; i < stop; i++) {
    if (strlen(argv[i]) == 0) {
      std::cerr << "Argument #" << i << " is empty." << std::endl;
      continue;
    }
    if (argv[i][0] != '-') {
      std::cerr << "Argument expected at position #" << i << " instead of '" << argv[i] << "'." << std::endl;
      continue;
    }
    if (strlen(argv[i]) == 1) {
      std::cerr << "Feature specifier missing at position #" << i << ": '" << argv[i] << "'" << std::endl;
      continue;
    }
    switch (argv[i][1]) {
      case 'a':
        if (strlen(argv[i]) == 2) {
          std::cerr << "Aminoacid type switch miss a definition of translation table" << std::endl;
          continue;
        }
        {
          inspire::backend::Feature<std::string>* feature = new inspire::backend::SingleAminoacidFeature(it, std::string(argv[i]).substr(2));
          features.push_back(feature);
        }
        break;
      case 'c':
      {
        inspire::backend::Feature<std::string>* feature = new inspire::backend::CoordinateFeature(it);
        features.push_back(feature);
      }
      break;
      case 'e':
      {
        inspire::backend::Feature<std::string>* feature = new inspire::backend::CompositionFeature(it);
        features.push_back(feature);
      }
      break;
      case 'i':
        if (strlen(argv[i]) == 2) {
          std::cerr << "Aminoacid type switch miss a definition of translation table" << std::endl;
          continue;
        }
        {
          std::string arg(argv[i]);
          size_t space = arg.find_last_of(' ');
          inspire::backend::Feature<std::string>* feature;
          if (space != arg.npos && space != arg.size()-1 && arg[arg.size()-1] != '"' && arg[arg.size()-1] != '\'') {
            feature = new inspire::backend::InterfaceLargeFeature(it, arg.substr(2, space-2), std::stod(arg.substr(space+1)));
          } else {
            feature = new inspire::backend::InterfaceLargeFeature(it, arg.substr(2), 0.5);
          }
          features.push_back(feature);
        }
        break;
#ifdef FREESASA
      case 'r':
      {
        std::string arg(argv[i]);
        size_t first = arg.find(';');
        size_t second = arg.rfind(';');
        if (first == arg.npos) {
          std::cerr << "Solvent Accessible Surface Area switch miss a definition of aminoacids composition" << std::endl;
          continue;
        }
        if (second == first) {
          std::cerr << "Solvent Accessible Surface Area switch miss a definition of maximal sasas" << std::endl;
          continue;
        }
        inspire::backend::Feature<float>* subfeature;
        subfeature = new common::sasa::SasaFeature(it, arg.substr(2, first-2), arg.substr(first+1, second-first-1));
        inner_features.push_back(subfeature);
        subfeature = new inspire::backend::RelativeAminoacidFeature(it, subfeature, arg.substr(second+1));
        inner_features.push_back(subfeature);
        inspire::backend::Feature<std::string>* feature;
        feature = new inspire::backend::ToStringFeature<float>(it, subfeature);
        features.push_back(feature);
      }
      break;
      case 's':
      {
        std::string arg(argv[i]);
        size_t first = arg.find(';');
        if (first == arg.npos) {
          std::cerr << "Solvent Accessible Surface Area switch miss a definition of aminoacids composition" << std::endl;
          continue;
        }
        inspire::backend::Feature<float>* subfeature;
        subfeature = new common::sasa::SasaFeature(it, arg.substr(2, first-2), arg.substr(first+1));
        inner_features.push_back(subfeature);
        inspire::backend::Feature<std::string>* feature;
        feature = new inspire::backend::ToStringFeature<float>(it, subfeature);
        features.push_back(feature);
      }
      break;
#endif // FREESASA
      case 't':
        {
          inspire::backend::Feature<std::string>* feature = new inspire::backend::TemperatureFeature(it);
          features.push_back(feature);
        }
        break;
      default:
        std::cerr << "Unknown feature specifier: '" << argv[i] << "'" << std::endl;
        continue;
        break;
    }
    // TODO: Consider an opposite approach extracting features parallel into multiple files as it saves reading of files
    if (separate) {
      std::string file = path + features[0]->title() + ".tur";
      extractor.extract(file, features);
      delete features[0];
      features.clear();
    }
  }
  if (!separate) {
    extractor.extract(path, features);
  }
  for (auto features_it = features.begin(); features_it != features.end(); ++features_it) {
    delete *features_it;
  }
  for (auto features_it = inner_features.begin(); features_it != inner_features.end(); ++features_it) {
    delete *features_it;
  }
  delete filter;
  delete it;

  return 0;
}
