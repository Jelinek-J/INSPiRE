// features.cpp : Defines the entry point for the console application.
//

//#define TESTING
#ifdef TESTING
//#define FREESASA
#endif // TESTING

#include "../backend/features.h"
#include "../common/exception.h"
#ifdef FREESASA
#include "../common/sasa.h"
#endif // FREESASA

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

    static backend::Feature<std::string>* StringFactory(const char **& argv, const size_t argc, size_t & argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator * it);
    static backend::Feature<float>* FloatFactory(const char **& argv, const size_t argc, size_t & argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator * it);
    static backend::Feature<size_t>* IntFactory(const char **& argv, const size_t argc, size_t & argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator * it);
    static void FactoryCheck(const size_t argi, const size_t argc, const char** &argv) {
      if (argi >= argc) {
        throw common::exception::TitledException("Argument #" + std::to_string(argi) + " is empty.");
      }
      if (strlen(argv[argi]) == 0) {
        throw common::exception::TitledException("Argument #" + std::to_string(argi) + " is empty.");
      }
      if (argv[argi][0] != '-') {
        throw common::exception::TitledException("Argument expected at position #" + std::to_string(argi) + " instead of '" + argv[argi] + "'.");
      }
      if (strlen(argv[argi]) == 1) {
        throw common::exception::TitledException("Feature specifier missing at position #" + std::to_string(argi) + ": '" + argv[argi] + "'");
      }
    }

    static backend::Feature<std::string>* StringFactory(const char** &argv, const size_t argc, size_t &argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator* it) {
      FactoryCheck(argi, argc, argv);
      switch (argv[argi][1]) {
        case 's':
        case 't':
        case 'F':
        case 'R':
          {
            inspire::backend::Feature<float>* subfeature = FloatFactory(argv, argc, argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::ToStringFeature<float>(subfeature);
          }
          break;
        case 'B':
          {
            inspire::backend::Feature<size_t>* subfeature = IntFactory(argv, argc, argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::ToStringFeature<size_t>(subfeature);
          }
          break;
        case 'a':
          return new inspire::backend::AminoacidFeature(it);
          break;
        case 'c':
          return new inspire::backend::CoordinateFeature(it);
          break;
        case 'e':
          return new inspire::backend::CompositionFeature(it);
          break;
        case 'i':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Interface type switch miss a definition of translation table");
          }
          {
            std::string arg(argv[argi]);
            size_t space = arg.find_last_of(';');
            if (space != arg.npos && space != arg.size()-1 && arg[arg.size()-1] != '"' && arg[arg.size()-1] != '\'') {
              return new inspire::backend::InterfaceLargeFeature(it, arg.substr(2, space-2), std::stod(arg.substr(space+1)));
            } else {
              return new inspire::backend::InterfaceLargeFeature(it, arg.substr(2), 0.5);
            }
          }
          break;
        case 'p':
          return new inspire::backend::ProteinIdFeature();
          break;
        case 'L':
          {
            std::string arg(argv[argi]);
            size_t first = arg.find(';');
            if (first == arg.npos) {
              throw common::exception::TitledException("String feature loader miss a definition of header to load");
            }
            return new inspire::backend::StringLoaderFeature(arg.substr(2, first-2), arg.substr(first+1));
          }
          break;
        case 'N':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Rename feature switch miss a definition of translation table");
          }
          {
            std::string name = std::string(argv[argi]).substr(2);
            inspire::backend::Feature<std::string>* subfeature = StringFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::RenameFeature<std::string>(name, subfeature);
          }
          break;
        case 'P':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Projection feature switch miss a definition of translation table");
          }
          {
            std::string mapping = std::string(argv[argi]).substr(2);
            inspire::backend::Feature<std::string>* subfeature = StringFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::StringProjectionFeature(mapping, subfeature);
          }
          break;
        case 'X':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Xenofeature switch miss a definition of file format");
          }
          switch (argv[argi][2]) {
            case 's':
              return new inspire::backend::BasicExternalLoaderFeature(it, std::string(argv[argi]).substr(3), "xenofeature");
              break;
            case 'f':
              return new inspire::backend::FullExternalLoaderFeature(it, std::string(argv[argi]).substr(3), "xenofeature");
              break;
            default:
              throw common::exception::TitledException("Unknown file format speciffier specifier: '" + std::string(1, argv[argi][2]) + "'");
              break;
          }
          break;
          break;
        default:
          throw common::exception::TitledException("Unknown feature specifier: '" + std::string(argv[argi]) + "'");
      }
    }

    static backend::Feature<float>* FloatFactory(const char** &argv, const size_t argc, size_t &argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator* it) {
      FactoryCheck(argi, argc, argv);
      switch (argv[argi][1]) {
        case 'a':
        case 'c':
        case 'e':
        case 'i':
        case 'p':
        case 'B':
        case 'L':
        case 'P':
          throw common::exception::TitledException("Unexpected feature that is not convertible to float at position #" + std::to_string(argi) + ": '" + argv[argi] + "'");
          break;
#ifdef FREESASA
        case 's':
          {
            std::string arg(argv[argi]);
            size_t first = arg.find(';');
            if (first == arg.npos) {
              throw common::exception::TitledException("Solvent Accessible Surface Area switch miss a definition of aminoacids composition");
            }
            inspire::backend::Feature<float>* subfeature;
            return new common::sasa::SasaFeature(it, arg.substr(2, first-2), arg.substr(first+1));
          }
          break;
#endif // FREESASA
        case 't':
          return new inspire::backend::TemperatureFeature(it);
          break;
        case 'F':
          {
            std::string arg(argv[argi]);
            size_t first = arg.find(';');
            if (first == arg.npos) {
              throw common::exception::TitledException("Float feature loader miss a definition of header to load");
            }
            return new inspire::backend::FloatLoaderFeature(arg.substr(2, first-2), arg.substr(first+1));
          }
          break;
        case 'N':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Rename feature switch miss a definition of translation table");
          }
          {
            std::string name = std::string(argv[argi]).substr(2);
            inspire::backend::Feature<float>* subfeature = FloatFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::RenameFeature<float>(name, subfeature);
          }
          break;
        case 'R':
          {
            if (strlen(argv[argi]) <= 2) {
              throw common::exception::TitledException("Relative feature miss a definition of reference values");
            }
            std::string arg(argv[argi]);
            inspire::backend::Feature<float>* floatfeature = FloatFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(floatfeature);
            inspire::backend::Feature<std::string>* stringfeature = StringFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(stringfeature);
            return new inspire::backend::RelativeAminoacidFeature(floatfeature, stringfeature, arg.substr(2));
          }
          break;
        default:
          throw common::exception::TitledException("Unknown feature specifier: '" + std::string(argv[argi]) + "'");
      }
    }

    static backend::Feature<size_t>* IntFactory(const char** &argv, const size_t argc, size_t &argi, std::vector<inspire::backend::IFeature*> inner_features, inspire::backend::ProteinIterator* it) {
      FactoryCheck(argi, argc, argv);
      switch (argv[argi][1]) {
        case 'a':
        case 'c':
        case 'e':
        case 'i':
        case 'p':
        case 's':
        case 't':
        case 'F':
        case 'L':
        case 'P':
        case 'R':
          throw common::exception::TitledException("Unexpected feature that is not convertible to size_t at position #" + std::to_string(argi) + ": '" + argv[argi] + "'");
          break;
        case 'B':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Binning feature switch miss a definition of intervals");
          }
          {
            std::string intervals = std::string(argv[argi]).substr(2);
            inspire::backend::Feature<float>* subfeature = FloatFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::BinFeature(subfeature, intervals);
          }
          break;
        case 'N':
          if (strlen(argv[argi]) == 2) {
            throw common::exception::TitledException("Rename feature switch miss a definition of translation table");
          }
          {
            std::string name = std::string(argv[argi]).substr(2);
            inspire::backend::Feature<size_t>* subfeature = IntFactory(argv, argc, ++argi, inner_features, it);
            inner_features.push_back(subfeature);
            return new inspire::backend::RenameFeature<size_t>(name, subfeature);
          }
          break;
        default:
          throw common::exception::TitledException("Unknown feature specifier: '" + std::string(argv[argi]) + "'");
      }
    }
  }
}

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Extracts required features from proteins indexed in a given index file.\n\n";

  std::cout << "Usage:\t[-b|-c|-bc|-w] <INDEX-FILE> <OUTPUT-PATH> [-] <FEATURE>+ (<PROTEINS-PATH>)+\n";
  std::cout << "      \t-h\n\n";

  std::cout << "Options:\t<INDEX-FILE>     \tPath to a index file\n";
  std::cout << "        \t-                \tSave each feature in a separate file (otherwise all features are stored together in one file)\n";
  std::cout << "        \t<OUTPUT-PATH>    \tWhere to store output file.\n";
  std::cout << "        \t                 \tIf '-' is typed, <OUTPUT-PATH> must be a directory;\n";
  std::cout << "        \t                 \tand each feature is stored in a separated file named according to the corresponding feature with '.tur' as an extension.\n";
  std::cout << "        \t                 \tOtherwise(if '-' is not typed), all features are stored in the same file.\n";
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
  std::cout << "        \t<FEATURE> can be following (lowercase switcher means the feature is loaded directly from proteins,\n";
  std::cout << "        \t                            while uppercase switcher means the feature transforms informations from other features/ files):\n";
  std::cout << "        \t-a  \tAmino acid type three letter codes\n";
  std::cout << "        \t-c  \tCoordinates of Carbon alpha of a residue separated by space\n";
  std::cout << "        \t-e  \tSorted atom names of all atoms in a residue except helium and deuterium separated by space\n";
  std::cout << "        \t-i<RADII-FILE>[;<DISTANCE>]\n";
  std::cout << "        \t    \tWhether a residue is an interfacial residue with <RADII-FILE> defining radiuses of chemical elements and\n";
  std::cout << "        \t    \t<DISTANCE> sets the maximal allowed distance of two van der Waals radiuses (0.5A is a default value).\n";
#ifdef FREESASA
  std::cout << "        \t-s<RADII-FILE>;<COMPOSITION-FILE>\n";
  std::cout << "        \t    \tRelative solvent accessible surface area with residues' composition defined in <COMPOSITION-FILE>, atomic radiuses defined in <RADII-FILE>.\n";
#endif // FREESASA
  std::cout << "        \t-p  \tProtein's identifier (usefull e.g. to make temperature factor comparable between proteins)\n";
  std::cout << "        \t-t  \tTemperature factor of a residuen\n";
  std::cout << "        \t-B<BOUNDARIES-FILE> <FLOAT-FEATURE>\n";
  std::cout << "        \t    \tBinned values of <FLOAT-FEATURE> based on splitting points in <BOUDARIES-FILE>.\n";
  std::cout << "        \t    \tI.e. in the case of sorted <BOUDARIES-FILE> it returns <i> if the value is greater than ith splitting point but lower than or equal to (<i>+1)th splitting point.\n";
  std::cout << "        \t    \t(And 0 if it is lower than to equal to the first splitting point and (count of splitting points)+1 if it is greater than the last splitting point.).\n";
  std::cout << "        \t-F<FEATURE-FILE>;<HEADER>\n";
  std::cout << "        \t    \tLoad a float feature <HEADER> from a file <FEATURE-FILE>\n";
  std::cout << "        \t-L<FEATURE-FILE>;<HEADER>\n";
  std::cout << "        \t    \tLoad a string feature <HEADER> from a file <FEATURE-FILE>\n";
  std::cout << "        \t-N<NEW-TITLE> <FEATURE>\n";
  std::cout << "        \t    \tChange title of the inner feature (does not change values of the feature)\n";
  std::cout << "        \t-P<PROJECTION-FILE> <STRING-FEATURE>\n";
  std::cout << "        \t    \tTransforms feature <STRING-FEATURE> based on a dictionary defined in a <PROJECTION-FILE> (missing keys are skipped).\n";
  std::cout << "        \t    \tE.g. to transform amino acids' three-letter codes to one-letter codes.\n";
  std::cout << "        \t    \tDictionary in the <PROJECTION-FILE> should be in format 'key<TAB>value'.\n";
  std::cout << "        \t-R<REFERENCE-VALUES> <FLOAT-FEATURE> <STRING-FEATURE>\n";
  std::cout << "        \t    \tFor each residue relativizes a value of <FLOAT-FEATURE> based on reference value in file <REFERENCE-VALUES> for corresponding value of <STRING-FEATURE>.\n";
  std::cout << "        \t    \tE.g. to transform solvent accessible surface area to relative solvent accessible surface area.\n";
  std::cout << "        \t-X(s|f)<FEATURE-FILE>\n";
  std::cout << "        \t    \tLoad an external feature from <FEATURE-FILE> file in format <residue_id>\t<value>.\n";
  std::cout << "        \t    \tThe feature gets 'xenofeature' as a title (it can be changed by '-N' feature).\n";
  std::cout << "        \t    \t's'/ 'f' defines a format of the <residue_id> - 's' means a simple format, while 'f' means a full format.\n";
  std::cout << "        \t    \tThe simple format has <residue_id> in a form '<protein_id>.<chain_id>.<residue_number><insertion_code>' and\n";
  std::cout << "        \t    \tvalues will by copied to all models, biomolecules, chains transformed by biomolecule/ crystallographic transformations, etc.\n";
  std::cout << "        \t    \tThe full format depends on the used iterator and its identifiers used in the index file:\n";
  std::cout << "        \t    \t    for the default iterator: '<protein_id>.<chain_id>+<assemblyTransformationID>.<residue_number><insertion_code>';\n";
  std::cout << "        \t    \t    for the iterator '-w':    '<protein_id>.<chain_id>.<residue_number><insertion_code>';\n";
  std::cout << "        \t    \t    for the iterator '-b':    '<protein_id>;<biomoleculeID>,<modelID>.<chain_id>+<assemblyTransformationID>.<residue_number><insertion_code>';\n";
  std::cout << "        \t    \t    for the iterator '-c':    '<protein_id>.<chain_id>+<assemblyTransformationID>*<crystallographicTransformationID>.<residue_number><insertion_code>';\n";
  std::cout << "        \t    \t    for the iterator '-bc':   '<protein_id>;<biomoleculeID>,<modelID>.<chain_id>+<assemblyTransformationID>*<crystallographicTransformationID>.<residue_number><insertion_code>'.\n\n";

  std::cout << "Feature Files Format:\n";
  std::cout << "\tHeader line:\tThe first line of each feature file; names of columns (features) are separated by a tabulator.\n";
  std::cout << "\tData line:  \tEach feature is separated by a tabulator, if it is not possible to extract an feature, the corresponding cell is empty.\n";
  std::cout << "\t            \tDefaultly, index of each line is an index of the previous line plus one (the first line has index '1').\n";
  std::cout << "\t            \tIf a line has not consecutive index, the index is specified in an extra(last) column.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  argc = 8;
  const char* arg[] = {argv[0], "-b", "C:\\Inspire\\test\\query\\residues.ind", "C:\\Inspire\\test\\query\\", "-", "-PC:\\Inspire\\test\\query\\aminoacid.nor", "-LC:\\Inspire\\test\\query\\aminoacid.tur;aminoacid", "C:\\Inspire\\pdb\\query\\"};
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
  try {
    for (size_t i = start; i < stop; i++) {
      features.push_back(inspire::frontend::StringFactory(argv, stop, i, inner_features, it));
      // TODO: Consider an opposite approach extracting features parallel into multiple files as it saves reading of files
      if (separate) {
        try {
          std::string file = path + features[0]->title() + ".tur";
          extractor.extract(file, features);
          delete features[0];
          features.clear();
        } catch (common::exception::TitledException e) {
          std::cerr << e.what() << std::endl;
        }
      }
    }
    if (!separate) {
      extractor.extract(path, features);
    }
  } catch (common::exception::TitledException e) {
    std::cerr << e.what() << std::endl;
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
