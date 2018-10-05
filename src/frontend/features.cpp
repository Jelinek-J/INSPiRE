// features.cpp : Defines the entry point for the console application.
//

#include "../backend/features.h"
#include "../elemental/exception.h"

namespace inspire {
  namespace frontend {
    class Features {
      backend::Features FEATURES;

      public:
      Features(std::string index, backend::ProteinIterator* it, backend::BasicFilter* filter) :
        FEATURES(index, it, filter) {
      }

      void index(std::string file) {
        try {
          FEATURES.index_pdb(file);
        } catch (const elemental::exception::TitledException& e) {
          std::cerr << "ERROR: " << e.what() << std::endl;
        } catch (const std::exception& e) {
          std::cerr << "ERROR: " << e.what() << std::endl;
        } catch (...) {
          std::cerr << "UNKNOWN ERROR" << std::endl;
        }
      }

      void extract(std::string &file, std::vector<backend::Feature*> features) {
        if (elemental::string::ends_with(file, ".tur")) {
          try {
            FEATURES.extract_features(file, features);
          } catch (const elemental::exception::TitledException& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
          } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
          } catch (...) {
            std::cerr << "UNKNOWN ERROR" << std::endl;
          }
        }
      }
    };
  }
}

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Extract a required features for the given index file.\n";
  std::cout << "[-(b|c|bc)] <index_path> <features_path> [-s] (-a<transformation_path>|-c|-i<radiuses_path>[<distance>]|-t)+ (<pdb_directory>|<pdb_file>)+\tIndex all features from proteins indexed in file located in <index_path> and\n";
  std::cout << "                                                                             \tstored in any *.pdb file from <pdb_directory> or listed as <pdb_file>.\n\n";

  std::cout << "Iterator used for indexing files:\n";
  std::cout << "-cb\tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "-b\tAll biomolecules and models, but only the first crystallographic transformation are used.\n";
  std::cout << "-c\tAll crystallographic transformations, but only the first biomolecule and model are used.\n";
  std::cout << "\tIf not switch is typed, only the first biomolecule, model and crystallographic transformations are used.\n\n";

  std::cout << "Recognized features:\n";
  std::cout << "-a<transformation_path>      \tAminoacid type is transformed with transformation defined in file located <transformation_path> in format 'key<TAB>value'";
  std::cout << "-c                           \tCoordinates of carbon_alpha of an aminoacid";
  std::cout << "-i<radiuses_path>[<distance>]\tWhether an aminoacid is an interfacial aminoacid with file <radiuses_path> defining radiuses of chemical elements and <distance> settings the maximal allowed distance of two van der Waals radiuses (0.5A is a default value)";
//  std::cout << "-r                           \tRelative solvent accessible surface area";
  std::cout << "-s                           \tStore each feature in a separated file, in that case, <features_path> must be an directory and each feature will be stored in file '<feature_name>.tur'";
  std::cout << "-t                           \tTemperature factor of an aminoacid";

  std::cout << "Format of lines in feature files is:\n";
  std::cout << "Header line: The first line of each feature file; names of columns (features) are separated by a tabulator.\n";
  std::cout << "Data line: Each feature is separated by a tabulator, if it is not possible to extract an feature, the corresponding cell is empty.\n";
  std::cout << "           Defaultly, index of each line is an index of the previous line plus one (the first line has index '1'),\n";
  std::cout << "           Except the line has one extra cell with a number - that set the index of that line.\n\n";
}

int main(int argc, const char** argv) {
  if (argc < 5 && (argc < 6 || std::strlen(argv[2]) == 0 || argv[2][0] == '-')) {
    help();
    return 0;
  }

  size_t start;
  inspire::backend::ProteinIterator* it;
  if (std::strlen(argv[1]) > 1 && argv[1][0] == '-') {
    if (argv[1] == "-c") {
      it = new inspire::backend::FirstModelCrystallographicIterator();
    } else if (argv[1] == "-b") {
      it = new inspire::backend::BiomoleculesIterator();
    } else if (argv[1] == "-bc") {
      it = new inspire::backend::AllExceptAltLocIterator();
    } else if (argv[1] == "-w") {
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
  if (argv[start+1] == "-s") {
    if (!elemental::filesystem::is_directory(argv[start])) {
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

    if (path.size() > 0 && path[path.size()-1] != elemental::filesystem::directory_separator) {
      path.push_back(elemental::filesystem::directory_separator);
    }
    separate = true;
    start += 2;
  } else {
    if (path.empty() || path.back() == elemental::filesystem::directory_separator) {
      path += "features.tur";
    } else if (!elemental::string::ends_with(path, ".tur")) {
      path += ".tur";
    }
    separate = false;
    ++start;
  }

  std::vector<inspire::backend::Feature*> features;
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
          inspire::backend::Feature* feature = new inspire::backend::BasicAminoacidFeature(it, std::string(argv[i]).substr(2));
          features.push_back(feature);
        }
        break;
      case 'c':
        {
          inspire::backend::Feature* feature = new inspire::backend::CoordinateFeature(it);
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
          inspire::backend::Feature* feature;
          if (space != arg.npos && space != arg.size()-1 && arg[arg.size()-1] != '"' && arg[arg.size()-1] != '\'') {
            feature = new inspire::backend::InterfaceLargeFeature(it, arg.substr(2, space-2), std::stod(arg.substr(space+1)));
          } else {
            feature = new inspire::backend::InterfaceLargeFeature(it, arg.substr(2), 0.5);
          }
          features.push_back(feature);
        }
        break;
      //case 'r':
      //  std::cerr << "Feature RASA is not implemented yet" << std::endl;
      //  continue;
      //  break;
      case 't':
      {
        inspire::backend::Feature* feature = new inspire::backend::TemperatureFeature(it);
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
      features.clear();
    }
  }
  if (!separate) {
    extractor.extract(path, features);
  }

  return 0;
}
