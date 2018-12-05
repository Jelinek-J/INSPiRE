#pragma once

#include "../common/multithread.h"
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <sstream>
#include <thread>

namespace inspire {
  namespace backend {
    class Mine {
      private:
      // Precomputed number of bits set to one.
      std::map<int, int> DISTANCES;

      // Number of parallel threads for data mining
      char THREADS;
      // Number of most similar elements to mine
      int LIMIT;

      // Order of features for prefilter fingerprints according to central residue
      const std::set<std::string> FILTERS;
      // Fingerprints and indices of coresponding structural elements binned according to ordered features of central residue
      std::unordered_map<std::string, std::unordered_map<std::string, std::vector<uint32_t> > > KNOWLEDGE_BASE;
      // What structural elements should be ignored to prevent overlearning during training (e.g. they are from the same protein)
      std::unordered_map<int, std::unordered_set<int>*> EXCLUDE;
      
      // List of queries
      // NOTE: Compare with sync_bounded_queue - should occupy less memory, on the other hand it require load file in temporary collection to find out its size and thus is slower.
      // NOTE: Of course, it could be JSON parser, but this approach allows a parallel parsing.
      // NOTE: Well, lines could be readed online when necessary, but it would result in delays caused by HDD characteristics.
      common::multithread::synchronized_queue TASKS;


      // Check whether file path contains all requered features and check its extension.
      bool kb_path_check(const std::string& file, const std::set<std::string> &filters) {
        if (!common::string::ends_with(file, ".fin")) {
          return false;
        }
        for (auto filters_it = filters.begin(); filters_it != filters.end(); ++filters_it) {
          if (file.find(std::string(1, boost::filesystem::path::preferred_separator) + *filters_it + std::string(1, '_')) == file.npos) {
            return false;
          }
        }
        return true;
      }

      // Create unique key from features to ensure it konsistency between knowledge-base loading and prediction
      std::string create_key(std::map<std::string, std::string>& features) {
        //? Will be some 'stringbuilder' e.g. std::stringstream significantly more effective?
        std::string ret;
        for (auto features_it = features.begin(); features_it != features.end(); ++features_it) {
          if (!ret.empty()) {
            ret += ".";
          }
          ret += features_it->second;
        }
        return ret;
      }

      // Load queries
      void load_tasks(std::string path) {
        std::ifstream input(path);
        std::string line;
        while (std::getline(input, line)) {
          TASKS.push(line);
        }
      }

      std::vector<uint32_t>* filter_siblings(int id, std::vector<uint32_t>* fingerprints) {
        // Sibilings are not defined.
        if (EXCLUDE.size() == 0) {
          return fingerprints;
        }

        //? Initialize buffer for fingerprints.size() length - some indices can be filtered out, but it is just temporary structure, so resizing could lead to worse performance
        std::vector<uint32_t>* ret = new std::vector<uint32_t>();
        //? 'for each' is not C++ construct - it is only in Visual C++, right?
        for (std::vector<uint32_t>::iterator it = fingerprints->begin(); it != fingerprints->end(); it++) {
          auto siblings = EXCLUDE.find(id);
          if (siblings == EXCLUDE.end() || siblings->second->find(*it) == siblings->second->end()) {
            ret->push_back(*it);
          }
        }
        //? How does it work, does it also create a new copy of it for the caller, or does it return exactly this one instance like in C# etc.?
        return ret;
      }

      void selectThread(std::string output) {
        std::ofstream stream(output);
        std::string line;
        while (TASKS.try_pull(line) == boost::queue_op_status::success) {
          int id;
          std::string key;
          //? How much effective is to initialize with size? It is expected that all fingerprints have the same size but in principle every fingerprint could have a different length.
          std::string fingerprint;
          {
            std::stringstream stream(line);

            // Load id of task
            if (!std::getline(stream, key, '\t')) { continue; }
            id = std::stoi(key);

            // Load features of central residue
            std::map<std::string, std::string> features;
            //? Will be "value" as constant for whole program or it will be created every iteration?
            while (std::getline(stream, key, ':') && !stream.eof()) {
              std::string value;
              std::getline(stream, value, '\t');
              std::set<std::string>::iterator it = FILTERS.find(key);
              if (it != FILTERS.end()) {
                // TODO: Binning of the feature. Maybe JSON file for creation of fingerprints can be reused. Just optional.
                features.insert({key, value});
              }
            }

            // Load fingerprint
            char current = 0;
            for (size_t i = 0; i < key.size(); i++) {
              if (key[i] == '1') {
                current += 1 << (i % 8);
              }
              if ((i & 7) == 7) {
                fingerprint += current;
                current = 0;
              }
            }
            size_t tail = key.size() % 8;
            if (tail > 0) {
              // Just aesthetic reasons to have fingerprint consistently compact form left to right
              //current <<= (8-tail);
              fingerprint += current;
            }

            // Create a key from features
            key = create_key(features);
          }

          //TODO: Short could be enough, but will be measurable difference in the effectivity?
          std::map<int, std::vector<uint32_t> > distances;

          auto group = KNOWLEDGE_BASE.find(key);
          if (group != KNOWLEDGE_BASE.end()) {
            //? As I read somewhere, if key does not exists, it return default value. Is it true? Does it return an empty vector, or some like null - i.e. is .size() valid? And is it effective?
            //? It often does not find the key even if it does exists!?
            std::vector<uint32_t>* templates = &group->second[fingerprint];
            //? 'templates = filter_siblings(id, templates)' create new instance of the output, so it is necessary to expand it into more complex branching?
            //?   Or extract it into function, but it is also paid, right? How about inline functions, should be some marked as it?
            if (templates->size() >= LIMIT && (templates = filter_siblings(id, templates))->size() >= LIMIT) {
              // TODO: A later filtering is useless for this branch.
              //? Actually, how this works, does not it create new instances of inserted value?
              distances.insert({0, *templates});
            } else {
              //? What is the map::local_iterator?
              for (std::unordered_map<std::string, std::vector<uint32_t> >::iterator it = group->second.begin(); it != group->second.end(); it++) {
                int distance = 0;
                // Consideration that fingerprints have the same size.
                for (size_t i = 0; i < fingerprint.length(); i++) {
                  //? Does this retype affect effectivity?
                  distance += DISTANCES[(unsigned char)(fingerprint[i] ^ it->first[i])];
                }
                std::vector<uint32_t>& list = distances[distance];
                // TODO: It probably will be more effective to have vector of vectors instead of merging them?
                list.insert(list.end(), it->second.begin(), it->second.end());
              }
            }
          }

          // Write results of query on the current element to the file.
          // Output strongly reduced in comparison to C# version. Stats (if required) will be computed later.
          stream << id << '\n';
          size_t count = 0;
          for (std::map<int, std::vector<uint32_t> >::iterator it = distances.begin(); count < LIMIT && it != distances.end(); it++) {
            //? Does 'templates = filter_siblings(id, templates)' create new instance of the output, so it is necessary to use a new variable?
            std::vector<uint32_t>* templates = filter_siblings(id, &it->second);
            count += templates->size();
            for (std::vector<uint32_t>::iterator i = templates->begin(); i != templates->end(); i++) {
              uint32_t t = *i;
              // Format: index of the element (to allow stats)  \t  distance
              stream << t << '\t' << it->first << '\n';
            }
            if (templates->size() > 0) {
              std::string info = std::to_string(id) + "\t" + std::to_string(it->first) + "\t" + std::to_string(count) + "\t\r";
            }
          }
          stream << std::endl;
        }
        stream.close();
      }


      public:
      // Load knowledge base in memory
      // <filters> Order of features for prefilter fingerprints according to central residue
      Mine(std::string knowledge_base, const std::set<std::string> filters, char threads, int limit) : FILTERS(filters), THREADS(threads), LIMIT(limit) {
        for (int i = 1; i < 256; i *= 2) {
          for (int j = 0; j < i; j++) {
            DISTANCES[i + j] = DISTANCES[j] + 1;
          }
        }

        // Basic check
        if (!common::filesystem::exists(knowledge_base)) {
          throw common::exception::TitledException("Path to the knowledge '" + knowledge_base + "' base does not exist");
        }
        if (!common::filesystem::is_directory(knowledge_base)) {
          throw common::exception::TitledException("Path to the knowledge base '" + knowledge_base + "' is not a directory");
        }

        // Previously loaded elements to reuse them in the case of recurrence and thus save RAM
        std::unordered_set<std::string> fingerprints;

        // Iterate through knowledgebase files
        common::filesystem::RecursiveDirectoryFileIterator file_iterator(knowledge_base);
        if (file_iterator.has_file()) {
          do {
            if (!common::filesystem::is_regular_file(file_iterator.filename())) continue;
            std::string relative = common::filesystem::relative(file_iterator.filename(), knowledge_base);
            if (kb_path_check(std::string(1, common::filesystem::directory_separator) + relative, filters)) {
              // Processing of knowledgebase file
              relative.erase(relative.length() - 4);

              // Creating a key according to considered features (ordered by feature name to ensure unambiguity)
              std::stringstream name(relative);
              std::string part;
              std::map<std::string, std::string> features;
              while (std::getline(name, part, common::filesystem::directory_separator)) {
                std::string::size_type position = part.find('_');
                if (position != std::string::npos) {
                  std::string key = part.substr(0, position);
                  std::set<std::string>::iterator it = filters.find(key);
                  if (it != filters.end()) {
                    // TODO: Binning of the feature. Maybe JSON file for creation of fingerprints can be reused. Just optional.
                    features.insert({key, part.substr(position + 1)});
                  }
                }
              }
              const std::string& key_feature = create_key(features);

              // Create and find group of fingerprints for the current file
              std::unordered_map<std::string, std::unordered_map<std::string, std::vector<uint32_t> > >::iterator iterator = KNOWLEDGE_BASE.find(key_feature);
              if (iterator == KNOWLEDGE_BASE.end()) {
                KNOWLEDGE_BASE.insert({key_feature, std::unordered_map<std::string, std::vector<uint32_t> >()});
              }
              std::unordered_map<std::string, std::vector<uint32_t> >& group = KNOWLEDGE_BASE.at(key_feature);

              // Read fingerprints from the current file
              std::ifstream stream(file_iterator.filename(), std::ios::in | std::ios::binary);
              if (stream.is_open()) {
                uint32_t counter = 0;
                uint32_t length;
                stream.read(reinterpret_cast<char *>(&length), sizeof(length));
                while (stream.peek() != EOF) {
                  uint32_t id;
                  stream.read(reinterpret_cast<char *>(&id), sizeof(id));
                  //? Must be initialized for that size?
                  std::string key_fingerprint(length, '\0');
                  stream.read(&key_fingerprint[0], length);

                  std::unordered_set<std::string>::iterator fp_iterator = fingerprints.find(key_fingerprint);
                  if (fp_iterator == fingerprints.end()) {
                    fingerprints.insert(key_fingerprint);
                  } else {
                    key_fingerprint = *fp_iterator;
                  }

                  std::unordered_map<std::string, std::vector<uint32_t> >::iterator kb_iterator = group.find(key_fingerprint);
                  if (kb_iterator == group.end()) {
                    //? How to initialize a vector with size 1? (e.g. in Java, the default initial size of List is 10, which is uneffective in the case most of Lists contain single item only and must be restricted explicit.)
                    std::vector<uint32_t> tmp;
                    tmp.push_back(id);
                    group.insert({key_fingerprint, tmp});
                  } else {
                    kb_iterator->second.push_back(id);
                  }
                  counter++;
                }
              }
            }
          } while (file_iterator.has_next());
        } else {
          throw common::exception::TitledException("There is no file in the directory '" + knowledge_base + "'");
        }
      }

      // Load what entries should be skipped for particular queries
      void load_excludes(std::string path) {
        std::ifstream stream(path, std::ios::in);
        if (stream.is_open()) {
          std::string line;
          while (std::getline(stream, line)) {
            size_t tab = line.find('\t');
            if (tab >= 0) {
              std::stringstream values(line.substr(tab+1));
              std::unordered_set<int>* exclude = new std::unordered_set<int>(std::istream_iterator<int>(values), std::istream_iterator<int>());
              std::stringstream keys(line.substr(0, tab));
              int id;
              while (keys >> id) {
                EXCLUDE.insert({id, exclude});
              }
            }
          }
        }
      }

      void clear_excludes() {
        EXCLUDE.clear();
      }

      void threads(char threads) {
        THREADS = threads;
      }

      void limit(int limit) {
        LIMIT = limit;
      }

      void select(std::string input, std::string output) {
        if (output.empty() || output.back() == common::filesystem::directory_separator) {
          size_t i = input.rfind(common::filesystem::directory_separator);
          std::string tmp = (i == input.npos ? input : input.substr(i+1));
          if (common::string::ends_with(tmp, ".fit")) {
            output += tmp.substr(0, tmp.size()-4);
          } else {
            output += tmp;
          }
        } else if (common::string::ends_with(output, ".med")) {
          output = output.substr(0, output.size()-4);
        }
        load_tasks(input);

        std::map<std::string, std::thread> threads;
        for (size_t i = 0; i < THREADS; i++) {
          std::string temp = output + "." + std::to_string(i) + ".med";
          threads[temp] = std::thread(&Mine::selectThread, this, temp);
        }
        for (auto threads_it = threads.begin(); threads_it != threads.end(); ++threads_it) {
          threads_it->second.join();
        }
        if (threads.size() >= 1) {
          common::filesystem::move(threads.begin()->first, output + ".med");
          if (threads.size() > 1) {
            std::ofstream out(output + ".med", std::ios_base::binary | std::ios_base::app);
            for (auto threads_it = ++threads.begin(); threads_it != threads.end(); ++threads_it) {
              std::ifstream in(threads_it->first, std::ios_base::binary);
              out << in.rdbuf();
              in.close();
              common::filesystem::remove_file(threads_it->first);
            }
            out.close();
          }
        }
      }

    };
  }
}