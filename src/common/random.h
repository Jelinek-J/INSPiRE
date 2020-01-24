#pragma once

#include <limits>
#include "../common/exception.h"

namespace common {
  namespace random {
    void set_seed(unsigned int seed) {
      srand(seed);
    }

    // Generates a random number in interval [0; count)
    size_t random_size_t(size_t count) {
      if (count == 0) {
        throw common::exception::TitledException("Interval [0; 0) is empty");
      }
      if (count-1 > RAND_MAX) {
        throw common::exception::TitledException("Random generator for so high numbers (" + std::to_string(count) + ") is not implemented yet");
      }
      size_t max;
      if (count-1 == RAND_MAX) {
        max = count-1;
      } else {
        size_t modulo = RAND_MAX % (count);
        max = modulo == count-1 ? RAND_MAX : RAND_MAX-modulo-1;
      }
      int random;
      do {
        random = rand();
      } while (random > max);
      return random % count;
    }
  }
}