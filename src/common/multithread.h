#pragma once

#include "boost/thread/sync_queue.hpp"

namespace common {
  namespace multithread {
    typedef boost::sync_queue<std::string> synchronized_queue;
  }
}