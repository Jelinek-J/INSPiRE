#pragma once

#include <string>

namespace elemental {
  namespace exception {
    // A general exception with a message that should be used to describe the exception
    class TitledException : std::exception {
      private:
      // explanatory message
      std::string TITLE;

      public:
      // 'title' will be returned in the what
      TitledException(const std::string& title) : TITLE(title) { }

      virtual const char* what() const throw() {
        return TITLE.c_str();
      }
    };
  }
}
