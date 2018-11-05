#pragma once

#include "../elemental/string.h"

namespace inspire {
  namespace backend {
    // Trivial example, filter nothing
    struct BasicFilter {
      // Tests, whether the line should be parsed (returns TRUE), or skipped (returns FALSE)
      virtual bool keep(const std::string& line) {
        return true;
      }
    };

    // Inverts test of the inner filter
    struct NotFilter : BasicFilter {
      private:
      // Filter that should be inverted
      BasicFilter FILTER;

      public:
      // This constructor accepts the filter whose answers should be inverted
      NotFilter(BasicFilter filter) : FILTER(filter) { }

      bool keep(const std::string& line) override {
        return !FILTER.keep(line);
      }
    };

    // All filters must agree to keep a line
    // NOTE: This filter should be used to conjoin filters with memory, where is necessary to propagate input to all filters.
    struct AndFilter : BasicFilter {
      private:
      // Filters that should be conjoined
      std::list<BasicFilter> FILTERS;

      public:
      // This constructor accepts a list of filters whose answers should be conjoined
      AndFilter(std::list<BasicFilter> filters) : FILTERS(filters) { }

      bool keep(const std::string& line) override {
        bool ret = true;
        for (auto it = FILTERS.begin(); it != FILTERS.end(); ++it) {
          ret &= it->keep(line);
        }
        return ret;
      }
    };

    // All filters must agree to keep a line
    // NOTE: This filter tests until first FALSE is returned, thus it should be used to conjoin flow-filters without memory only.
    struct LazyAndFilter : BasicFilter {
      private:
      // Filters that should be conjoined
      std::list<BasicFilter> FILTERS;

      public:
      // This constructor accepts a list of filters whose answers should be conjoined
      LazyAndFilter(std::list<BasicFilter> filters) : FILTERS(filters) { }

      bool keep(const std::string& line) override {
        for (auto it = FILTERS.begin(); it != FILTERS.end(); ++it) {
          if (!it->keep(line)) {
            return false;
          }
        }
        return true;
      }
    };

    // At least one filter must agree to keep a line
    // NOTE: This filter should be used to disjoin filters with memory, where is necessary to propagate input to all filters.
    struct OrFilter : BasicFilter {
      private:
      // Filters that should be conjoined
      std::list<BasicFilter> FILTERS;

      public:
      // This constructor accepts a list of filters whose answers should be conjoined
      OrFilter(std::list<BasicFilter> filters) : FILTERS(filters) { }

      bool keep(const std::string& line) override {
        bool ret = false;
        for (auto it = FILTERS.begin(); it != FILTERS.end(); ++it) {
          ret |= it->keep(line);
        }
        return ret;
      }
    };

    // At least one filter must agree to keep a line
    // NOTE: This filter tests until first TRUE is returned, thus it should be used to cdisjoin flow-filters without memory only.
    struct LazyOrFilter : BasicFilter {
      private:
      // Filters that should be conjoined
      std::list<BasicFilter> FILTERS;

      public:
      // This constructor accepts a list of filters whose answers should be conjoined
      LazyOrFilter(std::list<BasicFilter> filters) : FILTERS(filters) { }

      bool keep(const std::string& line) override {
        for (auto it = FILTERS.begin(); it != FILTERS.end(); ++it) {
          if (it->keep(line)) {
            return true;
          }
        }
        return false;
      }
    };

    // Skips all water defined in 'HETATM' line;
    // it is expected, that water is marked as 'HOH'.
    // TODO: Skip it outside chains only? It will require memory of the filter and duplicite parsing of a file though in a limited form.
    struct SkipH2OFilter : BasicFilter {
      bool keep(const std::string& line) override {
        return line.size() < 20 || !elemental::string::starts_with(line, "HETATM") || !elemental::string::contains_at(line, "HOH", 17);
      }
    };

    // Skips all atoms in 'ATOM' section except Carbon alphas
    // NOTE: If an aminoacid has no C_alpha specified, it is lost
    struct CalphaOnlyFilter : BasicFilter {
      bool keep(const std::string& line) override {
        return line.size() < 16 || !elemental::string::starts_with(line, "ATOM  ") || elemental::string::contains_at(line, " CA ", 12);
      }
    };

    // Skips all hydrogens in 'ATOM' section
    struct SkipHydrogenFilter : BasicFilter {
      bool keep(const std::string& line) override {
        return line.size() < 78 || !elemental::string::starts_with(line, "ATOM  ") || elemental::string::contains_at(line, " H", 76);
      }
    };

    // Skips header
    // NOTE: Usefull for basic atomic features, where interdistances between chains are not used, e.g. aminoacid type, charge or temperature factor.
    struct SkipRemarksFilter : BasicFilter {
      bool keep(const std::string& line) override {
        return !elemental::string::starts_with(line, "REMARK ");
      }
    };

  }
}