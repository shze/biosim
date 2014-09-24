#ifndef che_cchb_h
#define che_cchb_h

#include <map>
#include <list>
#include <memory>
#include "tools/enumerate.h"

namespace biosim {
  namespace che {
    // exception specific for cchb
    struct cchb_data_not_found : std::runtime_error {
      explicit cchb_data_not_found(std::string const &__s) : std::runtime_error("cchb data not found: " + __s) {}
    };

    // stores the hydrogen bond state for a cc; design see cc
    class cchb {
    public:
      // specificity defined as enum for easier comparison; only one specificity allowed, no combinations.
      enum specificity_type { specificity_defined, specificity_unknown, specificity_gap };
      using weight_map = std::map<std::string, double>; // simplify naming

      // ctor from id; id is assumed unique; throws if not found
      explicit cchb(std::string __id);
      // ctor from specificity and monomer_type; assumed not unique, constructs first
      explicit cchb(specificity_type __specificity);
      // ctor from id and weight set
      cchb(std::string __id, weight_map __weights);

      // get identifier
      std::string get_identifier() const;
      // get specificity
      specificity_type get_specificity() const;
      // return if an unknown (unspecified) compound; only for unknown, not for gap.
      bool is_unknown() const;
      // return if gap
      bool is_gap() const;
      // get profile weights
      weight_map get_weights() const;

      // define ordering
      bool operator<(cchb const &__rhs) const;
      // equality
      bool operator==(cchb const &__rhs) const;

      // get list of all identifiers; static public interface
      static std::list<std::string> get_id_list();

    private:
      // contains the cchb data; type that is enumerated
      struct data {
        // ctor for primary symbol
        data(std::string __id, specificity_type __specificity = specificity_defined);
        // ctor for profile symbol
        data(std::string __id, weight_map __weights);

        std::string _id; // identifier
        specificity_type _specificity; // specificity
        weight_map _weights; // weights of this data object based on compounds in the enumerated instance set
      }; // struct data

      // define identifier function for use with enumerate; friend and not part of the class
      friend std::string identifier(std::shared_ptr<data const> const &__data_ptr) { return __data_ptr->_id; }

      using data_enum = tools::enumerate<std::shared_ptr<data const>>; // simplify naming
      // find compound with id; id is assumed unique
      static std::shared_ptr<data const> find(std::string __id);
      // find compound with specificity; NOT assumed unique, only the first compound is returned
      static std::shared_ptr<data const> find(specificity_type __specificity);
      static bool _initialized; // static variable to initialize data_enum with a minimal set of data objects
      // constructs data objects of the enumerated instance set; subset of states from DSSP is used, see:
      // http://en.wikipedia.org/wiki/DSSP_%28protein%29
      static bool initialize();

      std::shared_ptr<data const> _impl; // ptr to implementation; shared_ptr is acceptable, b/c it points to const data
    }; // class cc
  } // namespace che
} // namespace biosim

#endif // che_cchb_h
