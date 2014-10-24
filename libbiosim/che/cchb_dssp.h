#ifndef che_cchb_dssp_h
#define che_cchb_dssp_h

#include <map>
#include <list>
#include <memory>
#include "tools/enumerate.h"

namespace biosim {
  namespace che {
    // exception specific for cchb_dssp
    struct cchb_dssp_data_not_found : std::runtime_error {
      explicit cchb_dssp_data_not_found(std::string const &__s)
          : std::runtime_error("cchb_dssp data not found: " + __s) {}
    }; // struct cchb_dssp_data_not_found

    // stores the hydrogen bond state for a cc using DSSP symbols; design see cc
    class cchb_dssp {
    public:
      // specificity defined as enum for easier comparison; only one specificity allowed, no combinations.
      enum specificity_type { specificity_defined, specificity_unknown, specificity_gap };
      using weight_map = std::map<char, double>; // simplify naming

      // ctor from id; id is assumed unique; throws if not found
      explicit cchb_dssp(char __id);
      // ctor from specificity; assumed not unique, constructs first
      explicit cchb_dssp(specificity_type __specificity);
      // ctor from id and weight set
      cchb_dssp(char __id, weight_map __weights);

      // get identifier
      char get_identifier() const;
      // get identifier char
      char get_identifier_char() const;
      // get specificity
      specificity_type get_specificity() const;
      // return if an unknown (unspecified) compound; only for unknown, not for gap.
      bool is_unknown() const;
      // return if gap
      bool is_gap() const;
      // get profile weights
      weight_map get_weights() const;

      // define ordering
      bool operator<(cchb_dssp const &__rhs) const;
      // equality
      bool operator==(cchb_dssp const &__rhs) const;

      // get list of all identifiers; static public interface
      static std::list<char> get_id_list();
      // create a string of identifier chars; static public interface
      static std::string get_identifier_char_string();

    private:
      // contains the cchb_dssp data; type that is enumerated
      struct data {
        // ctor for primary symbol
        data(char __id, specificity_type __specificity = specificity_defined);
        // ctor for profile symbol
        data(char __id, weight_map __weights);

        char _id; // identifier
        specificity_type _specificity; // specificity
        weight_map _weights; // weights of this data object based on compounds in the enumerated instance set
      }; // struct data

      // define identifier function for use with enumerate; friend and not part of the class
      friend std::string identifier(std::shared_ptr<data const> const &__data_ptr) {
        return std::string(1, __data_ptr->_id);
      } // identifier()

      using data_enum = tools::enumerate<std::shared_ptr<data const>>; // simplify naming
      // find compound with id; id is assumed unique
      static std::shared_ptr<data const> find(char __id);
      // find compound with specificity; NOT assumed unique, only the first compound is returned
      static std::shared_ptr<data const> find(specificity_type __specificity);
      static bool _initialized; // static variable to initialize data_enum with a minimal set of data objects
      // constructs data objects of the enumerated instance set; subset of states from DSSP is used, see:
      // http://en.wikipedia.org/wiki/DSSP_%28protein%29
      static bool initialize();

      std::shared_ptr<data const> _impl; // ptr to implementation; shared_ptr is acceptable, b/c it points to const data
    }; // class cchb_dssp
  } // namespace che
} // namespace biosim

#endif // che_cchb_dssp_h
