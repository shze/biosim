#include "che/bond.h"

namespace biosim {
  namespace che {
    // ctor from value type string; default ctor
    bond::bond(std::string __type_string) : _type(get_value_type_map().at(__type_string)) {}
    // get bond type
    bond::value_type bond::get_type() const { return _type; }

    // (static) mapping value type string -> value_type
    std::map<std::string, bond::value_type> const &bond::get_value_type_map() {
      static std::map<std::string, value_type> mapping = {{"AROM", value_type::arom},
                                                          {"DELO", value_type::delo},
                                                          {"DOUB", value_type::doub},
                                                          {"PI", value_type::pi},
                                                          {"POLY", value_type::poly},
                                                          {"QUAD", value_type::quad},
                                                          {"SING", value_type::sing},
                                                          {"TRIP", value_type::trip}};
      return mapping;
    } // get_value_type_map()
  } // namespace che
} // namespace biosim