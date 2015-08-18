#ifndef che_bond_h
#define che_bond_h

#include <map>

namespace biosim {
  namespace che {
    // stores a bond
    class bond {
    public:
      // bond types, see http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp_bond.value_order.html
      enum class value_type { arom, delo, doub, pi, poly, quad, sing, trip };

      // ctor from value type string; default ctor
      explicit bond(std::string __type_string = "SING");
      // get bond type
      value_type get_type() const;

      // mapping value type string -> value_type
      static std::map<std::string, value_type> const &get_value_type_map();

    private:
      value_type _type; // bond type
    }; // class bond
  } // namespace che
} // namespace biosim

#endif // che_bond_h
