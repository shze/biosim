#ifndef che_bond_h
#define che_bond_h

namespace biosim {
  namespace che {
    // stores a bond
    class bond {
    public:
      // bond types, see http://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_chem_comp_bond.value_order.html
      enum class value_type { arom, delo, doub, pi, poly, quad, sing, trip };

      // ctor from value_type; default ctor
      explicit bond(value_type __type = value_type::sing);

    private:
      value_type _type; // bond type
    }; // class bond
  } // namespace che
} // namespace biosim

#endif // che_bond_h
