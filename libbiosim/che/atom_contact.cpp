#include "che/atom_contact.h"

namespace biosim {
  namespace che {
    // ctor from all data members
    atom_contact::atom_contact(assembly_atom __atom1, assembly_atom __atom2, double __distance)
        : std::tuple<assembly_atom, assembly_atom, double>{__atom1, __atom2, __distance} {}
    // get first atom
    assembly_atom const &atom_contact::get_first_atom() const { return std::get<0>(*this); }
    // get second atom
    assembly_atom const &atom_contact::get_second_atom() const { return std::get<1>(*this); }
    // get distance
    double const &atom_contact::get_distance() const { return std::get<2>(*this); }
  } // namespace che
} // namespace biosim