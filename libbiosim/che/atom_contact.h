#ifndef che_atom_contact_h
#define che_atom_contact_h

#include <tuple>
#include <iostream>
#include "che/assembly_atom.h"

namespace biosim {
  namespace che {
    // stores an atom contact, i.e. the two atoms with information about their assembly, and the distance between them
    // design: derived from tuple to leverage tuple e.g. operator<
    struct atom_contact : public std::tuple<assembly_atom, assembly_atom, double> {
      // ctor from all data members
      atom_contact(assembly_atom __atom1, assembly_atom __atom2, double __distance);
      // get first atom
      assembly_atom const &get_first_atom() const;
      // get second atom
      assembly_atom const &get_second_atom() const;
      // get distance
      double const &get_distance() const;
    }; // struct atom_contact

    // output operator for atom_contact
    inline std::ostream &operator<<(std::ostream &__out, atom_contact const &__atom_contact) {
      assembly_atom const &atom1(__atom_contact.get_first_atom());
      assembly_atom const &atom2(__atom_contact.get_second_atom());
      __out << "atom_pair=(" << atom1 << ", " << atom2 << "), euclidean_distance=" << __atom_contact.get_distance();
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_atom_contact_h
