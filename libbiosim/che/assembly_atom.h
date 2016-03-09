#ifndef che_assembly_atom_h
#define che_assembly_atom_h

#include <iostream>
#include "che/atom.h"

namespace biosim {
  namespace che {
    // atom with data to locate it in an assembly; ensemble sample position is not part b/c samples cannot be averaged.
    // design: derive from atom, allows all code for atom to also process atom_labeled
    class assembly_atom : public atom {
    public:
      // ctor from all data members
      assembly_atom(std::string __chain_id, size_t __cc_pos, che::atom __atom);
      // get chain id
      std::string const &get_chain_id() const;
      // get residue id
      size_t const &get_chain_pos() const;
      // returns true if this identifier is less than the rhs identifier
      bool operator<(assembly_atom const &__rhs) const;

    private:
      std::string _chain_id; // chain id
      size_t _chain_pos; // position of cc in the chain, i.e. residue id
    }; // class assembly_atom

    // output operator for assembly_atom
    inline std::ostream &operator<<(std::ostream &__out, assembly_atom const &__a) {
      __out << __a.get_chain_id() << __a.get_chain_pos() << __a.get_identifier();
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_assembly_atom_h
