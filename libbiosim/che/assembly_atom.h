#ifndef che_assembly_atom_h
#define che_assembly_atom_h

#include "che/atom.h"

namespace biosim {
  namespace che {
    // atom with data to locate it in an assembly; ensemble sample position is not part b/c samples cannot be averaged.
    // design: derive from atom, allows all code for atom to also process atom_labeled
    class assembly_atom : public atom {
    public:
      // ctor from all data members
      assembly_atom(std::string __chain_id, size_t __cc_id, che::atom __atom);
      // get chain id
      std::string const &get_chain_id() const;
      // get residue_id
      size_t const &get_cc_id() const;
      // returns true if this identifier is less than the rhs identifier
      bool operator<(assembly_atom const &__rhs) const;

    private:
      std::string _chain_id; // chain id
      size_t _cc_id; // residue id
    }; // class assembly_atom
  } // namespace che
} // namespace biosim

#endif // che_assembly_atom_h
