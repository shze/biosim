#include "che/assembly_atom.h"

namespace biosim {
  namespace che {
    // ctor from all data members
    assembly_atom::assembly_atom(std::string __chain_id, size_t __cc_id, che::atom __atom)
        : atom(__atom), _chain_id(__chain_id), _cc_id(__cc_id) {}
    // get chain id
    std::string const &assembly_atom::get_chain_id() const { return _chain_id; }
    // get residue_id
    size_t const &assembly_atom::get_cc_id() const { return _cc_id; }
    // returns true if this identifier is less than the rhs identifier
    bool assembly_atom::operator<(assembly_atom const &__rhs) const {
      return _chain_id < __rhs._chain_id || (!(_chain_id < __rhs._chain_id) && _cc_id < __rhs._cc_id) ||
             (!(_chain_id < __rhs._chain_id) && !(_cc_id < __rhs._cc_id) && atom::operator<(__rhs));
    } // operator<()
  } // namespace che
} // namespace biosim