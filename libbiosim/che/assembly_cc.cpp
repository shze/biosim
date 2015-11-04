#include "che/assembly_cc.h"

namespace biosim {
  namespace che {
    // ctor from all data members
    assembly_cc::assembly_cc(std::string __chain_id, size_t __chain_pos, che::cc __cc)
        : cc(__cc), _chain_id(__chain_id), _chain_pos(__chain_pos) {}
    // get chain id
    std::string const &assembly_cc::get_chain_id() const { return _chain_id; }
    // get residue id
    size_t const &assembly_cc::get_chain_pos() const { return _chain_pos; }
  } // namespace che
} // namespace biosim