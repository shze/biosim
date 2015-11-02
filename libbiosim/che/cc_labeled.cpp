#include "che/cc_labeled.h"

namespace biosim {
  namespace che {
    // ctor from all data members
    cc_labeled::cc_labeled(std::string __chain_id, size_t __cc_id, che::cc __cc)
        : cc(__cc), _chain_id(__chain_id), _cc_id(__cc_id) {}
    // get chain id
    std::string const &cc_labeled::get_chain_id() const { return _chain_id; }
    // get residue id
    size_t const &cc_labeled::get_cc_id() const { return _cc_id; }
  } // namespace che
} // namespace biosim