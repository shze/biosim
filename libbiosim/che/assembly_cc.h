#ifndef che_assembly_cc_h
#define che_assembly_cc_h

#include "che/cc.h"

namespace biosim {
  namespace che {
    // cc with data to locate it in an assembly; ensemble sample position is not part b/c samples cannot be averaged.
    // design: derive from cc, allows all code for cc to also process cc_labeled
    class assembly_cc : public cc {
    public:
      // ctor from all data members
      assembly_cc(std::string __chain_id, size_t __chain_pos, che::cc __cc);
      // get chain id
      std::string const &get_chain_id() const;
      // get residue id
      size_t const &get_chain_pos() const;

    private:
      std::string _chain_id; // chain id
      size_t _chain_pos; // position of cc in the chain, i.e. residue id
    }; // class assembly_cc
  } // namespace che
} // namespace biosim

#endif // che_assembly_cc_h
