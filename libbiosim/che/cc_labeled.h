#ifndef che_cc_labeled_h
#define che_cc_labeled_h

#include "che/cc.h"

namespace biosim {
  namespace che {
    // cc with data to locate it in an assembly; ensemble sample position is not part b/c samples cannot be averaged.
    // design: derive from cc, allows all code for cc to also process cc_labeled
    class cc_labeled : public cc {
    public:
      // ctor from all data members
      cc_labeled(std::string __chain_id, size_t __cc_id, che::cc __cc);
      // get chain id
      std::string const &get_chain_id() const;
      // get residue id
      size_t const &get_cc_id() const;

    private:
      std::string _chain_id; // chain id
      size_t _cc_id; // residue id
    }; // class cc_labeled
  } // namespace che
} // namespace biosim

#endif // che_cc_labeled_h
