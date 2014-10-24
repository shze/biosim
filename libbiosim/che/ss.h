#ifndef che_ss_h
#define che_ss_h

#include "che/sequence.h"
#include "che/sequence_interval.h"
#include <vector>

namespace biosim {
  namespace che {
    // secondary structure
    class ss {
    public:
      // ctor from secondary structure sequence
      explicit ss(sequence<cchb_dssp> __sequence);
      // ctor from pool
      explicit ss(std::set<cchb_dssp_interval> __pool);
      // get secondary structure sequence
      sequence<cchb_dssp> const &get_sequence() const;
      // get secondary structure intervals
      std::set<cchb_dssp_interval> const &get_sses() const;

    private:
      sequence<cchb_dssp> _sequence; // secondary structure sequence
      std::set<cchb_dssp_interval> _sses; // secondary structure elements
    }; // class ss

    // output operator for ss
    inline std::ostream &operator<<(std::ostream &__out, ss const &__ss) {
      __out << __ss.get_sequence();
      for(auto const &sse : __ss.get_sses()) {
        __out << "\n" << sse;
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_ss_h
