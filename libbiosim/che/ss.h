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
      // ctor from pool and optional sequence length, default 0
      explicit ss(std::set<cchb_dssp_interval> __pool, size_t __sequence_length = 0);
      // get secondary structure sequence
      sequence<cchb_dssp> get_sequence() const;
      // get secondary structure intervals
      std::set<cchb_dssp_interval> const &get_sses() const;

    private:
      sequence<cchb_dssp> _sequence; // secondary structure sequence
      bool _sequence_generated; // if the sequence was generated from the sse intervals
      std::set<cchb_dssp_interval> _sses; // secondary structure elements
      size_t _length; // sequence length if created from sses, 0=ends with last sse; otherwise same as sequence length
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
