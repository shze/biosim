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
      // default ctor
      ss();
      // ctor from secondary structure sequence
      explicit ss(sequence<cchb_dssp> __sequence);
      // ctor from pool and optional sequence length, default 0
      explicit ss(std::set<cchb_dssp_interval> __pool, size_t __sequence_length = 0);
      // return if defined
      bool defined() const;
      // if the ss is empty, i.e. no sequence and no sses
      sequence<cchb_dssp> const &get_sequence() const;
      // get secondary structure intervals
      std::set<cchb_dssp_interval> const &get_sses() const;
      // return the length
      size_t get_length() const;

    private:
      enum class input_type { empty, sequence, pool }; // input defined as enum; only one value allowed, no combinations
      input_type _input; // the input type used to create this object
      sequence<cchb_dssp> _sequence; // secondary structure sequence
      std::set<cchb_dssp_interval> _sses; // secondary structure elements
      size_t _length; // sequence length if created from sses, 0=ends with last sse; otherwise same as sequence length
    }; // class ss

    // output operator for ss
    inline std::ostream &operator<<(std::ostream &__out, ss const &__ss) {
      __out << __ss.get_sequence() << "\n";
      for(auto const &sse : __ss.get_sses()) {
        __out << sse << "\n";
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_ss_h
