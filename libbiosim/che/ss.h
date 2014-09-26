#ifndef che_ss_h
#define che_ss_h

#include "che/sequence_interval.h"
#include <vector>

namespace biosim {
  namespace che {
    // secondary structure
    class ss {
    public:
      // ctor from secondary structure sequence
      explicit ss(std::vector<cchb> __sequence);
      // ctor from pool
      explicit ss(std::set<cchb_interval> __pool);
      // get secondary structure sequence
      std::vector<cchb> const &get_sequence() const;
      // get secondary structure intervals
      std::set<cchb_interval> const &get_sses() const;

    private:
      std::vector<cchb> _sequence; // secondary structure sequence
      std::set<cchb_interval> _sses; // secondary structure elements
    };
  } // namespace che
} // namespace biosim

#endif // che_ss_h
