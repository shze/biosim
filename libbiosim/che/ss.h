#ifndef che_ss_h
#define che_ss_h

#include "che/sequence_interval.h"

namespace biosim {
  namespace che {
    // secondary structure
    class ss {
    private:
      std::vector<cchb> _sequence; // secondary structure sequence
      std::set<cchb_interval> _sses; // secondary structure elements
    };
  } // namespace che
} // namespace biosim

#endif // che_ss_h
