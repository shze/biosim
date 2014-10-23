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
      explicit ss(std::vector<cchb_dssp> __sequence);
      // ctor from pool
      explicit ss(std::set<cchb_dssp_interval> __pool);
      // get secondary structure sequence
      std::vector<cchb_dssp> const &get_sequence() const;
      // get secondary structure intervals
      std::set<cchb_dssp_interval> const &get_sses() const;

    private:
      std::vector<cchb_dssp> _sequence; // secondary structure sequence
      std::set<cchb_dssp_interval> _sses; // secondary structure elements
    }; // class ss
  } // namespace che
} // namespace biosim

#endif // che_ss_h
