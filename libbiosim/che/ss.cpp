#include "che/ss.h"

namespace biosim {
  namespace che {
    // ctor from secondary structure sequence
    ss::ss(std::vector<cchb> __sequence) : _sequence(__sequence), _sses() {}
    // ctor from pool
    ss::ss(std::set<cchb_interval> __pool) : _sequence(), _sses(__pool) {}
    // get secondary structure sequence
    std::vector<cchb> const &ss::get_sequence() const { return _sequence; }
    // get secondary structure intervals
    std::set<cchb_interval> const &ss::get_sses() const { return _sses; }
  } // namespace che
} // namespace biosim