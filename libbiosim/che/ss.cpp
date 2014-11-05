#include "che/ss.h"
#include "math/interval_scheduler_maximize.h"

namespace biosim {
  namespace che {
    // ctor from secondary structure sequence
    ss::ss(sequence<cchb_dssp> __sequence)
        : _sequence(__sequence),
          _sequence_generated(false),
          _sses(cchb_dssp_interval::to_sequence_intervals(__sequence)) {}
    // ctor from pool
    ss::ss(std::set<cchb_dssp_interval> __pool) : _sequence(), _sequence_generated(true), _sses(__pool) {}
    // get secondary structure sequence
    sequence<cchb_dssp> ss::get_sequence() const {
      if(!_sequence_generated) {
        return _sequence;
      }

      math::interval_scheduler_maximize<cchb_dssp_interval> scheduler;
      std::set<cchb_dssp_interval> optimized_sses(scheduler.schedule(_sses));
      return sequence_interval<cchb_dssp>::to_sequence(optimized_sses, 0);
    } // get_sequence()
    // get secondary structure intervals
    std::set<cchb_dssp_interval> const &ss::get_sses() const { return _sses; }
  } // namespace che
} // namespace biosim
