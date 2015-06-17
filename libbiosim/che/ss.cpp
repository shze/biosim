#include "che/ss.h"
#include "math/algo/interval_scheduler_maximize.h"

namespace biosim {
  namespace che {
    // default ctor
    ss::ss() : _input(input_type::empty), _sequence(), _sses(), _length(0) {}
    // ctor from secondary structure sequence
    ss::ss(sequence<cchb_dssp> __sequence)
        : _input(input_type::sequence),
          _sequence(__sequence),
          _sses(cchb_dssp_interval::to_sequence_intervals(__sequence)),
          _length(__sequence.size()) {}
    // ctor from pool
    ss::ss(std::set<cchb_dssp_interval> __pool, size_t __sequence_length)
        : _input(input_type::pool), _sequence(), _sses(__pool), _length(__sequence_length) {
      // convert sses into sequence
      math::algo::interval_scheduler_maximize<cchb_dssp_interval> scheduler;
      std::set<cchb_dssp_interval> optimized_sses(scheduler.schedule(_sses));
      _sequence = sequence_interval<cchb_dssp>::to_sequence(optimized_sses, _length);
      _length = std::max(_length, _sequence.size());
    }
    // return if defined
    bool ss::defined() const { return _input != input_type::empty; }
    // get secondary structure sequence
    sequence<cchb_dssp> const &ss::get_sequence() const { return _sequence; }
    // get secondary structure intervals
    std::set<cchb_dssp_interval> const &ss::get_sses() const { return _sses; }
    // return the length
    size_t ss::get_length() const { return _length; }
  } // namespace che
} // namespace biosim
