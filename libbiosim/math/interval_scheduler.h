#ifndef math_interval_scheduler_h
#define math_interval_scheduler_h

#include "math/interval.h"
#include <set>

namespace biosim {
  namespace math {
    // interface for interval scheduling algorithms
    template <typename T>
    class interval_scheduler {
    public:
      // method to select a non-overlapping subset of intervals from the input set
      virtual std::set<interval<T>> schedule(std::set<interval<T>> __intervals) = 0;
    }; // class interval_scheduler
  } // namespace math
} // namespace biosim

#endif // math_interval_scheduler_h
