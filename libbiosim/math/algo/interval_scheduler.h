#ifndef math_interval_scheduler_h
#define math_interval_scheduler_h

#include "math/interval.h"
#include <set>

namespace biosim {
  namespace math {
    namespace algo {
      // interface for interval scheduling algorithms
      template <class T>
      class interval_scheduler {
      public:
        // method to select a non-overlapping subset of intervals from the input set
        virtual std::set<T> schedule(std::set<T> __intervals) = 0;
      }; // class interval_scheduler
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_interval_scheduler_h
