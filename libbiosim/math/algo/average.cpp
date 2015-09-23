#include "math/algo/average.h"

#include <limits>

namespace biosim {
  namespace math {
    namespace algo {
      // default ctor
      cma::cma() : _average(0.0), _n(0) {}
      // observe a given value
      void cma::observe(double __value) { _average += (__value - _average) / ++_n; }
      // get current number of observed values
      size_t cma::get_n() const { return _n; }
      // get current average
      double cma::get_average() const { return _n > 0 ? _average : std::numeric_limits<double>::quiet_NaN(); }
    } // namespace algo
  } // namespace math
} // namespace biosim