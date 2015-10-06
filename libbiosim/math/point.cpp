#include "math/point.h"
#include "math/floating_point.h"

namespace biosim {
  namespace math {
    // ctor from base class; default ctor
    point::point(std::array<double, 3> __data) : std::array<double, 3>(__data) {}
    // returns true if this point has the same position as the rhs (using almost equal)
    bool point::operator==(point const &__rhs) const {
      return math::almost_equal(operator[](0), __rhs[0]) && math::almost_equal(operator[](1), __rhs[1]) &&
             math::almost_equal(operator[](2), __rhs[2]);
    } // operator==()
  } // namespace math
} // namespace biosim
