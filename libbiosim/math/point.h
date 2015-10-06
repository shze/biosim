#ifndef math_point_h
#define math_point_h

#include <array>

namespace biosim {
  namespace math {
    // point is a array of fixed length; implement as class instead of using/typedef to ensure initialization
    class point : public std::array<double, 3> {
    public:
      // ctor from base class; default ctor
      explicit point(std::array<double, 3> __data = {});
      // returns true if this point has the same position as the rhs (using almost equal)
      bool operator==(point const &__rhs) const;
    }; // class point
  } // namespace math
} // namespace biosim

#endif // math_point_h
