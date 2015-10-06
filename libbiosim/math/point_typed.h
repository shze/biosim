#ifndef math_point_typed_h
#define math_point_typed_h

#include "math/point.h"
#include <iostream>

namespace biosim {
  namespace math {
    // extend point with a coordinate type for conversion and output;
    // design: not derived from point, point_typed should never be used in place of point
    class point_typed {
    public:
      // coordinate type
      enum class coordinate_type { cartesian, cylindrical };

      // ctor taking a point and a coordinate type (default cartesian)
      explicit point_typed(point __point, coordinate_type __type = coordinate_type::cartesian);
      // get point data
      point const &get_point() const;
      // get coordinate type
      coordinate_type const &get_type() const;

      // get the data description for the given coordinate_type
      static std::array<std::string, 3> get_data_description(coordinate_type __type);

      // returns a new point converted to given coordinate type
      point to_type(coordinate_type __type) const;

    private:
      // convert cartesian to cylindrical coordinates
      static point cartesian_to_cylindrical(point const &__point);
      // convert cylindrical to cartesian coordinates
      static point cylindrical_to_cartesian(point const &__point);

      point _point; // point data
      coordinate_type _type; // coordinate type
    }; // class point_typed

    // output operator for point_typed
    inline std::ostream &operator<<(std::ostream &__out, point_typed const &__p) {
      std::array<std::string, 3> d(point_typed::get_data_description(__p.get_type()));
      point const &p(__p.get_point());
      __out << d[0] << '=' << p[0] << ", " << d[1] << "=" << p[1] << ", " << d[2] << "=" << p[2];
      return __out;
    } // operator<<()
  } // namespace math
} // namespace biosim

#endif // math_point_typed_h
