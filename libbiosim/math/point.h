#ifndef math_point_h
#define math_point_h

#include <array>
#include <ostream>
#include "math/floating_point.h"

namespace biosim {
  namespace math {
    // stores a three-dimensional point in a coordinate system type
    class point {
    public:
      // coordinate type
      enum class coordinate_type { cartesian, cylindrical };

      // ctor taking array of three values (default zeros) and a coordinate type (default cartesian); default ctor
      explicit point(std::array<double, 3> __data = {0, 0, 0}, coordinate_type __type = coordinate_type::cartesian);
      // get coordinate type
      coordinate_type get_coordinate_type() const;
      // get first value: cartesian x, cylindrical rho
      double first() const;
      // get second value: cartesian y, cylindrical phi
      double second() const;
      // get third value: cartesian z, cylindrical z
      double third() const;
      // returns true if this point has the same position as the rhs
      bool operator==(point const &__rhs) const;

      // returns a new point converted to given coordinate type
      point to_coordinate_type(coordinate_type __type) const;
      // get the data description for the given coordinate_type
      static std::array<std::string, 3> get_data_description(coordinate_type __type);

    private:
      // convert cartesian to cylindrical coordinates
      static std::array<double, 3> cartesian_to_cylindrical(std::array<double, 3> __data);
      // convert cylindrical to cartesian coordinates
      static std::array<double, 3> cylindrical_to_cartesian(std::array<double, 3> __data);

      std::array<double, 3> _data; // three data values
      coordinate_type _type; // coordinate type
    }; // class point

    // returns if both points have the same type and position
    inline bool equal_position_and_type(point const &__p1, point const &__p2) {
      return __p1.get_coordinate_type() == __p2.get_coordinate_type() &&
             math::almost_equal(__p1.first(), __p2.first()) && math::almost_equal(__p1.second(), __p2.second()) &&
             math::almost_equal(__p1.third(), __p2.third());
    } // equal_position_and_type()
    // returns if both points have the position independ of their coordinate system type
    inline bool equal_position(point const &__p1, point const &__p2) {
      return equal_position_and_type(__p1, __p2.to_coordinate_type(__p1.get_coordinate_type()));
    } // equal_position()

    // output operator for point
    inline std::ostream &operator<<(std::ostream &__out, point const &__p) {
      std::array<std::string, 3> d(point::get_data_description(__p.get_coordinate_type()));
      __out << d[0] << "=" << __p.first() << ", " << d[1] << "=" << __p.second() << ", " << d[2] << "=" << __p.third();
      return __out;
    } // operator<<()
  } // namespace math
} // namespace biosim

#endif // math_point_h
