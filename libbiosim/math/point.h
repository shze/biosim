#ifndef math_point_h
#define math_point_h

#include <array>
#include <ostream>

namespace biosim {
  namespace math {
    // stores a three-dimensional point in a coordinate system type
    class point {
    public:
      // coordinate type
      enum class coordinate_type { cartesian, cylindrical };

      // ctor taking array of three values and a coordinate type (defaults to cartesian coordinates)
      point(std::array<double, 3> __data, coordinate_type __type = coordinate_type::cartesian);
      // get coordinate type
      coordinate_type get_coordinate_type() const;
      // get first value: cartesian x, cylindrical rho
      double first() const;
      // get second value: cartesian y, cylindrical phi
      double second() const;
      // get third value: cartesian z, cylindrical z
      double third() const;

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

    // output operator for point
    inline std::ostream &operator<<(std::ostream &__out, point const &__p) {
      std::array<std::string, 3> d(point::get_data_description(__p.get_coordinate_type()));
      __out << d[0] << "=" << __p.first() << ", " << d[1] << "=" << __p.second() << ", " << d[2] << "=" << __p.third();
      return __out;
    } // operator<<()
  } // namespace math
} // namespace biosim

#endif // math_point_h
