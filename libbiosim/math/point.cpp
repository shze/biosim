#include "math/point.h"
#include <cmath>

namespace biosim {
  namespace math {
    // ctor taking array of three values and a coordinate type (defaults to cartesian coordinates)
    point::point(std::array<double, 3> __data, coordinate_type __type) : _data(__data), _type(__type) {}
    // get coordinate type
    point::coordinate_type point::get_coordinate_type() const { return _type; }
    // get first value: cartesian x, cylindrical rho
    double point::first() const { return _data[0]; }
    // get second value: cartesian y, cylindrical phi
    double point::second() const { return _data[1]; }
    // get third value: cartesian z, cylindrical z
    double point::third() const { return _data[2]; }

    // returns a new point converted to given coordinate type
    point point::to_coordinate_type(coordinate_type __type) const {
      if(get_coordinate_type() == __type) { // if no conversion if needed, just return the data
        return point(*this);
      } // if

      std::array<double, 3> cartesian_data;
      switch(_type) { // first, convert all other coordinate types into cartesian coordinates
      case coordinate_type::cartesian:
        cartesian_data = _data;
      case coordinate_type::cylindrical:
        cartesian_data = cylindrical_to_cartesian(_data);
      } // switch
      switch(__type) { // then convert the cartesian coordinates into the requested coordinate type and return it
      case coordinate_type::cartesian:
        return point(cartesian_data, __type);
      case coordinate_type::cylindrical:
        return point(cartesian_to_cylindrical(cartesian_data), __type);
      } // switch
    } // to_coordinate_type()
    // (static) get the data description for the given coordinate_type
    std::array<std::string, 3> point::get_data_description(coordinate_type __type) {
      switch(__type) {
      case coordinate_type::cartesian:
        return {"x", "y", "z"};
      case coordinate_type::cylindrical:
        return {"rho", "phi", "z"};
      } // switch
    } // get_data_description()

    // (static) convert cartesian to cylindrical coordinates
    std::array<double, 3> point::cartesian_to_cylindrical(std::array<double, 3> __data) {
      double cartesian_x(__data[0]), cartesian_y(__data[1]), cartesian_z(__data[2]);
      double cylindrical_rho = sqrt(cartesian_x * cartesian_x + cartesian_y * cartesian_y);
      double cylindrical_phi = atan(cartesian_y / cartesian_x);
      double cylindrical_z = cartesian_z;
      return {cylindrical_rho, cylindrical_phi, cylindrical_z};
    } // cartesian_to_cylindrical()
    // (static) convert cylindrical to cartesian coordinates
    std::array<double, 3> point::cylindrical_to_cartesian(std::array<double, 3> __data) {
      double cylindrical_rho(__data[0]), cylindrical_phi(__data[1]), cylindrical_z(__data[2]);
      double cartesian_x = cylindrical_rho * cos(cylindrical_phi);
      double cartesian_y = cylindrical_rho * sin(cylindrical_phi);
      double cartesian_z = cylindrical_z;
      return {cartesian_x, cartesian_y, cartesian_z};
    } // cylindrical_to_cartesian()
  } // namespace math
} // namespace biosim
