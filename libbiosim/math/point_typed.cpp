#include "math/point_typed.h"
#include <cmath>

namespace biosim {
  namespace math {
    // ctor taking a point and a coordinate type (default cartesian)
    point_typed::point_typed(point __point, coordinate_type __type) : _point(__point), _type(__type) {}
    // get point data
    point const &point_typed::get_point() const { return _point; }
    // get coordinate type
    point_typed::coordinate_type const &point_typed::get_type() const { return _type; }

    // get the data description for the given coordinate_type
    std::array<std::string, 3> point_typed::get_data_description(coordinate_type __type) {
      switch(__type) {
      case coordinate_type::cartesian:
        return {"x", "y", "z"};
      case coordinate_type::cylindrical:
        return {"rho", "phi", "z"};
      } // switch
    } // get_data_description()

    // returns a new point converted to given coordinate type
    point point_typed::to_type(coordinate_type __type) const {
      if(get_type() == __type) { // if no conversion if needed, just return the data
        return get_point();
      } // if

      point cartesian_point;
      switch(_type) { // first, convert all other coordinate types into cartesian coordinates
      case coordinate_type::cartesian:
        cartesian_point = _point;
        break;
      case coordinate_type::cylindrical:
        cartesian_point = cylindrical_to_cartesian(_point);
        break;
      } // switch
      switch(__type) { // then convert the cartesian coordinates into the requested coordinate type and return it
      case coordinate_type::cartesian:
        return cartesian_point;
      case coordinate_type::cylindrical:
        return cartesian_to_cylindrical(cartesian_point);
      } // switch
    } // to_type()

    // convert cartesian to cylindrical coordinates
    point point_typed::cartesian_to_cylindrical(point const &__point) {
      double cartesian_x(__point[0]), cartesian_y(__point[1]), cartesian_z(__point[2]);
      double cylindrical_rho = sqrt(cartesian_x * cartesian_x + cartesian_y * cartesian_y);
      double cylindrical_phi = atan2(cartesian_y, cartesian_x);
      double cylindrical_z = cartesian_z;
      return point({cylindrical_rho, cylindrical_phi, cylindrical_z});
    } // cartesian_to_cylindrical()
    // convert cylindrical to cartesian coordinates
    point point_typed::cylindrical_to_cartesian(point const &__point) {
      double cylindrical_rho(__point[0]), cylindrical_phi(__point[1]), cylindrical_z(__point[2]);
      double cartesian_x = cylindrical_rho * cos(cylindrical_phi);
      double cartesian_y = cylindrical_rho * sin(cylindrical_phi);
      double cartesian_z = cylindrical_z;
      return point({cartesian_x, cartesian_y, cartesian_z});
    } // cylindrical_to_cartesian()
  } // namespace math
} // namespace biosim
