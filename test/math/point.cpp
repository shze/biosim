#include <boost/test/unit_test.hpp>

#include "math/point.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_point)

BOOST_AUTO_TEST_CASE(point_ctor) {
  math::point p;
  BOOST_CHECK(p.get_coordinate_type() == math::point::coordinate_type::cartesian);
  BOOST_CHECK(p.first() == 0.0);
  BOOST_CHECK(p.second() == 0.0);
  BOOST_CHECK(p.third() == 0.0);

  p = math::point({1.0, 2.0, 3.0});
  BOOST_CHECK(p.get_coordinate_type() == math::point::coordinate_type::cartesian);
  BOOST_CHECK(p.first() == 1.0);
  BOOST_CHECK(p.second() == 2.0);
  BOOST_CHECK(p.third() == 3.0);

  p = math::point({1.0, 2.0, 3.0}, math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p.first() == 1.0);
  BOOST_CHECK(p.second() == 2.0);
  BOOST_CHECK(p.third() == 3.0);
}

BOOST_AUTO_TEST_CASE(point_cmp_operator) {
  math::point p_cart({3.0, 4.0, 5.0});
  math::point p_cart2({3.0, 4.0, 6.0});
  math::point p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  math::point p_cyl2 = p_cart2.to_coordinate_type(math::point::coordinate_type::cylindrical);

  BOOST_CHECK(math::equal_position_and_type(p_cart, p_cart) == true);
  BOOST_CHECK(math::equal_position_and_type(p_cart, p_cart2) == false);
  BOOST_CHECK(math::equal_position_and_type(p_cart, p_cyl) == false);
  BOOST_CHECK(math::equal_position_and_type(p_cart, p_cyl2) == false);

  BOOST_CHECK(math::equal_position(p_cart, p_cart) == true);
  BOOST_CHECK(math::equal_position(p_cart, p_cart2) == false);
  BOOST_CHECK(math::equal_position(p_cart, p_cyl) == true);
  BOOST_CHECK(math::equal_position(p_cart, p_cyl2) == false);

  BOOST_CHECK(p_cart == p_cart);
  BOOST_CHECK(!(p_cart == p_cart2));
  BOOST_CHECK(p_cart == p_cyl);
  BOOST_CHECK(!(p_cart == p_cyl2));
}

BOOST_AUTO_TEST_CASE(point_convert) {
  math::point p_cart({3.0, 4.0, 5.0});
  math::point p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 5.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), 0.927295, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  math::point p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK(p_cart2.get_coordinate_type() == math::point::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2.first(), 3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2.second(), 4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({-3.0, 4.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 5.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), 2.214297, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK(p_cart2.get_coordinate_type() == math::point::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2.first(), -3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2.second(), 4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({-3.0, -4.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 5.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), -2.214297, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK(p_cart2.get_coordinate_type() == math::point::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2.first(), -3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2.second(), -4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({3.0, -4.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 5.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), -0.927295, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2.first(), 3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2.second(), -4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({0.0, 4.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 4.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), 1.5708, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2.first(), 1e-5);
  BOOST_CHECK_CLOSE(p_cart2.second(), 4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({0.0, -4.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 4.0);
  BOOST_CHECK_CLOSE(p_cyl.second(), -1.5708, 1e-3);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2.first(), 1e-5);
  BOOST_CHECK_CLOSE(p_cart2.second(), -4.0, 1e-3);
  BOOST_CHECK(p_cart2.third() == 5.0);

  p_cart = math::point({0.0, 0.0, 5.0});
  p_cyl = p_cart.to_coordinate_type(math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.get_coordinate_type() == math::point::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl.first() == 0.0);
  BOOST_CHECK(p_cyl.second() == 0.0);
  BOOST_CHECK(p_cyl.third() == 5.0);
  p_cart2 = p_cyl.to_coordinate_type(math::point::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2.first(), 1e-5);
  BOOST_CHECK_SMALL(p_cart2.second(), 1e-5);
  BOOST_CHECK(p_cart2.third() == 5.0);
}

BOOST_AUTO_TEST_SUITE_END()
