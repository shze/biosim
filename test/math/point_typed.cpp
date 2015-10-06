#include <boost/test/unit_test.hpp>

#include "math/point_typed.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_point_typed)

BOOST_AUTO_TEST_CASE(point_typed_ctor) {
  math::point p_data({1.0, 2.0, 3.0});
  math::point_typed p(p_data);
  BOOST_CHECK(p.get_type() == math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK(p.get_point() == p_data);

  p = math::point_typed(p_data, math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p.get_type() == math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p.get_point() == p_data);
}

BOOST_AUTO_TEST_CASE(point_typed_cmp_operator) {
  math::point p_cart({3.0, 4.0, 5.0});
  math::point p_cart2({3.0, 4.0, 6.0});
  math::point p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  math::point p_cart_again = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                                 .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK(p_cart == p_cart);
  BOOST_CHECK(p_cart == p_cart_again);
  BOOST_CHECK(p_cart != p_cart2);
  BOOST_CHECK(p_cart != p_cyl);
}

BOOST_AUTO_TEST_CASE(point_typed_convert) {
  math::point p_cart({3.0, 4.0, 5.0});
  math::point p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 5.0);
  BOOST_CHECK_CLOSE(p_cyl[1], 0.927295, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  math::point p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                            .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2[0], 3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2[1], 4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({-3.0, 4.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 5.0);
  BOOST_CHECK_CLOSE(p_cyl[1], 2.214297, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2[0], -3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2[1], 4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({-3.0, -4.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 5.0);
  BOOST_CHECK_CLOSE(p_cyl[1], -2.214297, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2[0], -3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2[1], -4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({3.0, -4.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 5.0);
  BOOST_CHECK_CLOSE(p_cyl[1], -0.927295, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_CLOSE(p_cart2[0], 3.0, 1e-3);
  BOOST_CHECK_CLOSE(p_cart2[1], -4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({0.0, 4.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 4.0);
  BOOST_CHECK_CLOSE(p_cyl[1], 1.5708, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2[0], 1e-5);
  BOOST_CHECK_CLOSE(p_cart2[1], 4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({0.0, -4.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 4.0);
  BOOST_CHECK_CLOSE(p_cyl[1], -1.5708, 1e-3);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2[0], 1e-5);
  BOOST_CHECK_CLOSE(p_cart2[1], -4.0, 1e-3);
  BOOST_CHECK(p_cart2[2] == 5.0);

  p_cart = math::point({0.0, 0.0, 5.0});
  p_cyl = math::point_typed(p_cart).to_type(math::point_typed::coordinate_type::cylindrical);
  BOOST_CHECK(p_cyl[0] == 0.0);
  BOOST_CHECK(p_cyl[1] == 0.0);
  BOOST_CHECK(p_cyl[2] == 5.0);
  p_cart2 = math::point_typed(p_cyl, math::point_typed::coordinate_type::cylindrical)
                .to_type(math::point_typed::coordinate_type::cartesian);
  BOOST_CHECK_SMALL(p_cart2[0], 1e-5);
  BOOST_CHECK_SMALL(p_cart2[1], 1e-5);
  BOOST_CHECK(p_cart2[2] == 5.0);
}

BOOST_AUTO_TEST_SUITE_END()
