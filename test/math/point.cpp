#include <boost/test/unit_test.hpp>

#include "math/point.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_point)

BOOST_AUTO_TEST_CASE(point_ctor) {
  math::point p;
  BOOST_CHECK(p[0] == 0.0);
  BOOST_CHECK(p[1] == 0.0);
  BOOST_CHECK(p[2] == 0.0);
  BOOST_CHECK(p == p);

  p = math::point({1.0, 2.0, 3.0});
  BOOST_CHECK(p[0] == 1.0);
  BOOST_CHECK(p[1] == 2.0);
  BOOST_CHECK(p[2] == 3.0);
  BOOST_CHECK(p == p);
}

BOOST_AUTO_TEST_SUITE_END()
