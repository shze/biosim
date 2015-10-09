#include <boost/test/unit_test.hpp>

#include "math/algo/point_distance.h" // header to test
#include "tools/log.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_point_distance)

BOOST_AUTO_TEST_CASE(point_distance) {
  math::point p0({0, 0, 0}), p3({3, 0, 0}), p32({3, 2, 0}), p321({3, 2, 1});

  BOOST_CHECK(math::algo::minkowski_distance(p0, p0, 3.0) == 0.0);
  BOOST_CHECK(math::algo::minkowski_distance(p0, p3, 3.0) == 3.0);
  BOOST_CHECK_CLOSE(math::algo::minkowski_distance(p0, p32, 3.0), 3.27107, 1e-3);
  BOOST_CHECK_CLOSE(math::algo::minkowski_distance(p0, p321, 3.0), 3.30193, 1e-3);

  BOOST_CHECK(math::algo::minkowski_distance_power(p0, p0, 3.0) == 0.0);
  BOOST_CHECK(math::algo::minkowski_distance_power(p0, p3, 3.0) == 27.0);
  BOOST_CHECK(math::algo::minkowski_distance_power(p0, p32, 3.0) == 35.0);
  BOOST_CHECK(math::algo::minkowski_distance_power(p0, p321, 3.0) == 36.0);

  BOOST_CHECK(math::algo::euclidean_distance(p0, p0) == 0.0);
  BOOST_CHECK(math::algo::euclidean_distance(p0, p3) == 3.0);
  BOOST_CHECK_CLOSE(math::algo::euclidean_distance(p0, p32), 3.60555, 1e-3);
  BOOST_CHECK_CLOSE(math::algo::euclidean_distance(p0, p321), 3.74166, 1e-3);

  BOOST_CHECK(math::algo::euclidean_distance_power(p0, p0) == 0.0);
  BOOST_CHECK(math::algo::euclidean_distance_power(p0, p3) == 9.0);
  BOOST_CHECK(math::algo::euclidean_distance_power(p0, p32) == 13.0);
  BOOST_CHECK(math::algo::euclidean_distance_power(p0, p321) == 14.0);

  BOOST_CHECK(math::algo::manhattan_distance(p0, p0) == 0.0);
  BOOST_CHECK(math::algo::manhattan_distance(p0, p3) == 3.0);
  BOOST_CHECK(math::algo::manhattan_distance(p0, p32) == 5.0);
  BOOST_CHECK(math::algo::manhattan_distance(p0, p321) == 6.0);

  std::array<double, 1> p1d0{0}, p1d2{2}, p1dinf{std::numeric_limits<double>::infinity()};
  BOOST_CHECK(math::algo::manhattan_distance(p1d0, p1d0) == 0);
  BOOST_CHECK(math::algo::manhattan_distance(p1d0, p1d2) == 2);
  BOOST_CHECK(math::algo::euclidean_distance(p1d0, p1d0) == 0);
  BOOST_CHECK(math::algo::euclidean_distance(p1d0, p1d2) == 2);
  BOOST_CHECK(std::isinf(math::algo::manhattan_distance(p1d0, p1dinf)) == true);
  BOOST_CHECK(std::isinf(math::algo::euclidean_distance(p1d0, p1dinf)) == true);

  std::array<double, 2> p2d0{0, 0}, p2d2{2, 0}, p2d21{2, 1};
  BOOST_CHECK(math::algo::manhattan_distance(p2d0, p2d0) == 0);
  BOOST_CHECK(math::algo::manhattan_distance(p2d0, p2d2) == 2);
  BOOST_CHECK(math::algo::manhattan_distance(p2d0, p2d21) == 3);
  BOOST_CHECK(math::algo::euclidean_distance(p2d0, p2d0) == 0);
  BOOST_CHECK(math::algo::euclidean_distance(p2d0, p2d2) == 2);
  BOOST_CHECK_CLOSE(math::algo::euclidean_distance(p2d0, p2d21), 2.23607, 1e-3);

  std::array<double, 4> p4d0{0, 0, 0, 0}, p4d4321{4, 3, 2, 1};
  BOOST_CHECK(math::algo::manhattan_distance(p4d0, p4d4321) == 10);
  BOOST_CHECK_CLOSE(math::algo::euclidean_distance(p4d0, p4d4321), 5.47723, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
