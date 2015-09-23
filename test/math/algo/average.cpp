#include <boost/test/unit_test.hpp>

#include "math/algo/average.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_average)

BOOST_AUTO_TEST_CASE(average) {
  math::algo::cma avg;
  BOOST_CHECK(avg.get_n() == 0);
  BOOST_CHECK(std::isnan(avg.get_average()) == true);

  double v1 = 1.5;
  avg.observe(v1);
  BOOST_CHECK(avg.get_n() == 1);
  BOOST_CHECK(avg.get_average() == v1);

  avg.observe(v1);
  BOOST_CHECK(avg.get_n() == 2);
  BOOST_CHECK(avg.get_average() == v1);

  double v2 = 4.5;
  avg.observe(v2);
  BOOST_CHECK(avg.get_n() == 3);
  BOOST_CHECK_CLOSE(avg.get_average(), ((1.5 * 2 + 4.5) / 3), 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
