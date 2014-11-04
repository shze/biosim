#include <boost/test/unit_test.hpp>

#include "math/interval_scheduler_maximize.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_interval_scheduler_maximize)

BOOST_AUTO_TEST_CASE(interval_scheduler_maximize) {
  std::set<math::interval<int>> intervals;
  math::interval_scheduler_maximize<math::interval<int>> s;
  s.schedule(intervals);
}

BOOST_AUTO_TEST_SUITE_END()
