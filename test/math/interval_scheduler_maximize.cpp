#include <boost/test/unit_test.hpp>

#include "math/interval_scheduler_maximize.h" // header to test
#include "che/sequence_interval.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_interval_scheduler_maximize)

BOOST_AUTO_TEST_CASE(interval_scheduler_maximize_interval) {
  std::set<math::interval<int>> intervals;
  math::interval_scheduler_maximize<math::interval<int>> s;
  BOOST_CHECK(s.schedule(intervals).empty());

  math::interval<int> iv1(1, 5), iv2(2, 4), iv3(2, 6), iv4(5, 10), iv5(6, 10), iv6(7, 10);

  intervals.insert(iv1);
  BOOST_CHECK(s.schedule(intervals).size() == 1);

  intervals.insert(iv2);
  intervals.insert(iv3);
  BOOST_CHECK(s.schedule(intervals).size() == 1);

  intervals.insert(iv4);
  BOOST_CHECK(s.schedule(intervals).size() == 2);

  intervals.insert(iv5);
  intervals.insert(iv6);
  BOOST_CHECK(s.schedule(intervals).size() == 2);
}

BOOST_AUTO_TEST_CASE(interval_scheduler_maximize_sequence_interval) {
  std::set<che::cchb_dssp_interval> intervals;
  math::interval_scheduler_maximize<che::cchb_dssp_interval> s;
  BOOST_CHECK(s.schedule(intervals).empty());

  che::cchb_dssp_interval iv1(1, 5, che::cchb_dssp('H')), iv2(2, 4, che::cchb_dssp('H')),
      iv3(4, 6, che::cchb_dssp('E'));
  intervals.insert(iv1);
  intervals.insert(iv2);
  BOOST_CHECK(s.schedule(intervals).size() == 1);
  
  intervals.insert(iv3);
  BOOST_CHECK(s.schedule(intervals).size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
