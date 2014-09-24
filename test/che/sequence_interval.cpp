#include <boost/test/unit_test.hpp>

#include "che/sequence_interval.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_sequence_interval)

BOOST_AUTO_TEST_CASE(sequence_interval_ctor) {
  che::cchb_interval iv(math::interval<size_t>(1, 4), che::cchb("H"));
  BOOST_CHECK(iv.get_dimension().get_min() == 1);
  BOOST_CHECK(iv.get_dimension().get_max() == 4);
  BOOST_CHECK(iv.get_type().get_identifier() == "H");
}

BOOST_AUTO_TEST_CASE(sequence_interval_overlapping_subset) {
  std::set<che::cchb_interval> sses;
  sses.insert(che::cchb_interval(math::interval<size_t>(1, 4), che::cchb("H")));
  che::cchb_interval iv1(math::interval<size_t>(1, 4), che::cchb("H")),
      iv2(math::interval<size_t>(1, 4), che::cchb("C")), iv3(math::interval<size_t>(5, 5), che::cchb("H"));

  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv1, sses.begin(), sses.end()).size() == 1);
  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv2, sses.begin(), sses.end()).empty());
  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv3, sses.begin(), sses.end()).empty());

  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv1, sses.begin(), sses.end(), true).size() == 1);
  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv2, sses.begin(), sses.end(), true).size() == 1);
  BOOST_CHECK(che::cchb_interval::get_overlapping_subset(iv3, sses.begin(), sses.end(), true).empty());

  BOOST_CHECK(che::cchb_interval::overlaps(iv1, sses.begin(), sses.end()));
  BOOST_CHECK(che::cchb_interval::overlaps(iv2, sses.begin(), sses.end()) == false);
  BOOST_CHECK(che::cchb_interval::overlaps(iv3, sses.begin(), sses.end()) == false);

  BOOST_CHECK(che::cchb_interval::overlaps(iv1, sses.begin(), sses.end(), true));
  BOOST_CHECK(che::cchb_interval::overlaps(iv2, sses.begin(), sses.end(), true));
  BOOST_CHECK(che::cchb_interval::overlaps(iv3, sses.begin(), sses.end(), true) == false);
}

BOOST_AUTO_TEST_CASE(sequence_interval_less_dimension_min_max) {
  che::cchb_interval iv1(math::interval<size_t>(1, 4), che::cchb("E")),
      iv2(math::interval<size_t>(2, 4), che::cchb("E")), iv3(math::interval<size_t>(0, 4), che::cchb("E")),
      iv4(math::interval<size_t>(1, 4), che::cchb("C"));

  BOOST_CHECK(!(iv1 < iv1));
  BOOST_CHECK(iv1 < iv2);
  BOOST_CHECK(!(iv1 < iv3));
  BOOST_CHECK(!(iv1 < iv4));
  BOOST_CHECK(!(iv4 < iv1));
}

BOOST_AUTO_TEST_CASE(sequence_interval_less_dimension_max) {
  che::cchb_interval iv1(math::interval<size_t>(1, 4), che::cchb("E")),
      iv2(math::interval<size_t>(1, 4), che::cchb("H")), iv3(math::interval<size_t>(1, 3), che::cchb("E")),
      iv4(math::interval<size_t>(1, 5), che::cchb("E"));

  BOOST_CHECK(!(iv1 < iv1));
  BOOST_CHECK(!(iv1 < iv2));
  BOOST_CHECK(!(iv2 < iv1));
  BOOST_CHECK(!(iv1 < iv3));
  BOOST_CHECK(iv3 < iv1);
  BOOST_CHECK(iv1 < iv4);
  BOOST_CHECK(!(iv4 < iv1));
}

BOOST_AUTO_TEST_SUITE_END()
