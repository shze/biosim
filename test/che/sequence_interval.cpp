#include <boost/test/unit_test.hpp>

#include "che/sequence_interval.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_sequence_interval)

BOOST_AUTO_TEST_CASE(sequence_interval_ctor) {
  che::cchb_dssp_interval iv(1, 4, che::cchb_dssp('H'));
  BOOST_CHECK(iv.get_min() == 1);
  BOOST_CHECK(iv.get_max() == 4);
  BOOST_CHECK(iv.get_type().get_identifier() == 'H');
}

BOOST_AUTO_TEST_CASE(sequence_interval_overlapping_subset) {
  std::set<che::cchb_dssp_interval> sses;
  sses.insert(che::cchb_dssp_interval(1, 4, che::cchb_dssp('H')));
  che::cchb_dssp_interval iv1(1, 4, che::cchb_dssp('H')), iv2(1, 4, che::cchb_dssp('C')),
      iv3(5, 5, che::cchb_dssp('H'));

  BOOST_CHECK(iv1.get_overlapping_subset(sses.begin(), sses.end()).size() == 1);
  BOOST_CHECK(iv2.get_overlapping_subset(sses.begin(), sses.end()).size() == 1);
  BOOST_CHECK(iv3.get_overlapping_subset(sses.begin(), sses.end()).empty());

  BOOST_CHECK(iv1.get_overlapping_subset(sses.begin(), sses.end(), true).size() == 1);
  BOOST_CHECK(iv2.get_overlapping_subset(sses.begin(), sses.end(), true).size() == 1);
  BOOST_CHECK(iv3.get_overlapping_subset(sses.begin(), sses.end(), true).empty());

  BOOST_CHECK(iv1.get_overlapping_subset(sses.begin(), sses.end(), false).size() == 1);
  BOOST_CHECK(iv2.get_overlapping_subset(sses.begin(), sses.end(), false).empty());
  BOOST_CHECK(iv3.get_overlapping_subset(sses.begin(), sses.end(), false).empty());

  BOOST_CHECK(iv1.overlaps(sses.begin(), sses.end()));
  BOOST_CHECK(iv2.overlaps(sses.begin(), sses.end()));
  BOOST_CHECK(iv3.overlaps(sses.begin(), sses.end()) == false);

  BOOST_CHECK(iv1.overlaps(sses.begin(), sses.end(), true));
  BOOST_CHECK(iv2.overlaps(sses.begin(), sses.end(), true));
  BOOST_CHECK(iv3.overlaps(sses.begin(), sses.end(), true) == false);

  BOOST_CHECK(iv1.overlaps(sses.begin(), sses.end(), false));
  BOOST_CHECK(iv2.overlaps(sses.begin(), sses.end(), false) == false);
  BOOST_CHECK(iv3.overlaps(sses.begin(), sses.end(), false) == false);
}

BOOST_AUTO_TEST_CASE(sequence_interval_less_dimension_min_max) {
  che::cchb_dssp_interval iv1(1, 4, che::cchb_dssp('E')), iv2(2, 4, che::cchb_dssp('E')),
      iv3(0, 4, che::cchb_dssp('E')), iv4(1, 4, che::cchb_dssp('C'));

  BOOST_CHECK(!(iv1 < iv1));
  BOOST_CHECK(iv1 < iv2);
  BOOST_CHECK(!(iv1 < iv3));
  BOOST_CHECK(!(iv1 < iv4));
  BOOST_CHECK(!(iv4 < iv1));
}

BOOST_AUTO_TEST_CASE(sequence_interval_less_dimension_max) {
  che::cchb_dssp_interval iv1(1, 4, che::cchb_dssp('E')), iv2(1, 4, che::cchb_dssp('H')),
      iv3(1, 3, che::cchb_dssp('E')), iv4(1, 5, che::cchb_dssp('E'));

  BOOST_CHECK(!(iv1 < iv1));
  BOOST_CHECK(!(iv1 < iv2));
  BOOST_CHECK(!(iv2 < iv1));
  BOOST_CHECK(!(iv1 < iv3));
  BOOST_CHECK(iv3 < iv1);
  BOOST_CHECK(iv1 < iv4);
  BOOST_CHECK(!(iv4 < iv1));
}

BOOST_AUTO_TEST_SUITE_END()
