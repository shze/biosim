#include <boost/test/unit_test.hpp>

#include "math/interval.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_interval)

BOOST_AUTO_TEST_CASE(int_interval) {
  math::interval<int> iv;
  BOOST_CHECK(iv.get_min() == -std::numeric_limits<int>::max());
  BOOST_CHECK(iv.get_max() == std::numeric_limits<int>::max());

  iv = math::interval<int>(-1, 1);
  BOOST_CHECK(iv.get_min() == -1);
  BOOST_CHECK(iv.get_max() == 1);
  BOOST_CHECK(iv.get_length() == 2);

  iv.set_min(3);
  BOOST_CHECK(iv.get_min() == 3);
  BOOST_CHECK(iv.get_max() == 3);

  iv.set_max(2);
  BOOST_CHECK(iv.get_min() == 2);
  BOOST_CHECK(iv.get_max() == 2);

  BOOST_CHECK(math::interval<int>::get_epsilon() == 1);
}

BOOST_AUTO_TEST_CASE(size_t_interval) {
  math::interval<size_t> iv;
  BOOST_CHECK(iv.get_min() == std::numeric_limits<size_t>::min());
  BOOST_CHECK(iv.get_max() == std::numeric_limits<size_t>::max());

  BOOST_CHECK(math::interval<size_t>::get_epsilon() == 1);
}

BOOST_AUTO_TEST_CASE(double_interval) {
  math::interval<double> iv;
  BOOST_CHECK(iv.get_min() == -std::numeric_limits<double>::max());
  BOOST_CHECK(iv.get_max() == std::numeric_limits<double>::max());

  iv = math::interval<double>(-1.0, 1.0);
  BOOST_CHECK(iv.get_min() == -1.0);
  BOOST_CHECK(iv.get_max() == 1.0);
  BOOST_CHECK(iv.get_length() == 2.0);

  iv.set_min(3.0);
  BOOST_CHECK(iv.get_min() == 3.0);
  BOOST_CHECK(iv.get_max() == 3.0);

  iv.set_max(2.0);
  BOOST_CHECK(iv.get_min() == 2.0);
  BOOST_CHECK(iv.get_max() == 2.0);

  BOOST_CHECK(math::interval<double>::get_epsilon() == std::numeric_limits<double>::epsilon());
}

BOOST_AUTO_TEST_CASE(interval_overlap) {
  math::interval<int> iv1(1, 3), iv2(2, 5), iv3(3, 5), iv4(4, 5), iv5(5, 5), iv6(1, 5), iv7(0, 5);

  BOOST_CHECK(iv1.is_continuous(iv2) == true);
  BOOST_CHECK(iv2.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv2) == true);
  BOOST_CHECK(iv2.overlaps(iv1) == true);

  BOOST_CHECK(iv1.is_continuous(iv3) == true);
  BOOST_CHECK(iv3.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv3) == true);
  BOOST_CHECK(iv3.overlaps(iv1) == true);

  BOOST_CHECK(iv1.is_continuous(iv4) == true);
  BOOST_CHECK(iv4.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv4) == false);
  BOOST_CHECK(iv4.overlaps(iv1) == false);

  BOOST_CHECK(iv1.is_continuous(iv5) == false);
  BOOST_CHECK(iv5.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv5) == false);
  BOOST_CHECK(iv5.overlaps(iv1) == false);

  BOOST_CHECK(iv1.is_continuous(iv6) == false);
  BOOST_CHECK(iv6.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv6) == true);
  BOOST_CHECK(iv6.overlaps(iv1) == true);

  BOOST_CHECK(iv1.is_continuous(iv7) == false);
  BOOST_CHECK(iv7.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.overlaps(iv7) == true);
  BOOST_CHECK(iv7.overlaps(iv1) == true);

  math::interval<int> iv8(-1, 2), iv9(-1, 1), iv10(-1, 0), iv11(-1, -1), iv12(-1, 3);

  BOOST_CHECK(iv8.is_continuous(iv1) == true);
  BOOST_CHECK(iv1.is_continuous(iv8) == false);
  BOOST_CHECK(iv1.overlaps(iv8) == true);
  BOOST_CHECK(iv8.overlaps(iv1) == true);

  BOOST_CHECK(iv9.is_continuous(iv1) == true);
  BOOST_CHECK(iv1.is_continuous(iv9) == false);
  BOOST_CHECK(iv1.overlaps(iv9) == true);
  BOOST_CHECK(iv9.overlaps(iv1) == true);

  BOOST_CHECK(iv10.is_continuous(iv1) == true);
  BOOST_CHECK(iv1.is_continuous(iv10) == false);
  BOOST_CHECK(iv1.overlaps(iv10) == false);
  BOOST_CHECK(iv10.overlaps(iv1) == false);

  BOOST_CHECK(iv11.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.is_continuous(iv11) == false);
  BOOST_CHECK(iv1.overlaps(iv11) == false);
  BOOST_CHECK(iv11.overlaps(iv1) == false);

  BOOST_CHECK(iv12.is_continuous(iv1) == false);
  BOOST_CHECK(iv1.is_continuous(iv12) == false);
  BOOST_CHECK(iv1.overlaps(iv12) == true);
  BOOST_CHECK(iv12.overlaps(iv1) == true);
}

BOOST_AUTO_TEST_CASE(interval_less_min_max) {
  math::interval<int> iv1(1, 3), iv2(1, 4), iv3(2, 3);

  BOOST_CHECK(math::interval<int>::less_min_max(iv1, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_min_max(iv1, iv2) == true);
  BOOST_CHECK(math::interval<int>::less_min_max(iv2, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_min_max(iv1, iv3) == true);
  BOOST_CHECK(math::interval<int>::less_min_max(iv3, iv1) == false);
}

BOOST_AUTO_TEST_CASE(interval_less_length) {
  math::interval<int> iv1(1, 3), iv2(1, 4), iv3(2, 4);

  BOOST_CHECK(math::interval<int>::less_length(iv1, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_length(iv1, iv2) == true);
  BOOST_CHECK(math::interval<int>::less_length(iv2, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_length(iv1, iv3) == false);
  BOOST_CHECK(math::interval<int>::less_length(iv3, iv1) == false);
}

BOOST_AUTO_TEST_CASE(interval_less_max) {
  math::interval<int> iv1(1, 3), iv2(1, 4), iv3(2, 4);

  BOOST_CHECK(math::interval<int>::less_max(iv1, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_max(iv1, iv2) == true);
  BOOST_CHECK(math::interval<int>::less_max(iv2, iv1) == false);

  BOOST_CHECK(math::interval<int>::less_max(iv2, iv3) == false);
  BOOST_CHECK(math::interval<int>::less_max(iv3, iv2) == false);
}

BOOST_AUTO_TEST_CASE(interval_equal_min_max) {
  math::interval<int> iv1(1, 3), iv2(1, 4), iv3(2, 4);
  BOOST_CHECK(math::interval<int>::equal_min_max(iv1, iv1));
  BOOST_CHECK(!math::interval<int>::equal_min_max(iv1, iv2));
  BOOST_CHECK(!math::interval<int>::equal_min_max(iv1, iv3));
}

BOOST_AUTO_TEST_CASE(interval_equal_length) {
  math::interval<int> iv1(1, 3), iv2(1, 4), iv3(2, 4);
  BOOST_CHECK(math::interval<int>::equal_length(iv1, iv1));
  BOOST_CHECK(!math::interval<int>::equal_length(iv1, iv2));
  BOOST_CHECK(math::interval<int>::equal_length(iv1, iv3));
}

BOOST_AUTO_TEST_CASE(interval_merge) {
  math::interval<int> iv1(5, 6), iv2(4, 7), iv3(4, 4), iv4(3, 3), iv5(6, 7), iv6(7, 7), iv7(8, 8);
  BOOST_CHECK(math::interval<int>::merge(iv1, iv2).size() == 1);
  BOOST_CHECK(math::interval<int>::merge(iv1, iv3).size() == 1);
  BOOST_CHECK(math::interval<int>::merge(iv1, iv4).empty());
  BOOST_CHECK(math::interval<int>::merge(iv1, iv5).size() == 1);
  BOOST_CHECK(math::interval<int>::merge(iv1, iv6).size() == 1);
  BOOST_CHECK(math::interval<int>::merge(iv1, iv7).empty());
}

BOOST_AUTO_TEST_SUITE_END()
