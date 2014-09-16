#include <boost/test/unit_test.hpp>
#include "tools/log.h"

#include "tools/incrementor.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_incrementor)

BOOST_AUTO_TEST_CASE(number_increment) {
  tools::incrementor inc("012");
  std::string s = "0";
  s = inc.next(s);
  BOOST_CHECK(s == "1");
  s = inc.next(s);
  BOOST_CHECK(s == "2");
  s = inc.next(s);
  BOOST_CHECK(s == "10");
  s = inc.next(s, 4);
  BOOST_CHECK(s == "21");
  s = inc.next(s, 4);
  BOOST_CHECK(s == "102");
}

BOOST_AUTO_TEST_CASE(letter_increment) {
  tools::incrementor inc;
  std::string s = "A";
  s = inc.next(s);
  BOOST_CHECK(s == "B");
  s = inc.next(s);
  BOOST_CHECK(s == "C");
  s = inc.next(s, 26);
  BOOST_CHECK(s == "AC");
}

BOOST_AUTO_TEST_CASE(increment_exception) {
  tools::incrementor inc;
  std::string s = "0";
  BOOST_REQUIRE_THROW(inc.next(s), std::range_error);
}

BOOST_AUTO_TEST_SUITE_END()
