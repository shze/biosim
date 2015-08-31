#include <boost/test/unit_test.hpp>

#include "tools/incrementor.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_incrementor)

BOOST_AUTO_TEST_CASE(incrementor_increment_fixed_length) {
  std::vector<std::vector<char>> alphabets;
  tools::incrementor<std::string> inc_empty(alphabets); // this should work and not throw any exception
  BOOST_REQUIRE_THROW(inc_empty.next("0"), std::out_of_range);
  BOOST_REQUIRE_THROW(inc_empty.next("A"), std::out_of_range);
  inc_empty.next("", 0);

  alphabets.emplace_back(std::vector<char>());
  BOOST_REQUIRE_THROW((tools::incrementor<std::string>(alphabets)), std::invalid_argument);

  tools::incrementor<std::string> inc({{'A', 'B'}, {'A', 'B'}});
  BOOST_REQUIRE_THROW(inc.next("A"), std::out_of_range);
  BOOST_CHECK(inc.next("AA") == "AB");
  BOOST_CHECK(inc.next("AB") == "BA");
  BOOST_CHECK(inc.next("BA") == "BB");
  BOOST_CHECK(inc.overflow() == false);
  inc.next("BB");
  BOOST_CHECK(inc.overflow() == true);

  BOOST_REQUIRE_THROW(inc.next("AC"), std::out_of_range);
  inc.next("CA"); // this should work as the last digit can be incremented and the 'C' digit is not checked
  BOOST_REQUIRE_THROW(inc.next("CB"), std::out_of_range);
  inc.next("CB", 0);
}

BOOST_AUTO_TEST_SUITE_END()
