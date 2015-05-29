#include <boost/test/unit_test.hpp>

#include "tools/incrementor.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_incrementor)

BOOST_AUTO_TEST_CASE(incrementor_increment_fixed_length) {
  std::vector<std::vector<char>> empty_alphabets;
  tools::incrementor<std::string> inc_empty(empty_alphabets);
  BOOST_REQUIRE_THROW(inc_empty.next("0"), std::overflow_error);
  BOOST_REQUIRE_THROW(inc_empty.next("A"), std::overflow_error);

  empty_alphabets.push_back({});
  BOOST_REQUIRE_THROW(tools::incrementor<std::string> inc_empty2(empty_alphabets), std::invalid_argument);

  tools::incrementor<std::string> inc({{'A', 'B'}, {'A', 'B'}});
  BOOST_REQUIRE_THROW(inc.next("A"), std::out_of_range);
  BOOST_CHECK(inc.next("AA") == "AB");
  BOOST_CHECK(inc.next("AB") == "BA");
  BOOST_CHECK(inc.next("BA") == "BB");
  BOOST_REQUIRE_THROW(inc.next("BB"), std::overflow_error);
  BOOST_REQUIRE_THROW(inc.next("CA"), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(incrementor_increment_flexible_length) {
  std::vector<char> empty_alphabet;
  BOOST_REQUIRE_THROW(tools::incrementor<std::string> inc_empty(empty_alphabet), std::invalid_argument);

  tools::incrementor<std::string> inc_letter;
  BOOST_CHECK(inc_letter.next("A") == "B");
  BOOST_CHECK(inc_letter.next("B") == "C");
  BOOST_CHECK(inc_letter.next("C", 26) == "AC");
  BOOST_REQUIRE_THROW(inc_letter.next("0"), std::out_of_range);

  tools::incrementor<std::string> inc_letter2({'0', '1', '2'});
  BOOST_CHECK(inc_letter2.next("0") == "1");
  BOOST_CHECK(inc_letter2.next("1") == "2");
  BOOST_CHECK(inc_letter2.next("2") == "10");
  BOOST_CHECK(inc_letter2.next("10", 4) == "21");
  BOOST_CHECK(inc_letter2.next("21", 4) == "102");
  BOOST_CHECK(inc_letter2.next("0", 11) == "102");

  tools::incrementor<std::vector<size_t>> inc_digit({0, 1, 2});
  BOOST_CHECK(inc_digit.next({0}) == std::vector<size_t>({1}));
  BOOST_CHECK(inc_digit.next({1}) == std::vector<size_t>({2}));
  BOOST_CHECK(inc_digit.next({2}) == std::vector<size_t>({1, 0}));
  BOOST_CHECK(inc_digit.next({1, 0}, 4) == std::vector<size_t>({2, 1}));
  BOOST_CHECK(inc_digit.next({2, 1}, 4) == std::vector<size_t>({1, 0, 2}));
  BOOST_CHECK(inc_digit.next({0}, 11) == std::vector<size_t>({1, 0, 2}));
}

BOOST_AUTO_TEST_SUITE_END()
