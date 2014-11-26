#include <boost/test/unit_test.hpp>

#include "tools/string.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_string)

BOOST_AUTO_TEST_CASE(string_unique) {
  std::string s1, s2("a a a a");
  tools::char_function char_isspace = (int (*)(int)) & std::isspace;

  BOOST_CHECK(tools::get_unique_char_set_ignore(s1, char_isspace).empty());
  BOOST_CHECK(tools::get_unique_char_set_ignore(s2, char_isspace).size() == 1);

  BOOST_CHECK(tools::get_unique_char_set(s1).empty());
  BOOST_CHECK(tools::get_unique_char_set(s2).size() == 2);

  BOOST_CHECK(tools::get_unique_char_string_ignore(s1, char_isspace).empty());
  BOOST_CHECK(tools::get_unique_char_string_ignore(s2, char_isspace).size() == 1);

  BOOST_CHECK(tools::get_unique_char_string(s1).empty());
  BOOST_CHECK(tools::get_unique_char_string(s2).size() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
