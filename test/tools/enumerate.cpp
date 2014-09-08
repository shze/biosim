#include <boost/test/unit_test.hpp>

// need to define this before include
std::string identifier(std::string const &__s) { return __s; }
std::string identifier(int const &__i) { return std::to_string(__i); }

#include "tools/enumerate.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_enumerate)

BOOST_AUTO_TEST_CASE(string_enum) {
  using string_enum = tools::enumerate<std::string>;
  string_enum::add("x");
  string_enum::add("second");
  string_enum::add("3");
  string_enum::add("3");
  string_enum::add("3");
  BOOST_CHECK(string_enum::get_instances().size() == 3);
  for(auto p : string_enum::get_instances()) {
    BOOST_CHECK(p.first == p.second.get_object());
  }
}

BOOST_AUTO_TEST_CASE(int_enum) {
  using int_enum = tools::enumerate<int>;
  int_enum::add(-3);
  int_enum::add(0);
  int_enum::add(546657468);
  BOOST_CHECK(int_enum::get_instances().size() == 3);
}

struct fancystring {
  fancystring(std::string __s) : _s(__s) {}
  std::string _s;
};

std::string identifier(fancystring const &__f) { return __f._s; }

BOOST_AUTO_TEST_CASE(struct_enum) {
  using fancystring_enum = tools::enumerate<fancystring>;
  fancystring_enum::add(fancystring("f"));
  std::string s = "f2";
  fancystring_enum::add(s);
  BOOST_CHECK(fancystring_enum::get_instances().size() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
