#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/ss.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ss)

BOOST_AUTO_TEST_CASE(ss_ctor) {
  std::vector<che::cchb> v;
  che::ss s(v);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  v.emplace_back("H");
  s = che::ss(v);
  BOOST_CHECK(s.get_sequence().size() == 1);
  BOOST_CHECK(s.get_sses().empty());
}

BOOST_AUTO_TEST_SUITE_END()
