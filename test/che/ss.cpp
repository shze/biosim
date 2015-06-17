#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/ss.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ss)

BOOST_AUTO_TEST_CASE(ss_ctor) {
  che::ss s(0);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  s = che::ss(100);
  BOOST_CHECK(s.get_sequence().size() == 100);
  BOOST_CHECK(s.get_sses().empty());

  che::sequence<che::cchb_dssp> v;
  s = che::ss(v);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  v.emplace_back('H');
  s = che::ss(v);
  BOOST_CHECK(s.get_sequence().size() == 1);
  BOOST_CHECK(s.get_sses().size() == 1);

  std::set<che::cchb_dssp_interval> pool;
  s = che::ss(pool);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  pool.insert(che::cchb_dssp_interval(5, 10, che::cchb_dssp('H')));
  s = che::ss(pool);
  BOOST_CHECK(s.get_sequence().size() == 11); // sse position is 0-based, sequence is 0-4 C, 5-10 H, length is 11
  BOOST_CHECK(s.get_sses().size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
