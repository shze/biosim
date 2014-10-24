#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/ss.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ss)

BOOST_AUTO_TEST_CASE(ss_ctor) {
  che::sequence<che::cchb_dssp> v;
  che::ss s(v);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  v.emplace_back('H');
  s = che::ss(v);
  BOOST_CHECK(s.get_sequence().size() == 1);
  BOOST_CHECK(s.get_sses().empty());

  std::set<che::cchb_dssp_interval> pool;
  s = che::ss(pool);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().empty());

  pool.insert(che::cchb_dssp_interval(math::interval<size_t>(5, 10), che::cchb_dssp('H')));
  s = che::ss(pool);
  BOOST_CHECK(s.get_sequence().empty());
  BOOST_CHECK(s.get_sses().size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
