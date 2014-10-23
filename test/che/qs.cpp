#include <boost/test/unit_test.hpp>

#include "che/qs.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_qs)

BOOST_AUTO_TEST_CASE(qs_ctor) {
  che::qs q;
  BOOST_CHECK(q.get_chain_id_list().empty());
  BOOST_REQUIRE_THROW(q.get_ts("A"), std::out_of_range);

  che::ts t;
  std::vector<che::cchb_dssp> v;
  che::ss s(v);
  q = che::qs(t, s);
  BOOST_CHECK(q.get_chain_id_list().size() == 1);
  BOOST_CHECK(q.get_ts("A").get_ps().empty());
  BOOST_CHECK(q.get_ss("A").get_sequence().empty());
}

BOOST_AUTO_TEST_CASE(qs_add) {
  che::ts t;
  che::qs q;
  BOOST_CHECK(q.add(t) == "A"); // check returned chain_id
  BOOST_CHECK(q.get_chain_id_list().size() == 1);
  BOOST_CHECK(q.get_ts("A").get_length() == 0);
  BOOST_REQUIRE_THROW(q.get_ss("A"), std::out_of_range); // no ss created by add(ts)

  std::vector<che::cchb_dssp> v;
  che::ss s(v);
  BOOST_CHECK(q.add(t, s) == "B");
  BOOST_CHECK(q.get_chain_id_list().size() == 2);
  BOOST_CHECK(q.get_ts("B").get_length() == 0);
  BOOST_CHECK(q.get_ss("B").get_sequence().empty());
}

BOOST_AUTO_TEST_SUITE_END()
