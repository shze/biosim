#include <boost/test/unit_test.hpp>

#include "che/qs.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_qs)

BOOST_AUTO_TEST_CASE(qs_ctor) {
  che::qs q;
  BOOST_CHECK(q.get_chain_id_list().empty());
}

BOOST_AUTO_TEST_CASE(qs_add) {
  che::ts t;
  che::qs q;
  BOOST_CHECK(q.add(t) == "A"); // check returned chain_id
  BOOST_CHECK(q.get_chain_id_list().size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
