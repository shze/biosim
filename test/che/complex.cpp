#include <boost/test/unit_test.hpp>

#include "che/complex.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_complex)

BOOST_AUTO_TEST_CASE(complex_ctor) {
  che::complex c;
  BOOST_CHECK(c.get_chain_id_list().empty());
}

BOOST_AUTO_TEST_CASE(complex_add) {
  che::molecule m;
  che::complex c;
  BOOST_CHECK(c.add(m) == "A"); // check returned chain_id
  BOOST_CHECK(c.get_chain_id_list().size() == 1);
}

BOOST_AUTO_TEST_SUITE_END()
