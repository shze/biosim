#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/bond.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_bond)

BOOST_AUTO_TEST_CASE(bond_ctor) {
  che::bond b;
  BOOST_CHECK(b.get_type() == che::bond::value_type::sing);

  for(auto p : che::bond::get_value_type_map()) {
    che::bond b(p.first);
    BOOST_CHECK(b.get_type() == p.second);
  }
}

BOOST_AUTO_TEST_SUITE_END()
