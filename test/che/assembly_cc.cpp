#include <boost/test/unit_test.hpp>

#include "che/assembly_cc.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly_cc)

BOOST_AUTO_TEST_CASE(assembly_cc_ctor) {
  std::string chain_id("A");
  size_t chain_pos(10);
  std::string cc_identifier("ALA");
  che::assembly_cc c(chain_id, chain_pos, che::cc(cc_identifier));
  BOOST_CHECK(c.get_chain_id() == chain_id);
  BOOST_CHECK(c.get_chain_pos() == chain_pos);
  BOOST_CHECK(c.get_identifier() == cc_identifier);
}

BOOST_AUTO_TEST_SUITE_END()
