#include <boost/test/unit_test.hpp>

#include "che/assembly_atom.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly_atom)

BOOST_AUTO_TEST_CASE(assembly_atom_ctor) {
  std::string chain_id("A");
  size_t chain_pos(10);
  che::assembly_atom a(chain_id, chain_pos, che::atom());
  BOOST_CHECK(a.get_chain_id() == chain_id);
  BOOST_CHECK(a.get_chain_pos() == chain_pos);
  BOOST_CHECK(a.get_identifier().empty());
}

BOOST_AUTO_TEST_SUITE_END()
