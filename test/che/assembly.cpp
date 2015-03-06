#include <boost/test/unit_test.hpp>

#include "che/assembly.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly)

BOOST_AUTO_TEST_CASE(assembly_ctor) {
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().empty());
  BOOST_REQUIRE_THROW(a.get_molecule("A"), std::out_of_range);

  che::molecule m;
  che::sequence<che::cchb_dssp> v;
  che::ss s(v);
  a = che::assembly(m, s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_molecule("A").get_ps().empty());
  BOOST_CHECK(a.get_ss("A").get_sequence().empty());
}

BOOST_AUTO_TEST_CASE(assembly_add) {
  che::molecule m;
  che::assembly a;
  BOOST_CHECK(a.add(m) == "A"); // check returned chain_id
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_molecule("A").get_length() == 0);
  BOOST_REQUIRE_THROW(a.get_ss("A"), std::out_of_range); // no ss created by add(ts)

  che::sequence<che::cchb_dssp> v;
  che::ss s(v);
  BOOST_CHECK(a.add(m, s) == "B");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_molecule("B").get_length() == 0);
  BOOST_CHECK(a.get_ss("B").get_sequence().empty());
}

BOOST_AUTO_TEST_SUITE_END()
