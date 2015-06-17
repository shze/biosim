#include <boost/test/unit_test.hpp>

#include "che/assembly.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly)

BOOST_AUTO_TEST_CASE(assembly_ctor) {
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().empty());
  BOOST_CHECK(a.has_molecule("A") == false);
  BOOST_REQUIRE_THROW(a.get_molecule("A"), std::out_of_range);

  che::molecule m;
  a = che::assembly(m);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.has_molecule("A") == true);
  BOOST_CHECK(a.get_molecule("A").get_ps().empty());
}

BOOST_AUTO_TEST_CASE(assembly_add_set) {
  che::molecule m;
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().size() == 0);
  BOOST_CHECK(a.add(m) == "A"); // check returned chain_id
  BOOST_CHECK(a.get_chain_id_list().size() == 1);

  a.set("A", m);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.has_molecule("B") == false);

  a.set("C", m);
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.has_molecule("B") == false);
  BOOST_CHECK(a.has_molecule("C") == true);
  BOOST_CHECK(a.has_molecule("D") == false);
}

BOOST_AUTO_TEST_SUITE_END()
