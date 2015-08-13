#include <boost/test/unit_test.hpp>

#include "che/assembly.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly)

BOOST_AUTO_TEST_CASE(assembly_ctor) {
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().empty());
  BOOST_CHECK(a.has_structure("A") == false);
  BOOST_REQUIRE_THROW(a.get_structure("A"), std::out_of_range);

  che::structure s;
  a = che::assembly(s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.has_structure("A") == true);
  BOOST_CHECK(a.get_structure("A").get_ps().empty());
}

BOOST_AUTO_TEST_CASE(assembly_add_set) {
  che::structure s;
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().size() == 0);
  BOOST_CHECK(a.add(s) == "A"); // check returned chain_id
  BOOST_CHECK(a.get_chain_id_list().size() == 1);

  a.set("A", s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.has_structure("B") == false);

  a.set("C", s);
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.has_structure("B") == false);
  BOOST_CHECK(a.has_structure("C") == true);
  BOOST_CHECK(a.has_structure("D") == false);
}

BOOST_AUTO_TEST_SUITE_END()
