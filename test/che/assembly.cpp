#include <boost/test/unit_test.hpp>

#include "che/assembly.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_assembly)

BOOST_AUTO_TEST_CASE(assembly_ctor) {
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().empty());
  BOOST_CHECK(a.get_ensemble_size("A") == 0);
  BOOST_CHECK(a.has_ensemble("A") == false);
  BOOST_REQUIRE_THROW(a.get_ensemble("A"), std::out_of_range);
  BOOST_REQUIRE_THROW(a.get_first_structure("A"), std::out_of_range);
  BOOST_CHECK(a.get_first_structures().empty() == true);

  std::string storage = "some_storage";
  std::string id = "some_id";
  che::structure s(storage, id);
  a = che::assembly(s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble_size("A") == 1);
  BOOST_CHECK(a.has_ensemble("A") == true);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);
  BOOST_CHECK(a.get_first_structure("A").get_identifier() == id);
  BOOST_CHECK(a.get_first_structures().size() == 1);

  std::list<che::structure> structures{s, che::structure()};
  che::structure_ensemble e(structures.begin(), structures.end());
  a = che::assembly(e);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble_size("A") == 2);
  BOOST_CHECK(a.has_ensemble("A") == true);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 2);
  BOOST_CHECK(a.get_first_structure("A").get_identifier() == id);
  BOOST_CHECK(a.get_first_structures().size() == 1);
}

BOOST_AUTO_TEST_CASE(assembly_add_set) {
  std::string storage = "some_storage";
  std::string id = "some_id";
  che::structure s(storage, id);
  che::assembly a;
  BOOST_CHECK(a.get_chain_id_list().size() == 0);
  BOOST_CHECK(a.add(s) == "A"); // check returned chain_id
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);

  a.set("A", s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);
  BOOST_CHECK(a.has_ensemble("B") == false);

  a.set("C", s);
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);
  BOOST_CHECK(a.get_ensemble("C").get_samples().size() == 1);
  BOOST_CHECK(a.has_ensemble("B") == false);
  BOOST_CHECK(a.has_ensemble("C") == true);
  BOOST_CHECK(a.has_ensemble("D") == false);

  std::list<che::structure> structures{s, che::structure()};
  che::structure_ensemble e(structures.begin(), structures.end());
  a = che::assembly();
  BOOST_CHECK(a.get_chain_id_list().size() == 0);
  BOOST_CHECK(a.add(e) == "A"); // check returned chain_id
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 2);

  a.set("A", s);
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);
  BOOST_CHECK(a.has_ensemble("B") == false);

  a.set("C", e);
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_ensemble("A").get_samples().size() == 1);
  BOOST_CHECK(a.get_ensemble("C").get_samples().size() == 2);
  BOOST_CHECK(a.has_ensemble("B") == false);
  BOOST_CHECK(a.has_ensemble("C") == true);
  BOOST_CHECK(a.has_ensemble("D") == false);
}

BOOST_AUTO_TEST_SUITE_END()
