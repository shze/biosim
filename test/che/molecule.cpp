#include <boost/test/unit_test.hpp>

#include "che/molecule.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_molecule)

BOOST_AUTO_TEST_CASE(molecule_ctor) {
  che::molecule m;
  BOOST_CHECK(m.get_storage().empty());
  BOOST_CHECK(m.get_identifier().empty());
  BOOST_CHECK(m.get_ps().size() == 0);
  BOOST_CHECK(m.get_ps().size() == m.get_length());

  m = che::molecule("storage", "id");
  BOOST_CHECK(m.get_storage() == "storage");
  BOOST_CHECK(m.get_identifier() == "id");

  che::ps p;
  p.emplace_back('A');
  m = che::molecule("storage", "id", p);
  BOOST_CHECK(m.get_storage() == "storage");
  BOOST_CHECK(m.get_identifier() == "id");
  BOOST_CHECK(m.get_ps().size() == 1);
  BOOST_CHECK(m.get_ps().size() == m.get_length());
  BOOST_CHECK(m.get_ps().front().get_identifier_char() == 'A');
}

BOOST_AUTO_TEST_CASE(molecule_ps) {
  che::molecule m;
  che::ps p;
  p.emplace_back('A');
  m.set_ps(p);
  BOOST_CHECK(m.get_ps().size() == 1);
  BOOST_CHECK(m.get_ps().size() == m.get_length());
  BOOST_CHECK(m.get_ps().front().get_identifier_char() == 'A');
}

BOOST_AUTO_TEST_SUITE_END()
