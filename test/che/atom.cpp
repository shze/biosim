#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/atom.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_atom)

BOOST_AUTO_TEST_CASE(atom_ctor) {
  che::atom a;
  BOOST_CHECK(a.get_identifier().empty());
  BOOST_CHECK(a.get_position() == math::point());

  std::string atom_id("xx");
  a = che::atom(atom_id);
  BOOST_CHECK(a.get_identifier() == atom_id);
  BOOST_CHECK(a.get_position() == math::point());

  math::point p({1, 2, 3}, math::point::coordinate_type::cylindrical);
  a = che::atom(atom_id, p);
  BOOST_CHECK(a.get_identifier() == atom_id);
  BOOST_CHECK(a.get_position() == p);
}

BOOST_AUTO_TEST_CASE(atom_cmp_operator) {
  che::atom a;
  che::atom a2("x");
  BOOST_CHECK(!(a < a));
  BOOST_CHECK(a < a2);
}

BOOST_AUTO_TEST_SUITE_END()
