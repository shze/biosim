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
  BOOST_CHECK(m.get_ss().defined() == false);

  m = che::molecule("storage", "id");
  BOOST_CHECK(m.get_storage() == "storage");
  BOOST_CHECK(m.get_identifier() == "id");
  BOOST_CHECK(m.get_length() == 0);
  BOOST_CHECK(m.get_ss().defined() == false);

  che::ps p;
  p.emplace_back('A');
  m = che::molecule("storage", "id", p);
  BOOST_CHECK(m.get_storage() == "storage");
  BOOST_CHECK(m.get_identifier() == "id");
  BOOST_CHECK(m.get_ps().size() == 1);
  BOOST_CHECK(m.get_ps().size() == m.get_length());
  BOOST_CHECK(m.get_ps().front().get_identifier_char() == 'A');
  BOOST_CHECK(m.get_ss().defined() == false);

  std::set<che::cchb_dssp_interval> pool;
  pool.insert(che::cchb_dssp_interval(5, 10, che::cchb_dssp('H')));
  BOOST_REQUIRE_THROW(che::molecule("storage", "id", p, che::ss(pool)), std::invalid_argument);

  che::sequence<che::cchb_dssp> v;
  v.emplace_back('H');
  m = che::molecule("storage", "id", p, che::ss(v));
  BOOST_CHECK(m.get_storage() == "storage");
  BOOST_CHECK(m.get_identifier() == "id");
  BOOST_CHECK(m.get_ps().size() == 1);
  BOOST_CHECK(m.get_ps().size() == m.get_length());
  BOOST_CHECK(m.get_ps().front().get_identifier_char() == 'A');
  BOOST_CHECK(m.get_ss().defined());
  BOOST_CHECK(m.get_ss().get_length() == m.get_length());
  BOOST_CHECK(m.get_ss().get_sequence().front().get_identifier_char() == 'H');
}

BOOST_AUTO_TEST_SUITE_END()
