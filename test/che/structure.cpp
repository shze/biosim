#include <boost/test/unit_test.hpp>

#include "che/structure.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_structure)

BOOST_AUTO_TEST_CASE(structure_ctor) {
  che::structure s;
  BOOST_CHECK(s.get_storage().empty());
  BOOST_CHECK(s.get_identifier().empty());
  BOOST_CHECK(s.get_ps().size() == 0);
  BOOST_CHECK(s.get_ps().size() == s.get_length());
  BOOST_CHECK(s.get_ss().defined() == false);

  s = che::structure("storage", "id");
  BOOST_CHECK(s.get_storage() == "storage");
  BOOST_CHECK(s.get_identifier() == "id");
  BOOST_CHECK(s.get_length() == 0);
  BOOST_CHECK(s.get_ss().defined() == false);

  che::ps p;
  p.emplace_back('A');
  s = che::structure("storage", "id", p);
  BOOST_CHECK(s.get_storage() == "storage");
  BOOST_CHECK(s.get_identifier() == "id");
  BOOST_CHECK(s.get_ps().size() == 1);
  BOOST_CHECK(s.get_ps().size() == s.get_length());
  BOOST_CHECK(s.get_ps().front().get_identifier_char() == 'A');
  BOOST_CHECK(s.get_ss().defined() == false);

  std::set<che::cchb_dssp_interval> pool;
  pool.insert(che::cchb_dssp_interval(5, 10, che::cchb_dssp('H')));
  BOOST_REQUIRE_THROW(che::structure("storage", "id", p, che::ss(pool)), std::invalid_argument);

  che::sequence<che::cchb_dssp> v;
  v.emplace_back('H');
  s = che::structure("storage", "id", p, che::ss(v));
  BOOST_CHECK(s.get_storage() == "storage");
  BOOST_CHECK(s.get_identifier() == "id");
  BOOST_CHECK(s.get_ps().size() == 1);
  BOOST_CHECK(s.get_ps().size() == s.get_length());
  BOOST_CHECK(s.get_ps().front().get_identifier_char() == 'A');
  BOOST_CHECK(s.get_ss().defined());
  BOOST_CHECK(s.get_ss().get_length() == s.get_length());
  BOOST_CHECK(s.get_ss().get_sequence().front().get_identifier_char() == 'H');
}

BOOST_AUTO_TEST_SUITE_END()
