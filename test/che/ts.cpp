#include <boost/test/unit_test.hpp>

#include "che/ts.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ts)

BOOST_AUTO_TEST_CASE(ts_ctor) {
  che::ts t;
  BOOST_CHECK(t.get_storage().empty());
  BOOST_CHECK(t.get_identifier().empty());
  BOOST_CHECK(t.get_ps().size() == 0);

  t = che::ts("storage", "id");
  BOOST_CHECK(t.get_storage() == "storage");
  BOOST_CHECK(t.get_identifier() == "id");
}

BOOST_AUTO_TEST_CASE(ts_ps) {
  che::ts t;
  che::ps p;
  p.emplace_back('A');
  t.set_ps(p);
  BOOST_CHECK(t.get_ps().size() == 1);
  BOOST_CHECK(t.get_ps().front().get_identifier_char() == 'A');
}

BOOST_AUTO_TEST_SUITE_END()
