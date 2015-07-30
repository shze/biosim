#include <boost/test/unit_test.hpp>

#include "che/score/cm_cc_identity.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cm_cc_identity)

BOOST_AUTO_TEST_CASE(identity_default_ctor) {
  che::score::cm_cc_identity i;

  BOOST_CHECK(i.get_identifier() == "identity");
  BOOST_CHECK(i.compare(che::cc('A'), che::cc('A')) == 1);
  BOOST_CHECK(i.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(i.compare(che::cc('A'), che::cc('R')) == -1);
}

BOOST_AUTO_TEST_CASE(identity_min_unknown) {
  che::score::cm_cc_identity i;
  BOOST_CHECK_CLOSE(i.get_min_score(), -1.0, 1e-3);
  BOOST_CHECK_CLOSE(i.get_unknown_score(), -0.9, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
