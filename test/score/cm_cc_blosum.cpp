#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "score/cm_cc_blosum.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cm_cc_blosum)

BOOST_AUTO_TEST_CASE(blosum_default_ctor) {
  score::cm_cc_blosum b62;

  BOOST_CHECK(b62.get_identifier() == "blosum62");
  score::cm_cc_blosum::dbl_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  score::cm_cc_blosum::int_triangular_matrix i(b62.get_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 4);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_bitscore_fraction() == 0.5);

  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 4);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);

  BOOST_CHECK(b62.frequency_compare(che::cc('A'), che::cc('A')) == 26782.81);
  BOOST_CHECK(b62.frequency_compare(che::cc('R'), che::cc('A')) == 5848.09);
  BOOST_CHECK(b62.frequency_compare(che::cc('A'), che::cc('R')) == 5848.09);
}

BOOST_AUTO_TEST_CASE(blosum_frequency_ctor) {
  score::cm_cc_blosum b62_default;
  score::cm_cc_blosum b62("blosum62copy", b62_default.get_frequency_matrix());

  BOOST_CHECK(b62.get_identifier() == "blosum62copy");
  score::cm_cc_blosum::dbl_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  score::cm_cc_blosum::int_triangular_matrix i(b62.get_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 2);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_bitscore_fraction() == 1.0);

  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 2);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);

  BOOST_CHECK(b62.frequency_compare(che::cc('A'), che::cc('A')) == 26782.81);
  BOOST_CHECK(b62.frequency_compare(che::cc('R'), che::cc('A')) == 5848.09);
  BOOST_CHECK(b62.frequency_compare(che::cc('A'), che::cc('R')) == 5848.09);
}

BOOST_AUTO_TEST_CASE(blosum_bitscore_ctor) {
  score::cm_cc_blosum b62_default;
  score::cm_cc_blosum b62("blosum62copy", b62_default.get_bitscore_matrix(), b62_default.get_bitscore_fraction());

  BOOST_CHECK(b62.get_identifier() == "blosum62copy");
  score::cm_cc_blosum::dbl_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 0.0);
  BOOST_CHECK(f(1, 0) == 0.0);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  score::cm_cc_blosum::int_triangular_matrix i(b62.get_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 4);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_bitscore_fraction() == 0.5);

  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 4);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);

  BOOST_REQUIRE_THROW(b62.frequency_compare(che::cc('A'), che::cc('A')), score::compare_error);
  BOOST_REQUIRE_THROW(b62.frequency_compare(che::cc('R'), che::cc('A')), score::compare_error);
  BOOST_REQUIRE_THROW(b62.frequency_compare(che::cc('A'), che::cc('R')), score::compare_error);
}

BOOST_AUTO_TEST_SUITE_END()
