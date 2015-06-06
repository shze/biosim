#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/score/cm_cc_blosum.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cm_cc_blosum)

BOOST_AUTO_TEST_CASE(blosum_default_ctor) {
  che::score::cm_cc_blosum b62;

  BOOST_CHECK(b62.get_identifier() == "blosum62");
  che::score::cm_cc_blosum::d_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_use_dbl_bitscore() == false);
  BOOST_CHECK(b62.get_bitscore_fraction() == 0.5);
  che::score::cm_cc_blosum::d_triangular_matrix s(b62.get_dbl_bitscore_matrix());
  BOOST_CHECK_CLOSE(s(0, 0), 3.92912, 1e-3);
  BOOST_CHECK_CLOSE(s(1, 0), -1.413500, 1e-3);
  BOOST_REQUIRE_THROW(s(0, 1), std::out_of_range);
  che::score::cm_cc_blosum::i_triangular_matrix i(b62.get_int_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 4);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 4);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);

  che::score::cm_cc_blosum b62_double(true);

  BOOST_CHECK(b62_double.get_identifier() == "blosum62_dbl");
  f = b62_double.get_frequency_matrix();
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  BOOST_CHECK(b62_double.get_use_dbl_bitscore() == true);
  BOOST_CHECK(b62_double.get_bitscore_fraction() == 0.5);
  s = b62_double.get_dbl_bitscore_matrix();
  BOOST_CHECK_CLOSE(s(0, 0), 3.92912, 1e-3);
  BOOST_CHECK_CLOSE(s(1, 0), -1.413500, 1e-3);
  BOOST_REQUIRE_THROW(s(0, 1), std::out_of_range);
  i = b62_double.get_int_bitscore_matrix();
  BOOST_CHECK(i(0, 0) == 4);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('A'), che::cc('A')), 3.92912, 1e-3);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('R'), che::cc('A')), -1.413500, 1e-3);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('A'), che::cc('R')), -1.413500, 1e-3);
}

BOOST_AUTO_TEST_CASE(blosum_frequency_ctor) {
  che::score::cm_cc_blosum b62_default;
  che::score::cm_cc_blosum b62("blosum62copy", b62_default.get_frequency_matrix());

  BOOST_CHECK(b62.get_identifier() == "blosum62copy");
  che::score::cm_cc_blosum::d_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_use_dbl_bitscore() == false);
  BOOST_CHECK(b62.get_bitscore_fraction() == 1.0);
  che::score::cm_cc_blosum::d_triangular_matrix s(b62.get_dbl_bitscore_matrix());
  BOOST_CHECK_CLOSE(s(0, 0), 1.96456, 1e-3);
  BOOST_CHECK_CLOSE(s(1, 0), -0.706750, 1e-3);
  BOOST_REQUIRE_THROW(s(0, 1), std::out_of_range);
  che::score::cm_cc_blosum::i_triangular_matrix i(b62.get_int_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 2);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 2);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);

  che::score::cm_cc_blosum b62_double("blosum62copy", b62_default.get_frequency_matrix(), true);

  BOOST_CHECK(b62_double.get_identifier() == "blosum62copy_dbl");
  f = b62_double.get_frequency_matrix();
  BOOST_CHECK(f(0, 0) == 26782.81);
  BOOST_CHECK(f(1, 0) == 5848.09);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  BOOST_CHECK(b62_double.get_use_dbl_bitscore() == true);
  BOOST_CHECK(b62_double.get_bitscore_fraction() == 1.0);
  s = b62_double.get_dbl_bitscore_matrix();
  BOOST_CHECK_CLOSE(s(0, 0), 1.96456, 1e-3);
  BOOST_CHECK_CLOSE(s(1, 0), -0.706750, 1e-3);
  BOOST_REQUIRE_THROW(s(0, 1), std::out_of_range);
  i = b62_double.get_int_bitscore_matrix();
  BOOST_CHECK(i(0, 0) == 2);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('A'), che::cc('A')), 1.96456, 1e-3);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('R'), che::cc('A')), -0.706750, 1e-3);
  BOOST_CHECK_CLOSE(b62_double.compare(che::cc('A'), che::cc('R')), -0.706750, 1e-3);
}

BOOST_AUTO_TEST_CASE(blosum_bitscore_ctor) {
  che::score::cm_cc_blosum b62_default;
  che::score::cm_cc_blosum b62("blosum62copy", b62_default.get_int_bitscore_matrix(),
                               b62_default.get_bitscore_fraction());

  BOOST_CHECK(b62.get_identifier() == "blosum62copy");
  che::score::cm_cc_blosum::d_triangular_matrix f(b62.get_frequency_matrix());
  BOOST_CHECK(f(0, 0) == 0.0);
  BOOST_CHECK(f(1, 0) == 0.0);
  BOOST_REQUIRE_THROW(f(0, 1), std::out_of_range);
  BOOST_CHECK(b62.get_use_dbl_bitscore() == false);
  BOOST_CHECK(b62.get_bitscore_fraction() == 0.5);
  che::score::cm_cc_blosum::d_triangular_matrix s(b62.get_dbl_bitscore_matrix());
  BOOST_CHECK_CLOSE(s(0, 0), 0.0, 1e-3);
  BOOST_CHECK_CLOSE(s(1, 0), 0.0, 1e-3);
  BOOST_REQUIRE_THROW(s(0, 1), std::out_of_range);
  che::score::cm_cc_blosum::i_triangular_matrix i(b62.get_int_bitscore_matrix());
  BOOST_CHECK(i(0, 0) == 4);
  BOOST_CHECK(i(1, 0) == -1);
  BOOST_REQUIRE_THROW(i(0, 1), std::out_of_range);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('A')) == 4);
  BOOST_CHECK(b62.compare(che::cc('R'), che::cc('A')) == -1);
  BOOST_CHECK(b62.compare(che::cc('A'), che::cc('R')) == -1);
}

BOOST_AUTO_TEST_SUITE_END()
