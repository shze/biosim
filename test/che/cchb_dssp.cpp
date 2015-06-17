#include <boost/test/unit_test.hpp>

#include "che/cchb_dssp.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cchb_dssp)

BOOST_AUTO_TEST_CASE(cchb_dssp_core_symbol) {
  BOOST_CHECK(che::cchb_dssp::get_id_list().size() == 9);
  che::cchb_dssp('H');
  che::cchb_dssp('E');
  che::cchb_dssp('C');
  che::cchb_dssp('-');
  che::cchb_dssp('G');
  che::cchb_dssp('I');
  che::cchb_dssp('B');
  che::cchb_dssp('T');
  che::cchb_dssp('S');
}

BOOST_AUTO_TEST_CASE(cchb_dssp_ctor_from_id) {
  che::cchb_dssp hb_symbol('H');
  BOOST_CHECK(hb_symbol.get_identifier() == 'H');
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb_dssp::specificity_type::defined);
  BOOST_CHECK(hb_symbol.is_gap() == false);
  BOOST_CHECK(hb_symbol.is_unknown() == false);
  che::cchb_dssp::weight_map weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == 'H');
  BOOST_CHECK(weights.begin()->second == 1.0);

  BOOST_REQUIRE_THROW(che::cchb_dssp compound('@'), che::cchb_dssp_data_not_found);
}

BOOST_AUTO_TEST_CASE(cchb_dssp_ctor_from_specificity) {
  che::cchb_dssp hb_symbol(che::cchb_dssp::specificity_type::gap);
  BOOST_CHECK(hb_symbol.get_identifier() == '-');
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb_dssp::specificity_type::gap);
  BOOST_CHECK(hb_symbol.is_gap() == true);
  BOOST_CHECK(hb_symbol.is_unknown() == false);
  che::cchb_dssp::weight_map weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == '-');
  BOOST_CHECK(weights.begin()->second == 1.0);

  hb_symbol = che::cchb_dssp(che::cchb_dssp::specificity_type::unknown);
  BOOST_CHECK(hb_symbol.get_identifier() == 'C');
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb_dssp::specificity_type::unknown);
  BOOST_CHECK(hb_symbol.is_gap() == false);
  BOOST_CHECK(hb_symbol.is_unknown() == true);
  weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == 'C');
  BOOST_CHECK(weights.begin()->second == 1.0);
}

BOOST_AUTO_TEST_CASE(cchb_dssp_ctor_from_weights) {
  che::cchb_dssp::weight_map empty_weights;
  BOOST_REQUIRE_THROW(che::cchb_dssp hb_symbol('C', empty_weights), che::cchb_dssp_data_not_found);

  che::cchb_dssp::weight_map multiple_weights = {{'H', 0.3}, {'E', 0.4}, {'C', 0.3}};
  che::cchb_dssp hb_symbol('H', multiple_weights);
  BOOST_CHECK(hb_symbol.get_identifier() == 'H');
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb_dssp::specificity_type::profile);
  BOOST_CHECK(hb_symbol.is_gap() == false);
  BOOST_CHECK(hb_symbol.is_unknown() == false);
  che::cchb_dssp::weight_map weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 3);
  BOOST_CHECK(weights['H'] == 0.3);
  BOOST_CHECK(weights['E'] == 0.4);
  BOOST_CHECK(weights['C'] == 0.3);
}

BOOST_AUTO_TEST_SUITE_END()
