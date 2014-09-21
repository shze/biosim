#include <boost/test/unit_test.hpp>

#include "che/cchb.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cchb)

BOOST_AUTO_TEST_CASE(cchb_core_symbol) {
  BOOST_CHECK(che::cchb::get_id_list().size() == 4);
  che::cchb::cchb("H");
  che::cchb::cchb("E");
  che::cchb::cchb("C");
  che::cchb::cchb("-");
}

BOOST_AUTO_TEST_CASE(cchb_ctor_from_id) {
  che::cchb hb_symbol("H");
  BOOST_CHECK(hb_symbol.get_identifier() == "H");
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb::specificity_defined);
  BOOST_CHECK(hb_symbol.is_gap() == false);
  BOOST_CHECK(hb_symbol.is_unknown() == false);
  che::cchb::weight_map weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "H");
  BOOST_CHECK(weights.begin()->second == 1.0);

  BOOST_REQUIRE_THROW(che::cchb compound("this_does_not_exist"), che::cchb_data_not_found);
}

BOOST_AUTO_TEST_CASE(cchb_ctor_from_specificity) {
  che::cchb hb_symbol(che::cchb::specificity_gap);
  BOOST_CHECK(hb_symbol.get_identifier() == "-");
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb::specificity_gap);
  BOOST_CHECK(hb_symbol.is_gap() == true);
  BOOST_CHECK(hb_symbol.is_unknown() == false);
  che::cchb::weight_map weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "-");
  BOOST_CHECK(weights.begin()->second == 1.0);

  hb_symbol = che::cchb(che::cchb::specificity_unknown);
  BOOST_CHECK(hb_symbol.get_identifier() == "C");
  BOOST_CHECK(hb_symbol.get_specificity() == che::cchb::specificity_unknown);
  BOOST_CHECK(hb_symbol.is_gap() == false);
  BOOST_CHECK(hb_symbol.is_unknown() == true);
  weights = hb_symbol.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "C");
  BOOST_CHECK(weights.begin()->second == 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
