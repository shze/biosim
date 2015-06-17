#include <boost/test/unit_test.hpp>

#include "che/io/file_pssm.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_pssm)

BOOST_AUTO_TEST_CASE(file_pssm_read) {
  che::assembly a(che::io::file_pssm::read("../test/data/P01236-short.pssm_ascii"));
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_molecule("A").get_length() == 47);
  BOOST_CHECK(a.get_molecule("A").get_ps().size() == 47);
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_identifier_char() == 'L');
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_specificity() == che::cc::specificity_type::profile);
  BOOST_CHECK(a.get_molecule("A").get_ss().defined() == false);

  a = che::io::file_pssm::read("../test/data/P01236-short-old-blast.pssm_ascii");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_molecule("A").get_length() == 47);
  BOOST_CHECK(a.get_molecule("A").get_ps().size() == 47);
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_identifier_char() == 'L');
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_specificity() == che::cc::specificity_type::profile);
  BOOST_CHECK(a.get_molecule("A").get_ss().defined() == false);
}

BOOST_AUTO_TEST_SUITE_END()
