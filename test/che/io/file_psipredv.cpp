#include <boost/test/unit_test.hpp>

#include "che/io/file_psipredv.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_psipredv)

BOOST_AUTO_TEST_CASE(file_psipredv_read) {
  che::assembly a = che::io::file_psipredv::read("../test/data/T0666.psipred_ss2");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_first_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_first_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_first_structure("A").get_ps()[3].get_identifier_char() == 'L');
  BOOST_CHECK(a.get_first_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_first_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_first_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_first_structure("A").get_ss().get_sses().size() == 10);

  BOOST_CHECK(che::io::file_psipredv::read("../test/data/1a00_a.fasta").get_chain_id_list().size() == 0);
}

BOOST_AUTO_TEST_SUITE_END()
