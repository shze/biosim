#include <boost/test/unit_test.hpp>

#include "che/io/file_pdb.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_pdb)

BOOST_AUTO_TEST_CASE(file_pdb_read) {
  che::assembly a = che::io::file_pdb::read("../test/data/T0666_3UX4A_fixedbcl_dssp.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ps()[1].get_identifier_char() == 'L');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[0].get_identifier_char() == 'C');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[1].get_identifier_char() == 'H');

  a = che::io::file_pdb::read("../test/data/2MQW.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 163);
  BOOST_CHECK(a.get_structure("A").get_ps()[5].get_identifier_char() == 'N');
  BOOST_CHECK(a.get_structure("A").get_ps()[6].get_identifier_char() == 'N');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[5].get_identifier_char() == 'C');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[6].get_identifier_char() == 'H');

  a = che::io::file_pdb::read("../test/data/2MQW-test2a.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 163);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'N');
  BOOST_CHECK(a.get_structure("A").get_ps()[1].get_identifier_char() == 'N');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[0].get_identifier_char() == 'C');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[1].get_identifier_char() == 'H');

  a = che::io::file_pdb::read("../test/data/2MQW-test2b.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 168);
  BOOST_CHECK(a.get_structure("A").get_ps()[2].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ps()[3].get_identifier_char() == 'N');
  BOOST_CHECK(a.get_structure("A").get_ps()[4].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ps()[5].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[2].get_identifier_char() == 'C');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');

  a = che::io::file_pdb::read("../test/data/2MQW-test3.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 153);
  BOOST_CHECK(a.get_structure("A").get_ps()[150].get_identifier_char() == 'D');
  BOOST_CHECK(a.get_structure("A").get_ps()[151].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ps()[152].get_identifier_char() == 'L');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[150].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[151].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[152].get_identifier_char() == 'H');

  a = che::io::file_pdb::read("../test/data/2MQW-wrong-helix-sheet.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 0);

  a = che::io::file_pdb::read("../test/data/2MUR.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_structure("A").get_length() == 44);
  BOOST_CHECK(a.get_structure("A").get_ps()[9].get_identifier_char() == 'S');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[9].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_structure("B").get_length() == 78);
  BOOST_CHECK(a.get_structure("B").get_ps()[2].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence()[2].get_identifier_char() == 'E');
}

BOOST_AUTO_TEST_SUITE_END()
