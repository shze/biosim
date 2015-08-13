#include <boost/test/unit_test.hpp>

#include "che/io/file_assembly.h" // header to test
#include "tools/file.h" // for the exception
#include "tools/log.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_assembly)

BOOST_AUTO_TEST_CASE(file_assembly_read) {
  che::io::file_assembly reader;

  che::assembly a(reader.read("../test/data/1UBI_A_DSSP.fasta"));
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.has_structure("A"));
  BOOST_CHECK(a.get_structure("A").get_length() == 76);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[7].get_identifier_char() == 'T');
  BOOST_CHECK(a.has_structure("B"));
  BOOST_CHECK(a.get_structure("B").get_length() == 76);
  BOOST_CHECK(a.get_structure("B").get_ps()[0].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("B").get_ss().defined());
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence()[7].get_identifier_char() == 'C');

  a = reader.read("../test/data/T0666.psipred_ss");
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 17);

  a = reader.read("../test/data/T0666.psipred_ss2");
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 10);

  a = reader.read("../test/data/T0666.rdbProf");
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'E');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 13);

  a = reader.read("../test/data/T0666.jufo9d_ss");
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 10);

  a = reader.read("../test/data/T0666_3UX4A_fixedbcl_dssp.pdb");
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 7);

  BOOST_REQUIRE_THROW(reader.read("../test/data/3IM3-biomolecule-bcl.pool"), tools::unknown_file_format);

  reader.add_sse_pool_reader(a);
  a = reader.read("../test/data/3IM3-biomolecule-bcl.pool");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 3);
  BOOST_CHECK(a.get_structure("B").get_length() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps()[14].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("B").get_ss().defined());
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence()[14].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("B").get_ss().get_sses().size() == 3);
}

BOOST_AUTO_TEST_SUITE_END()
