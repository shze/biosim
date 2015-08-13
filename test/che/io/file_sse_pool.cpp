#include <boost/test/unit_test.hpp>

#include "che/io/file_sse_pool.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_sse_pool)

BOOST_AUTO_TEST_CASE(file_sse_pool_read) {
  che::assembly a = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS.pool");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[0].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 20);

  che::assembly a2 = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS_short.pool");
  BOOST_CHECK(a2.get_chain_id_list().size() == 1);
  BOOST_CHECK(a2.get_structure("A").get_length() == 163);
  BOOST_CHECK(a2.get_structure("A").get_ps().size() == 163);
  BOOST_CHECK(a2.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a2.get_structure("A").get_ss().defined());
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sequence().size() == 163);
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sequence()[0].get_identifier_char() == 'H');
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sses().size() == 17);
  a2 = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS_short.pool", a);
  BOOST_CHECK(a2.get_chain_id_list().size() == 1);
  BOOST_CHECK(a2.get_structure("A").get_length() == 195);
  BOOST_CHECK(a2.get_structure("A").get_ps().size() == 195);
  BOOST_CHECK(a2.get_structure("A").get_ps()[0].get_identifier_char() == 'M');
  BOOST_CHECK(a2.get_structure("A").get_ss().defined());
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sequence().size() == 195);
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sequence()[0].get_identifier_char() == 'H');
  BOOST_CHECK(a2.get_structure("A").get_ss().get_sses().size() == 17);

  a = che::io::file_sse_pool::read("../test/data/3IM3-biomolecule-bcl.pool");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_structure("A").get_length() == 49);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 49);
  BOOST_CHECK(a.get_structure("A").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 49);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 3);
  BOOST_CHECK(a.get_structure("B").get_length() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("B").get_ss().defined());
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("B").get_ss().get_sses().size() == 3);

  a = che::io::file_sse_pool::read("../test/data/T0666_3UX4A_fixedbcl_dssp.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 187);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 187);
  BOOST_CHECK(a.get_structure("A").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 187);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 7);

  a = che::io::file_sse_pool::read("../test/data/3IM3.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.get_structure("A").get_length() == 60);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 60);
  BOOST_CHECK(a.get_structure("A").get_ps()[13].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 60);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[13].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 3);

  a = che::io::file_sse_pool::read("../test/data/3IM3-biomolecule-bcl.pdb");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.get_structure("A").get_length() == 49);
  BOOST_CHECK(a.get_structure("A").get_ps().size() == 49);
  BOOST_CHECK(a.get_structure("A").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("A").get_ss().defined());
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence().size() == 49);
  BOOST_CHECK(a.get_structure("A").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("A").get_ss().get_sses().size() == 3);
  BOOST_CHECK(a.get_structure("B").get_length() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ps()[3].get_identifier_char() == 'X');
  BOOST_CHECK(a.get_structure("B").get_ss().defined());
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence().size() == 49);
  BOOST_CHECK(a.get_structure("B").get_ss().get_sequence()[3].get_identifier_char() == 'H');
  BOOST_CHECK(a.get_structure("B").get_ss().get_sses().size() == 3);

  BOOST_CHECK(che::io::file_sse_pool::read("../test/data/1a00_a.fasta").get_chain_id_list().size() == 0);
}

BOOST_AUTO_TEST_SUITE_END()
