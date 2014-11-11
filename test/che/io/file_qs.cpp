#include <boost/test/unit_test.hpp>

#include "che/io/file_qs.h" // header to test
#include "tools/file.h" // for the exception

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_qs)

BOOST_AUTO_TEST_CASE(file_qs_read) {
  che::io::file_qs reader;
  
  che::qs q(reader.read("../test/data/1UBI_A_DSSP.fasta"));
  BOOST_CHECK(q.get_chain_id_list().size() == 2);
  BOOST_CHECK(q.has_ts("A"));
  BOOST_CHECK(q.has_ts("B"));
  BOOST_CHECK(q.has_ss("A"));
  BOOST_CHECK(q.has_ss("B"));

  q = reader.read("../test/data/T0666.psipred_ss");
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 17);

  q = reader.read("../test/data/T0666.psipred_ss2");
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 10);

  q = reader.read("../test/data/T0666.rdbProf");
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 13);

  q = reader.read("../test/data/T0666.jufo9d_ss");
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 10);

  q = reader.read("../test/data/T0666_3UX4A_fixedbcl_dssp.pdb");
  BOOST_CHECK(q.get_ts("A").get_length() == 187);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 187);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 7);

  BOOST_REQUIRE_THROW(reader.read("../test/data/3IM3-biomolecule-bcl.pool"), tools::unknown_file_format);
  
  reader.add_sse_pool_reader(q);
  q = reader.read("../test/data/3IM3-biomolecule-bcl.pool");
  BOOST_CHECK(q.get_chain_id_list().size() == 2);
  BOOST_CHECK(q.get_ts("A").get_length() == 187);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 187);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 3);
  BOOST_CHECK(q.get_ts("B").get_length() == 49);
  BOOST_CHECK(q.get_ss("B").get_sequence().size() == 49);
  BOOST_CHECK(q.get_ss("B").get_sses().size() == 3);
}

BOOST_AUTO_TEST_SUITE_END()
