#include <boost/test/unit_test.hpp>

#include "che/io/file_sse_pool.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_sse_pool)

BOOST_AUTO_TEST_CASE(file_sse_pool_read) {
  che::qs q = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS.pool");
  BOOST_CHECK(q.get_chain_id_list().size() == 1);
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 20);

  che::qs q2 = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS_short.pool");
  BOOST_CHECK(q2.get_chain_id_list().size() == 1);
  BOOST_CHECK(q2.get_ts("A").get_length() == 163);
  BOOST_CHECK(q2.get_ss("A").get_sequence().size() == 163);
  BOOST_CHECK(q2.get_ss("A").get_sses().size() == 17);
  q2 = che::io::file_sse_pool::read("../test/data/T0666A.SSPredHighest_PSIPRED_JUFO_OCTOPUS_short.pool", q);
  BOOST_CHECK(q2.get_chain_id_list().size() == 1);
  BOOST_CHECK(q2.get_ts("A").get_length() == 195);
  BOOST_CHECK(q2.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q2.get_ss("A").get_sses().size() == 17);

  q = che::io::file_sse_pool::read("../test/data/T0666_3UX4A_fixedbcl_dssp.pdb");
  BOOST_CHECK(q.get_chain_id_list().size() == 1);
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 7);

  BOOST_CHECK(che::io::file_sse_pool::read("../test/data/1a00_a.fasta").get_chain_id_list().size() == 0);
}

BOOST_AUTO_TEST_SUITE_END()
