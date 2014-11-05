#include <boost/test/unit_test.hpp>

#include "che/io/file_qs.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_qs)

BOOST_AUTO_TEST_CASE(file_qs_read) {
  che::io::file_qs reader;
  
  che::qs q(reader.read("../test/data/T0666.psipred_ss2"));
  BOOST_CHECK(q.get_ts("A").get_length() == 195);
  BOOST_CHECK(q.get_ss("A").get_sequence().size() == 195);
  BOOST_CHECK(q.get_ss("A").get_sses().size() == 10);
}

BOOST_AUTO_TEST_SUITE_END()
