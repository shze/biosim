#include <boost/test/unit_test.hpp>

#include "che/algo/aligner_dp.h" // header to test
#include "che/io/file_fasta.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_aligner_dp)

BOOST_AUTO_TEST_CASE(aligner_dp_align_pair) {
  che::assembly a(che::io::file_fasta::read("../test/data/P01236-shorter.fasta"));
  che::structure s(a.get_structure("A"));

  che::algo::aligner_dp aligner;
  std::list<che::scored_alignment> alignment_list(aligner.align_pair(che::alignment(s), che::alignment(s)));
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 2);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 5);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 24.2329, 1e-3);

  che::alignment two_m({s, s}, {0, 0}, {s.get_length(), s.get_length()});
  alignment_list = aligner.align_pair(che::alignment(s), two_m);
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 3);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 5);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 24.2329, 1e-3);

  alignment_list = aligner.align_pair(two_m, two_m);
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 4);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 5);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 24.2329, 1e-3);

  che::ps profile;
  profile.insert(profile.end(), 5, che::cc("ASX"));
  che::structure s2("st", "id", profile);
  alignment_list = aligner.align_pair(che::alignment(s), che::alignment(s2));
  BOOST_CHECK(alignment_list.size() == 5);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 2);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 1);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 0.169935, 1e-3);

  alignment_list = aligner.align_pair(che::alignment(s2), che::alignment(s2));
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 2);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 5);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 17.4636, 1e-3);

  a = che::assembly(che::io::file_fasta::read("../test/data/P01942-shorter.fasta"));
  s = a.get_structure("A");
  a = che::assembly(che::io::file_fasta::read("../test/data/P02088-shorter.fasta"));
  s2 = a.get_structure("A");
  alignment_list = aligner.align_pair(che::alignment(s), che::alignment(s2));
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 2);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 26);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 50.0037, 1e-3);
}

BOOST_AUTO_TEST_CASE(aligner_dp_align_multiple) {
  che::assembly a(che::io::file_fasta::read("../test/data/P01236-shorter.fasta"));
  che::structure s(a.get_structure("A"));

  che::algo::aligner_dp aligner;
  std::list<che::scored_alignment> alignment_list(
      aligner.align_multiple({che::alignment(s), che::alignment(s), che::alignment(s)}));
  BOOST_CHECK(alignment_list.size() == 1);
  BOOST_CHECK(alignment_list.front().get_alignment().get_depth() == 3);
  BOOST_CHECK(alignment_list.front().get_alignment().get_length() == 5);
  BOOST_CHECK_CLOSE(alignment_list.front().get_score(), 24.2329, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
