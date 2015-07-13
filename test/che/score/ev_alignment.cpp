#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/score/ev_alignment.h" // header to test
#include "che/io/file_fasta.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ev_alignment)

BOOST_AUTO_TEST_CASE(ev_alignment_evaluate) {
  che::score::ev_alignment s;
  BOOST_CHECK(s.evaluate(che::alignment(che::molecule())) == 0);

  che::molecule m(che::io::file_fasta::read("../test/data/1a00_a.fasta").get_molecule("A"));
  che::ps seq(m.get_ps());
  che::ps shortened_seq;
  shortened_seq.insert(shortened_seq.end(), seq.begin(), seq.begin() + 6);
  che::molecule shortened_m(m.get_identifier(), m.get_storage(), shortened_seq);
  BOOST_CHECK(s.evaluate(che::alignment(shortened_m)) == 0);

  std::vector<che::molecule> m_vec2(2, shortened_m);
  BOOST_CHECK_CLOSE(s.evaluate(che::alignment(m_vec2)), 33.5659, 1e-3);

  std::vector<che::molecule> m_vec3(3, shortened_m);
  BOOST_CHECK_CLOSE(s.evaluate(che::alignment(m_vec3)), 33.5659, 1e-3);

  che::ps profile_seq;
  profile_seq.emplace_back(che::cc("ASX"));
  che::molecule profile_m(m.get_storage(), m.get_identifier(), profile_seq);
  std::vector<che::molecule> profile_vec(2, profile_m);
  che::alignment a(profile_vec);
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.49272, 1e-3);

  che::ps non_profile_seq;
  non_profile_seq.insert(non_profile_seq.end(), seq.begin(), seq.begin() + 1);
  che::molecule non_profile_m(m.get_storage(), m.get_identifier(), non_profile_seq);
  std::vector<che::molecule> sequence_and_profile;
  sequence_and_profile.emplace_back(profile_m);
  sequence_and_profile.emplace_back(non_profile_m);
  a = che::alignment(sequence_and_profile);
  BOOST_CHECK_CLOSE(s.evaluate(a), -2.60469, 1e-3);

  che::ps gap_seq;
  gap_seq.emplace_back(che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
  che::molecule gap_m(m.get_storage(), m.get_identifier(), gap_seq);
  std::vector<che::molecule> sequence_and_gap;
  sequence_and_gap.emplace_back(non_profile_m);
  sequence_and_gap.emplace_back(gap_m);
  a = che::alignment(sequence_and_gap);
  BOOST_CHECK_CLOSE(s.evaluate(a), -4.2143, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
