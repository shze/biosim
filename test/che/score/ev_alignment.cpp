#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "che/score/ev_alignment.h" // header to test
#include "che/io/file_fasta.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_ev_alignment)

BOOST_AUTO_TEST_CASE(ev_alignment_evaluate) {
  che::score::ev_alignment score;
  BOOST_CHECK(score.evaluate(che::alignment(che::structure())) == 0);

  che::structure s(che::io::file_fasta::read("../test/data/1a00_a.fasta").get_structure("A"));
  che::ps seq(s.get_ps());
  che::ps seq_shortened;
  seq_shortened.insert(seq_shortened.end(), seq.begin(), seq.begin() + 6);
  che::structure struct_shortened(s.get_identifier(), s.get_storage(), seq_shortened);
  BOOST_CHECK(score.evaluate(che::alignment(struct_shortened)) == 0);

  std::vector<che::structure> s_vec2(2, struct_shortened);
  che::alignment a2(s_vec2, {0, 0}, {struct_shortened.get_length(), struct_shortened.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a2), 33.5659, 1e-3);

  std::vector<che::structure> s_vec3(3, struct_shortened);
  che::alignment a3(s_vec3, {0, 0, 0}, std::vector<size_t>(3, struct_shortened.get_length()));
  BOOST_CHECK_CLOSE(score.evaluate(a3), 33.5659, 1e-3);

  che::ps seq_profile;
  seq_profile.emplace_back(che::cc("ASX"));
  che::structure struct_profile(s.get_storage(), s.get_identifier(), seq_profile);
  std::vector<che::structure> profile_vec(2, struct_profile);
  che::alignment a(profile_vec, {0, 0}, {struct_profile.get_length(), struct_profile.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 3.49272, 1e-3);

  che::ps seq_non_profile;
  seq_non_profile.insert(seq_non_profile.end(), seq.begin(), seq.begin() + 1);
  che::structure struct_non_profile(s.get_storage(), s.get_identifier(), seq_non_profile);
  std::vector<che::structure> sequence_and_profile{struct_profile, struct_non_profile};
  a = che::alignment(sequence_and_profile, {0, 0}, {struct_profile.get_length(), struct_non_profile.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), -2.60469, 1e-3);

  che::ps seq_gap;
  seq_gap.emplace_back(che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
  che::structure struct_gap(s.get_storage(), s.get_identifier(), seq_gap);
  std::vector<che::structure> sequence_and_gap{struct_non_profile, struct_gap};
  a = che::alignment(sequence_and_gap, {0, 0}, {struct_non_profile.get_length(), struct_gap.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), -4.2143, 1e-3);

  che::ps seq_n;
  seq_n.emplace_back(che::cc("ASN"));
  che::structure struct_n(s.get_storage(), s.get_identifier(), seq_n);
  std::vector<che::structure> structs{struct_profile, struct_n};
  a = che::alignment(structs, {0, 0}, {struct_profile.get_length(), struct_n.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 3.46247, 1e-3);

  structs = std::vector<che::structure>(2, struct_n);
  a = che::alignment(structs, {0, 0}, {struct_n.get_length(), struct_n.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 5.65324, 1e-3);

  che::ps seq_x;
  seq_x.emplace_back(che::cc('X'));
  che::structure struct_x(s.get_storage(), s.get_identifier(), seq_x);
  structs = std::vector<che::structure>({struct_n, struct_x});
  a = che::alignment(structs, {0, 0}, {struct_n.get_length(), struct_x.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), -0.991768, 1e-3);

  structs = std::vector<che::structure>(2, struct_x);
  a = che::alignment(structs, {0, 0}, {struct_x.get_length(), struct_x.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), -1.10221, 1e-3);

  che::ps seq_i;
  seq_i.emplace_back(che::cc('I'));
  che::structure struct_i(s.get_storage(), s.get_identifier(), seq_i);
  structs = std::vector<che::structure>(2, struct_i);
  a = che::alignment(structs, {0, 0}, {struct_i.get_length(), struct_i.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 3.99851, 1e-3);

  che::ps seq_iv;
  seq_iv.emplace_back(che::cc(che::cc::weight_map({{"ILE", 1.0}, {"VAL", 1.0}})));
  che::structure struct_iv(s.get_storage(), s.get_identifier(), seq_iv);
  structs = std::vector<che::structure>({struct_i, struct_iv});
  a = che::alignment(structs, {0, 0}, {struct_i.get_length(), struct_iv.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 3.27278, 1e-3);

  structs = std::vector<che::structure>(2, struct_iv);
  a = che::alignment(structs, {0, 0}, {struct_iv.get_length(), struct_iv.get_length()});
  BOOST_CHECK_CLOSE(score.evaluate(a), 3.21538, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
