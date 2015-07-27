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
  che::alignment a2(m_vec2, {0, 0}, {shortened_m.get_length(), shortened_m.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a2), 33.5659, 1e-3);

  std::vector<che::molecule> m_vec3(3, shortened_m);
  che::alignment a3(m_vec3, {0, 0, 0}, {shortened_m.get_length(), shortened_m.get_length(), shortened_m.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a3), 33.5659, 1e-3);

  che::ps profile_seq;
  profile_seq.emplace_back(che::cc("ASX"));
  che::molecule profile_m(m.get_storage(), m.get_identifier(), profile_seq);
  std::vector<che::molecule> profile_vec(2, profile_m);
  che::alignment a(profile_vec, {0, 0}, {profile_m.get_length(), profile_m.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.49272, 1e-3);

  che::ps non_profile_seq;
  non_profile_seq.insert(non_profile_seq.end(), seq.begin(), seq.begin() + 1);
  che::molecule non_profile_m(m.get_storage(), m.get_identifier(), non_profile_seq);
  std::vector<che::molecule> sequence_and_profile;
  sequence_and_profile.emplace_back(profile_m);
  sequence_and_profile.emplace_back(non_profile_m);
  a = che::alignment(sequence_and_profile, {0, 0}, {profile_m.get_length(), non_profile_m.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), -2.60469, 1e-3);

  che::ps gap_seq;
  gap_seq.emplace_back(che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
  che::molecule gap_m(m.get_storage(), m.get_identifier(), gap_seq);
  std::vector<che::molecule> sequence_and_gap;
  sequence_and_gap.emplace_back(non_profile_m);
  sequence_and_gap.emplace_back(gap_m);
  a = che::alignment(sequence_and_gap, {0, 0}, {non_profile_m.get_length(), gap_m.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), -4.2143, 1e-3);

  che::ps seq_n;
  seq_n.emplace_back(che::cc("ASN"));
  che::molecule mol_n(m.get_storage(), m.get_identifier(), seq_n);
  std::vector<che::molecule> mols;
  mols.emplace_back(profile_m);
  mols.emplace_back(mol_n);
  a = che::alignment(mols, {0, 0}, {profile_m.get_length(), mol_n.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.46247, 1e-3);

  mols = std::vector<che::molecule>(2, mol_n);
  a = che::alignment(mols, {0, 0}, {mol_n.get_length(), mol_n.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 5.65324, 1e-3);

  che::ps seq_x;
  seq_x.emplace_back(che::cc('X'));
  che::molecule mol_x(m.get_storage(), m.get_identifier(), seq_x);
  mols = std::vector<che::molecule>();
  mols.emplace_back(mol_n);
  mols.emplace_back(mol_x);
  a = che::alignment(mols, {0, 0}, {mol_n.get_length(), mol_x.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), -0.991768, 1e-3);

  mols = std::vector<che::molecule>(2, mol_x);
  a = che::alignment(mols, {0, 0}, {mol_x.get_length(), mol_x.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), -1.10221, 1e-3);

  che::ps seq_i;
  seq_i.emplace_back(che::cc('I'));
  che::molecule mol_i(m.get_storage(), m.get_identifier(), seq_i);
  mols = std::vector<che::molecule>(2, mol_i);
  a = che::alignment(mols, {0, 0}, {mol_i.get_length(), mol_i.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.99851, 1e-3);

  che::ps seq_iv;
  seq_iv.emplace_back(che::cc(che::cc::weight_map({{"ILE", 1.0}, {"VAL", 1.0}})));
  che::molecule mol_iv(m.get_storage(), m.get_identifier(), seq_iv);
  mols = std::vector<che::molecule>();
  mols.emplace_back(mol_i);
  mols.emplace_back(mol_iv);
  a = che::alignment(mols, {0, 0}, {mol_i.get_length(), mol_iv.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.27278, 1e-3);

  mols = std::vector<che::molecule>(2, mol_iv);
  a = che::alignment(mols, {0, 0}, {mol_iv.get_length(), mol_iv.get_length()});
  BOOST_CHECK_CLOSE(s.evaluate(a), 3.21538, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
