#include <boost/test/unit_test.hpp>

#include "che/cc.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cc)

BOOST_AUTO_TEST_CASE(cc_core_compounds) {
  BOOST_CHECK(che::cc::get_id_list().size() == 38);
  che::cc("ALA");
  che::cc("CYS");
  che::cc("ASP");
  che::cc("GLU");
  che::cc("PHE");
  che::cc("GLY");
  che::cc("HIS");
  che::cc("ILE");
  che::cc("LYS");
  che::cc("LEU");
  che::cc("MET");
  che::cc("ASN");
  che::cc("PRO");
  che::cc("GLN");
  che::cc("ARG");
  che::cc("SER");
  che::cc("THR");
  che::cc("VAL");
  che::cc("TRP");
  che::cc("TYR");
  che::cc("ASX");
  che::cc("GLX");
  che::cc("XLE");
  che::cc("UNK");
  che::cc("---");
  che::cc("PYH");
  che::cc("CSE");
  che::cc("MSE");
  che::cc("DA");
  che::cc("DC");
  che::cc("DG");
  che::cc("DT");
  che::cc("DI");
  che::cc("A");
  che::cc("C");
  che::cc("G");
  che::cc("U");
  che::cc("I");
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_id) {
  che::cc compound("ALA");
  BOOST_CHECK(compound.get_identifier() == "ALA");
  BOOST_CHECK(compound.get_identifier_char() == 'A');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::defined);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  che::cc::weight_map weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "ALA");
  BOOST_CHECK(weights.begin()->second == 1.0);

  BOOST_REQUIRE_THROW(che::cc compound("this_does_not_exist"), che::cc_data_not_found);
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_char) {
  che::cc compound('C', che::cc::monomer_type::l_peptide_linking);
  BOOST_CHECK(compound.get_identifier() == "CYS");
  BOOST_CHECK(compound.get_identifier_char() == 'C');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::defined);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  che::cc::weight_map weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "CYS");
  BOOST_CHECK(weights.begin()->second == 1.0);

  BOOST_REQUIRE_THROW(che::cc compound('B', che::cc::monomer_type::dna_linking), che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound('B', che::cc::monomer_type::rna_linking), che::cc_data_not_found);
  // maybe this should be done for every single char
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_specificity) {
  che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking);
  BOOST_CHECK(compound.get_identifier() == "---");
  BOOST_CHECK(compound.get_identifier_char() == '-');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::gap);
  BOOST_CHECK(compound.is_gap() == true);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  che::cc::weight_map weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "---");
  BOOST_CHECK(weights.begin()->second == 1.0);

  compound = che::cc(che::cc::specificity_type::unknown, che::cc::monomer_type::l_peptide_linking);
  BOOST_CHECK(compound.get_identifier() == "UNK");
  BOOST_CHECK(compound.get_identifier_char() == 'X');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::unknown);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == true);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "UNK");
  BOOST_CHECK(weights.begin()->second == 1.0);

  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::profile, che::cc::monomer_type::l_peptide_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::profile, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::unknown, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::profile, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::unknown, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::profile, che::cc::monomer_type::rna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::unknown, che::cc::monomer_type::rna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::rna_linking),
                      che::cc_data_not_found);
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_weights) {
  che::cc::weight_map empty_weights;
  BOOST_REQUIRE_THROW(che::cc compound(empty_weights), che::cc_data_not_found);

  che::cc::weight_map single_weight = {{"ASP", 1.0}};
  che::cc compound(single_weight);
  BOOST_CHECK(compound.get_identifier() == "ASP");
  BOOST_CHECK(compound.get_identifier_char() == 'D');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::profile);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  che::cc::weight_map weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 1);
  BOOST_CHECK(weights.begin()->first == "ASP");
  BOOST_CHECK(weights.begin()->second == 1.0);

  che::cc::weight_map multiple_weights = {{"GLU", 0.5}, {"PHE", 0.5}};
  compound = che::cc(multiple_weights);
  BOOST_CHECK(compound.get_identifier() == "GLU");
  BOOST_CHECK(compound.get_identifier_char() == 'E');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::profile);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
