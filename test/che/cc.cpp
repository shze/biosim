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

  BOOST_CHECK(che::cc('A').get_identifier() == "ALA");
  BOOST_CHECK(che::cc('C').get_identifier() == "CYS");
  BOOST_CHECK(che::cc('D').get_identifier() == "ASP");
  BOOST_CHECK(che::cc('E').get_identifier() == "GLU");
  BOOST_CHECK(che::cc('F').get_identifier() == "PHE");
  BOOST_CHECK(che::cc('G').get_identifier() == "GLY");
  BOOST_CHECK(che::cc('H').get_identifier() == "HIS");
  BOOST_CHECK(che::cc('I').get_identifier() == "ILE");
  BOOST_CHECK(che::cc('K').get_identifier() == "LYS");
  BOOST_CHECK(che::cc('L').get_identifier() == "LEU");
  BOOST_CHECK(che::cc('M').get_identifier() == "MET");
  BOOST_CHECK(che::cc('N').get_identifier() == "ASN");
  BOOST_CHECK(che::cc('P').get_identifier() == "PRO");
  BOOST_CHECK(che::cc('Q').get_identifier() == "GLN");
  BOOST_CHECK(che::cc('R').get_identifier() == "ARG");
  BOOST_CHECK(che::cc('S').get_identifier() == "SER");
  BOOST_CHECK(che::cc('T').get_identifier() == "THR");
  BOOST_CHECK(che::cc('V').get_identifier() == "VAL");
  BOOST_CHECK(che::cc('W').get_identifier() == "TRP");
  BOOST_CHECK(che::cc('Y').get_identifier() == "TYR");

  BOOST_CHECK(che::cc('B').get_identifier() == "ASX");
  BOOST_CHECK(che::cc('J').get_identifier() == "XLE");
  BOOST_CHECK(che::cc('X').get_identifier() == "UNK");
  BOOST_CHECK(che::cc('Z').get_identifier() == "GLX");
  BOOST_CHECK(che::cc('-').get_identifier() == "---");

  BOOST_REQUIRE_THROW(che::cc('O'), che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc('U'), che::cc_data_not_found);
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_id) {
  che::cc compound("ALA");
  BOOST_CHECK(compound.get_identifier() == "ALA");
  BOOST_CHECK(compound.get_identifier_char() == 'A');
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
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
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
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
  BOOST_CHECK(weights.size() == 20);

  BOOST_CHECK(che::cc(che::cc::specificity_type::primary, che::cc::monomer_type::l_peptide_linking).get_identifier() ==
              "ALA");
  BOOST_CHECK(che::cc(che::cc::specificity_type::secondary, che::cc::monomer_type::l_peptide_linking)
                  .get_identifier() == "CSE");

  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::primary, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::secondary, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::unknown, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::non_polymer),
                      che::cc_data_not_found);

  BOOST_CHECK(che::cc(che::cc::specificity_type::primary, che::cc::monomer_type::dna_linking).get_identifier() == "DA");
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::secondary, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::unknown, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::gap, che::cc::monomer_type::dna_linking),
                      che::cc_data_not_found);

  BOOST_CHECK(che::cc(che::cc::specificity_type::primary, che::cc::monomer_type::rna_linking).get_identifier() == "A");
  BOOST_REQUIRE_THROW(che::cc compound(che::cc::specificity_type::secondary, che::cc::monomer_type::rna_linking),
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
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
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
  BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::secondary);
  BOOST_CHECK(compound.is_gap() == false);
  BOOST_CHECK(compound.is_unknown() == false);
  BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
  weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
