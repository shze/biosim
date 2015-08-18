#include <boost/test/unit_test.hpp>

#include "che/cc.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_cc)

BOOST_AUTO_TEST_CASE(cc_core_compounds) {
  BOOST_CHECK(che::cc::get_id_list().size() == 37);
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
  che::cc("UNK");
  che::cc("---");
  che::cc("PYL");
  che::cc("SEC");
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
  std::map<std::string, std::pair<size_t, size_t>> cc_molecule_vertices_and_edges = {{"ALA", {13, 12}},
                                                                                     {"CYS", {14, 13}},
                                                                                     {"ASP", {16, 15}},
                                                                                     {"GLU", {19, 18}},
                                                                                     {"PHE", {23, 23}},
                                                                                     {"GLY", {10, 9}},
                                                                                     {"HIS", {21, 21}},
                                                                                     {"ILE", {22, 21}},
                                                                                     {"LYS", {25, 24}},
                                                                                     {"LEU", {22, 21}},
                                                                                     {"MET", {20, 19}},
                                                                                     {"ASN", {17, 16}},
                                                                                     {"PRO", {17, 17}},
                                                                                     {"GLN", {20, 19}},
                                                                                     {"ARG", {27, 26}},
                                                                                     {"SER", {14, 13}},
                                                                                     {"THR", {17, 16}},
                                                                                     {"VAL", {19, 18}},
                                                                                     {"TRP", {27, 28}},
                                                                                     {"TYR", {24, 24}}};
  for(auto p : cc_molecule_vertices_and_edges) {
    che::cc compound(p.first);
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == p.second.first);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == p.second.second);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 1);
    BOOST_CHECK(weights.begin()->first == p.first);
    BOOST_CHECK(weights.begin()->second == 1.0);
  } // for

  cc_molecule_vertices_and_edges =
      std::map<std::string, std::pair<size_t, size_t>>{{"ASX", {15, 14}}, {"GLX", {18, 17}}};
  for(auto p : cc_molecule_vertices_and_edges) {
    che::cc compound(p.first);
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == p.second.first);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == p.second.second);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 2);
  } // for

  {
    che::cc compound("UNK");
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::unknown);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == true);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == 16);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == 15);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 20);
  }

  {
    che::cc compound("---");
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::gap);
    BOOST_CHECK(compound.is_gap() == true);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == 0);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == 0);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 1);
    BOOST_CHECK(weights.begin()->first == compound.get_identifier());
    BOOST_CHECK(weights.begin()->second == 1.0);
  }

  cc_molecule_vertices_and_edges =
      std::map<std::string, std::pair<size_t, size_t>>{{"PYL", {39, 39}}, {"SEC", {14, 13}}, {"MSE", {20, 19}}};
  for(auto p : cc_molecule_vertices_and_edges) {
    che::cc compound(p.first);
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::secondary);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::l_peptide_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == p.second.first);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == p.second.second);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 1);
    BOOST_CHECK(weights.begin()->first != p.first);
    BOOST_CHECK(weights.begin()->second == 1.0);
  } // for

  cc_molecule_vertices_and_edges = std::map<std::string, std::pair<size_t, size_t>>{
      {"DA", {36, 38}}, {"DC", {34, 35}}, {"DG", {37, 39}}, {"DT", {36, 37}}, {"DI", {35, 37}}};
  for(auto p : cc_molecule_vertices_and_edges) {
    che::cc compound(p.first);
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::dna_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == p.second.first);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == p.second.second);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 1);
    BOOST_CHECK(weights.begin()->first == p.first);
    BOOST_CHECK(weights.begin()->second == 1.0);
  } // for

  cc_molecule_vertices_and_edges = std::map<std::string, std::pair<size_t, size_t>>{
      {"A", {37, 39}}, {"C", {35, 36}}, {"G", {38, 40}}, {"U", {34, 35}}, {"I", {36, 38}}};
  for(auto p : cc_molecule_vertices_and_edges) {
    che::cc compound(p.first);
    BOOST_CHECK(compound.get_specificity() == che::cc::specificity_type::primary);
    BOOST_CHECK(compound.is_gap() == false);
    BOOST_CHECK(compound.is_unknown() == false);
    BOOST_CHECK(compound.get_monomer_type() == che::cc::monomer_type::rna_linking);
    BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == p.second.first);
    BOOST_CHECK(boost::num_edges(compound.get_molecule()) == p.second.second);
    BOOST_CHECK(compound.get_determined_atoms().empty());
    che::cc::weight_map weights = compound.get_weights();
    BOOST_CHECK(weights.size() == 1);
    BOOST_CHECK(weights.begin()->first == p.first);
    BOOST_CHECK(weights.begin()->second == 1.0);
  } // for

  BOOST_REQUIRE_THROW(che::cc compound("this_does_not_exist"), che::cc_data_not_found);
}

BOOST_AUTO_TEST_CASE(cc_ctor_from_char) {
  std::string all_id_chars("ABCDEFGHIJKLMNOPQRSTUVWXYZ-");

  std::map<char, std::string> id_char_to_id_map = {{'A', "ALA"},
                                                   {'B', "ASX"},
                                                   {'C', "CYS"},
                                                   {'D', "ASP"},
                                                   {'E', "GLU"},
                                                   {'F', "PHE"},
                                                   {'G', "GLY"},
                                                   {'H', "HIS"},
                                                   {'I', "ILE"},
                                                   {'K', "LYS"},
                                                   {'L', "LEU"},
                                                   {'M', "MET"},
                                                   {'N', "ASN"},
                                                   {'P', "PRO"},
                                                   {'Q', "GLN"},
                                                   {'R', "ARG"},
                                                   {'S', "SER"},
                                                   {'T', "THR"},
                                                   {'V', "VAL"},
                                                   {'W', "TRP"},
                                                   {'X', "UNK"},
                                                   {'Y', "TYR"},
                                                   {'Z', "GLX"},
                                                   {'-', "---"}};
  for(auto c : all_id_chars) {
    if(id_char_to_id_map.find(c) == id_char_to_id_map.end()) {
      BOOST_REQUIRE_THROW(che::cc compound(c), che::cc_data_not_found);
      BOOST_REQUIRE_THROW(che::cc compound(c, che::cc::monomer_type::l_peptide_linking), che::cc_data_not_found);
    } // if
    else {
      BOOST_CHECK(che::cc(c).get_identifier() == id_char_to_id_map.at(c));
      BOOST_CHECK(che::cc(c, che::cc::monomer_type::l_peptide_linking).get_identifier() == id_char_to_id_map.at(c));
    } // else
  } // for

  id_char_to_id_map = std::map<char, std::string>{{'A', "DA"}, {'C', "DC"}, {'G', "DG"}, {'I', "DI"}, {'T', "DT"}};
  for(auto c : all_id_chars) {
    if(id_char_to_id_map.find(c) == id_char_to_id_map.end()) {
      BOOST_REQUIRE_THROW(che::cc compound(c, che::cc::monomer_type::dna_linking), che::cc_data_not_found);
    } // if
    else {
      BOOST_CHECK(che::cc(c, che::cc::monomer_type::dna_linking).get_identifier() == id_char_to_id_map.at(c));
    } // else
  } // for

  id_char_to_id_map = std::map<char, std::string>{{'A', "A"}, {'C', "C"}, {'G', "G"}, {'I', "I"}, {'U', "U"}};
  for(auto c : all_id_chars) {
    if(id_char_to_id_map.find(c) == id_char_to_id_map.end()) {
      BOOST_REQUIRE_THROW(che::cc compound(c, che::cc::monomer_type::rna_linking), che::cc_data_not_found);
    } // if
    else {
      BOOST_CHECK(che::cc(c, che::cc::monomer_type::rna_linking).get_identifier() == id_char_to_id_map.at(c));
    } // else
  } // for

  id_char_to_id_map = std::map<char, std::string>{};
  for(auto c : all_id_chars) {
    if(id_char_to_id_map.find(c) == id_char_to_id_map.end()) {
      BOOST_REQUIRE_THROW(che::cc compound(c, che::cc::monomer_type::non_polymer), che::cc_data_not_found);
    } // if
    else {
      BOOST_CHECK(che::cc(c, che::cc::monomer_type::non_polymer).get_identifier() == id_char_to_id_map.at(c));
    } // else
  } // for
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
  BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == 16);
  BOOST_CHECK(boost::num_edges(compound.get_molecule()) == 15);
  BOOST_CHECK(compound.get_determined_atoms().empty());
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
  BOOST_CHECK(boost::num_vertices(compound.get_molecule()) == 0);
  BOOST_CHECK(boost::num_edges(compound.get_molecule()) == 0);
  BOOST_CHECK(compound.get_determined_atoms().empty());
  weights = compound.get_weights();
  BOOST_CHECK(weights.size() == 2);
}

BOOST_AUTO_TEST_SUITE_END()
