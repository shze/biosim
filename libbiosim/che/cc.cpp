#include "che/cc.h"
#include "tools/log.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace biosim {
  namespace che {
    // construct cc with id; id is assumed unique; throws if not found
    cc::cc(std::string const &__id) : _impl(find(__id)), _atoms() {}
    // ctor from id_char and monomer_type, assumed unique; use const ref for args b/c they're no sink args
    cc::cc(char const &__id_char, monomer_type const &__monomer) : _impl(find(__id_char, __monomer)), _atoms() {}
    // ctor from weight set
    cc::cc(weight_map __weights) : _impl(), _atoms() {
      if(__weights.empty()) { // we don't know specificity nor monomer_type, so we cannot even set it to unknown
        throw cc_data_not_found("weight_map empty");
      } // if

      std::shared_ptr<cc::data const> first_compound(find(__weights.begin()->first)); // to get the monomer_type
      double highest_weight(0.0);
      std::shared_ptr<cc::data const> highest_weighted_compound(
          find('X', first_compound->_monomer)); // initialize with unknown
      for(auto const &p : __weights) {
        std::shared_ptr<cc::data const> compound(find(p.first));
        if(p.second > highest_weight) {
          highest_weight = p.second;
          highest_weighted_compound = compound;
        } // if
      } // for

      // if only one cc is in the weight_map, use it and do not construct a profile cc
      _impl = __weights.size() == 1
                  ? first_compound
                  : std::make_shared<data>(highest_weighted_compound->_id, highest_weighted_compound->_id_char,
                                           cc::specificity_type::secondary, first_compound->_monomer, __weights);
    } // ctor

    // get identifier
    std::string cc::get_identifier() const { return _impl->_id; }
    // get identifier char, i.e. single letter code
    char cc::get_identifier_char() const { return _impl->_id_char; }
    // get specificity
    cc::specificity_type cc::get_specificity() const { return _impl->_specificity; }
    // return if an unknown (unspecified) compound; only for unknown, not for gap.
    bool cc::is_unknown() const { return _impl->_specificity == specificity_type::unknown; }
    // return if gap
    bool cc::is_gap() const { return _impl->_specificity == specificity_type::gap; }
    // get monomer type
    cc::monomer_type cc::get_monomer_type() const { return _impl->_monomer; }
    // get profile weights
    cc::weight_map cc::get_weights() const { return _impl->_weights; }
    // get molecule
    molecule const &cc::get_molecule() const { return _impl->_molecule; }
    // get determined atoms
    std::set<atom> const &cc::get_determined_atoms() const { return _atoms; }

    // (static) get list of all identifiers; static public interface
    std::list<std::string> cc::get_id_list() {
      std::list<std::string> ids;
      for(auto p : data_enum::get_instances()) {
        ids.push_back(p.first);
      } // for
      return ids;
    } // get_id_list()
    // (static) create a string of identifier chars; static public interface
    std::string cc::get_identifier_char_string() {
      std::string id_chars;
      for(auto p : data_enum::get_instances()) {
        id_chars.push_back(p.second.get_object()->_id_char);
      } // for
      return id_chars;
    } // get_identifier_char_string

    // ctor for primary unique-single-letter AA; __monomer_type has a default
    cc::data::data(std::string __id, char __id_char, molecule __molecule, cc::monomer_type __monomer)
        : _id(__id),
          _id_char(__id_char),
          _specificity(specificity_type::primary),
          _monomer(__monomer),
          _weights(),
          _molecule(__molecule) {
      _weights = {{_id, 1.0}};
    } // data ctor
    // ctor for profile unique-single-letter AA; ctor does not check if base ids have same specificity and
    // monomer values; it assumes they are consistent and uses the values of the first base id.
    cc::data::data(std::string __id, char __id_char, std::string __base_id1, std::string __base_id2,
                   molecule __molecule)
        : _id(__id), _id_char(__id_char), _specificity(), _monomer(), _weights(), _molecule(__molecule) {
      std::shared_ptr<data const> base_ptr1(find(__base_id1));
      std::shared_ptr<data const> base_ptr2(find(__base_id2));
      _specificity = specificity_type::primary;
      _monomer = base_ptr1->_monomer;
      _weights = {{__base_id1, 0.5}, {__base_id2, 0.5}};
    } // data ctor
    // ctor for derived AA
    cc::data::data(std::string __id, std::string __base_id, molecule __molecule)
        : _id(__id), _id_char(), _specificity(), _monomer(), _weights(), _molecule(__molecule) {
      std::shared_ptr<data const> base_ptr(find(__base_id)); // find base id, use it to set the following values
      _id_char = base_ptr->_id_char;
      _specificity = specificity_type::secondary;
      _monomer = base_ptr->_monomer;
      _weights = {{__base_id, 1.0}};
    } // data ctor
    // ctor for profile, unknown, and gap AA
    cc::data::data(std::string __id, char __id_char, cc::specificity_type __specificity, cc::monomer_type __monomer,
                   weight_map __weights, molecule __molecule)
        : _id(__id),
          _id_char(__id_char),
          _specificity(__specificity),
          _monomer(__monomer),
          _weights(normalize_weights(__weights)),
          _molecule(__molecule) {}

    // (static) convert an atom-bond definition into a molecule; the atom-bond definition must be in CCD bond format:
    // comp_id atom_id_1 atom_id_2 value_order pdbx_aromatic_flag pdbx_stereo_config pdbx_ordinal
    molecule cc::data::to_molecule(std::list<std::string> __atom_bond_definition) {
      std::map<std::string, che::molecule::vertex_descriptor> molecule_map; // map atom to vertex, clean duplicate atoms
      std::set<std::pair<std::string, std::string>> bond_set; // to remove duplicate bonds
      che::molecule m; // final molecule
      for(auto bond_string : __atom_bond_definition) {
        std::vector<std::string> bond_string_parts;
        boost::algorithm::split(bond_string_parts, bond_string, boost::is_any_of(" \t"), boost::token_compress_on);
        std::string first_atom(bond_string_parts[1]);
        std::string second_atom(bond_string_parts[2]);
        std::string bond_value(bond_string_parts[3]);

        if(molecule_map.find(first_atom) == molecule_map.end()) { // if atom is not found
          molecule_map[first_atom] = m.add_vertex(che::atom(first_atom)); // add atom to molecule, store index in map
        } // if
        if(molecule_map.find(second_atom) == molecule_map.end()) { // if atom is not found
          molecule_map[second_atom] = m.add_vertex(che::atom(second_atom)); // add atom to molecule, store index in map
        } // if

        if(bond_set.find(std::pair<std::string, std::string>(first_atom, second_atom)) != bond_set.end()) {
          DEBUG << "Found bond id '(" << first_atom << ", " << second_atom << ")' multiple times, ignoring.";
          continue;
        } // if
        bond_set.insert(std::pair<std::string, std::string>(first_atom, second_atom)); // add bond to set
        m.add_edge(molecule_map[first_atom], molecule_map[second_atom], che::bond(bond_value)); // and graph
      } // for

      return m;
    } // to_molecule()
    // (static) normalize the given weight map (so that the sum of all weights is 1)
    cc::weight_map cc::data::normalize_weights(weight_map __weights) {
      double sum(0.0);
      for(auto const &p : __weights) {
        sum += p.second;
      } // for
      for(auto &p : __weights) {
        p.second /= sum;
      } // for
      return __weights;
    } // normalize_weights()

    // (static) find compound with id; id is assumed unique; throws if no object was found
    std::shared_ptr<cc::data const> cc::find(std::string const &__id) {
      try {
        return data_enum::get_instances().at(__id).get_object();
      } catch(std::out_of_range &e) {
        throw cc_data_not_found("map out of range, no data found with id=" + __id);
      }
    } // find()
    // (static) find compound with id_char and monomer_type; combination is assumed unique
    std::shared_ptr<cc::data const> cc::find(char const &__id_char, monomer_type const &__monomer) {
      auto itr = std::find_if(data_enum::get_instances().begin(), data_enum::get_instances().end(),
                              [&](std::pair<std::string, tools::enumerate<std::shared_ptr<data const>>> p) {
        // include all cc with the correct id_char && include all cc with the correct monomer_type && include all cc
        // that are not secondary (i.e. that have a unique single letter)
        return p.second.get_object()->_id_char == __id_char && p.second.get_object()->_monomer == __monomer &&
               p.second.get_object()->_specificity != specificity_type::secondary;
      });

      if(itr == data_enum::get_instances().end()) {
        throw cc_data_not_found(std::string("no data found with id_char=") + __id_char + " and monomer_type=" +
                                std::to_string((int)__monomer));
      } // if

      return itr->second.get_object();
    } // find()
    // (static) constructs data objects of the enumerated instance set, i.e. all single letter AA, see:
    // http://www.ebi.ac.uk/pdbe-srv/pdbechem/
    // http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bia01as00/full
    // http://www.chem.qmul.ac.uk/iupac/AminoAcid/
    // http://www.rcsb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/sequence.html
    // http://proline.bic.nus.edu.sg/~asif/AminoAcids.ppt
    bool cc::initialize() {
      // 20 primary single-letter AA
      data_enum::add(std::make_shared<data>(
          "ALA", 'A',
          data::to_molecule({"ALA N   CA  SING N N 1", "ALA N   H   SING N N 2", "ALA N   H2  SING N N 3",
                             "ALA CA  C   SING N N 4", "ALA CA  CB  SING N N 5", "ALA CA  HA  SING N N 6",
                             "ALA C   O   DOUB N N 7", "ALA C   OXT SING N N 8", "ALA CB  HB1 SING N N 9",
                             "ALA CB  HB2 SING N N 10", "ALA CB  HB3 SING N N 11", "ALA OXT HXT SING N N 12"})));
      data_enum::add(std::make_shared<data>(
          "CYS", 'C', //
          data::to_molecule({"CYS N   CA  SING N N 1", "CYS N   H   SING N N 2", "CYS N   H2  SING N N 3",
                             "CYS CA  C   SING N N 4", "CYS CA  CB  SING N N 5", "CYS CA  HA  SING N N 6",
                             "CYS C   O   DOUB N N 7", "CYS C   OXT SING N N 8", "CYS CB  SG  SING N N 9",
                             "CYS CB  HB2 SING N N 10", "CYS CB  HB3 SING N N 11", "CYS SG  HG  SING N N 12",
                             "CYS OXT HXT SING N N 13"})));
      data_enum::add(std::make_shared<data>(
          "ASP", 'D',
          data::to_molecule({"ASP N   CA  SING N N 1", "ASP N   H   SING N N 2", "ASP N   H2  SING N N 3",
                             "ASP CA  C   SING N N 4", "ASP CA  CB  SING N N 5", "ASP CA  HA  SING N N 6",
                             "ASP C   O   DOUB N N 7", "ASP C   OXT SING N N 8", "ASP CB  CG  SING N N 9",
                             "ASP CB  HB2 SING N N 10", "ASP CB  HB3 SING N N 11", "ASP CG  OD1 DOUB N N 12",
                             "ASP CG  OD2 SING N N 13", "ASP OD2 HD2 SING N N 14", "ASP OXT HXT SING N N 15"})));
      data_enum::add(std::make_shared<data>(
          "GLU", 'E',
          data::to_molecule({"GLU N   CA  SING N N 1", "GLU N   H   SING N N 2", "GLU N   H2  SING N N 3",
                             "GLU CA  C   SING N N 4", "GLU CA  CB  SING N N 5", "GLU CA  HA  SING N N 6",
                             "GLU C   O   DOUB N N 7", "GLU C   OXT SING N N 8", "GLU CB  CG  SING N N 9",
                             "GLU CB  HB2 SING N N 10", "GLU CB  HB3 SING N N 11", "GLU CG  CD  SING N N 12",
                             "GLU CG  HG2 SING N N 13", "GLU CG  HG3 SING N N 14", "GLU CD  OE1 DOUB N N 15",
                             "GLU CD  OE2 SING N N 16", "GLU OE2 HE2 SING N N 17", "GLU OXT HXT SING N N 18"})));
      data_enum::add(std::make_shared<data>(
          "PHE", 'F',
          data::to_molecule({"PHE N   CA  SING N N 1", "PHE N   H   SING N N 2", "PHE N   H2  SING N N 3",
                             "PHE CA  C   SING N N 4", "PHE CA  CB  SING N N 5", "PHE CA  HA  SING N N 6",
                             "PHE C   O   DOUB N N 7", "PHE C   OXT SING N N 8", "PHE CB  CG  SING N N 9",
                             "PHE CB  HB2 SING N N 10", "PHE CB  HB3 SING N N 11", "PHE CG  CD1 DOUB Y N 12",
                             "PHE CG  CD2 SING Y N 13", "PHE CD1 CE1 SING Y N 14", "PHE CD1 HD1 SING N N 15",
                             "PHE CD2 CE2 DOUB Y N 16", "PHE CD2 HD2 SING N N 17", "PHE CE1 CZ  DOUB Y N 18",
                             "PHE CE1 HE1 SING N N 19", "PHE CE2 CZ  SING Y N 20", "PHE CE2 HE2 SING N N 21",
                             "PHE CZ  HZ  SING N N 22", "PHE OXT HXT SING N N 23"})));
      data_enum::add(std::make_shared<data>(
          "GLY", 'G',
          data::to_molecule({"GLY N   CA  SING N N 1", "GLY N   H   SING N N 2", "GLY N   H2  SING N N 3",
                             "GLY CA  C   SING N N 4", "GLY CA  HA2 SING N N 5", "GLY CA  HA3 SING N N 6",
                             "GLY C   O   DOUB N N 7", "GLY C   OXT SING N N 8", "GLY OXT HXT SING N N 9"})));
      data_enum::add(std::make_shared<data>(
          "HIS", 'H',
          data::to_molecule({"HIS N   CA  SING N N 1", "HIS N   H   SING N N 2", "HIS N   H2  SING N N 3",
                             "HIS CA  C   SING N N 4", "HIS CA  CB  SING N N 5", "HIS CA  HA  SING N N 6",
                             "HIS C   O   DOUB N N 7", "HIS C   OXT SING N N 8", "HIS CB  CG  SING N N 9",
                             "HIS CB  HB2 SING N N 10", "HIS CB  HB3 SING N N 11", "HIS CG  ND1 SING Y N 12",
                             "HIS CG  CD2 DOUB Y N 13", "HIS ND1 CE1 DOUB Y N 14", "HIS ND1 HD1 SING N N 15",
                             "HIS CD2 NE2 SING Y N 16", "HIS CD2 HD2 SING N N 17", "HIS CE1 NE2 SING Y N 18",
                             "HIS CE1 HE1 SING N N 19", "HIS NE2 HE2 SING N N 20", "HIS OXT HXT SING N N 21"})));
      data_enum::add(std::make_shared<data>(
          "ILE", 'I',
          data::to_molecule({"ILE N   CA   SING N N 1", "ILE N   H    SING N N 2", "ILE N   H2   SING N N 3",
                             "ILE CA  C    SING N N 4", "ILE CA  CB   SING N N 5", "ILE CA  HA   SING N N 6",
                             "ILE C   O    DOUB N N 7", "ILE C   OXT  SING N N 8", "ILE CB  CG1  SING N N 9",
                             "ILE CB  CG2  SING N N 10", "ILE CB  HB   SING N N 11", "ILE CG1 CD1  SING N N 12",
                             "ILE CG1 HG12 SING N N 13", "ILE CG1 HG13 SING N N 14", "ILE CG2 HG21 SING N N 15",
                             "ILE CG2 HG22 SING N N 16", "ILE CG2 HG23 SING N N 17", "ILE CD1 HD11 SING N N 18",
                             "ILE CD1 HD12 SING N N 19", "ILE CD1 HD13 SING N N 20", "ILE OXT HXT  SING N N 21"})));
      data_enum::add(std::make_shared<data>(
          "LYS", 'K',
          data::to_molecule({"LYS N   CA  SING N N 1", "LYS N   H   SING N N 2", "LYS N   H2  SING N N 3",
                             "LYS CA  C   SING N N 4", "LYS CA  CB  SING N N 5", "LYS CA  HA  SING N N 6",
                             "LYS C   O   DOUB N N 7", "LYS C   OXT SING N N 8", "LYS CB  CG  SING N N 9",
                             "LYS CB  HB2 SING N N 10", "LYS CB  HB3 SING N N 11", "LYS CG  CD  SING N N 12",
                             "LYS CG  HG2 SING N N 13", "LYS CG  HG3 SING N N 14", "LYS CD  CE  SING N N 15",
                             "LYS CD  HD2 SING N N 16", "LYS CD  HD3 SING N N 17", "LYS CE  NZ  SING N N 18",
                             "LYS CE  HE2 SING N N 19", "LYS CE  HE3 SING N N 20", "LYS NZ  HZ1 SING N N 21",
                             "LYS NZ  HZ2 SING N N 22", "LYS NZ  HZ3 SING N N 23", "LYS OXT HXT SING N N 24"})));
      data_enum::add(std::make_shared<data>(
          "LEU", 'L',
          data::to_molecule({"LEU N   CA   SING N N 1", "LEU N   H    SING N N 2", "LEU N   H2   SING N N 3",
                             "LEU CA  C    SING N N 4", "LEU CA  CB   SING N N 5", "LEU CA  HA   SING N N 6",
                             "LEU C   O    DOUB N N 7", "LEU C   OXT  SING N N 8", "LEU CB  CG   SING N N 9",
                             "LEU CB  HB2  SING N N 10", "LEU CB  HB3  SING N N 11", "LEU CG  CD1  SING N N 12",
                             "LEU CG  CD2  SING N N 13", "LEU CG  HG   SING N N 14", "LEU CD1 HD11 SING N N 15",
                             "LEU CD1 HD12 SING N N 16", "LEU CD1 HD13 SING N N 17", "LEU CD2 HD21 SING N N 18",
                             "LEU CD2 HD22 SING N N 19", "LEU CD2 HD23 SING N N 20", "LEU OXT HXT  SING N N 21"})));
      data_enum::add(std::make_shared<data>(
          "MET", 'M',
          data::to_molecule({"MET N   CA  SING N N 1", "MET N   H   SING N N 2", "MET N   H2  SING N N 3",
                             "MET CA  C   SING N N 4", "MET CA  CB  SING N N 5", "MET CA  HA  SING N N 6",
                             "MET C   O   DOUB N N 7", "MET C   OXT SING N N 8", "MET CB  CG  SING N N 9",
                             "MET CB  HB2 SING N N 10", "MET CB  HB3 SING N N 11", "MET CG  SD  SING N N 12",
                             "MET CG  HG2 SING N N 13", "MET CG  HG3 SING N N 14", "MET SD  CE  SING N N 15",
                             "MET CE  HE1 SING N N 16", "MET CE  HE2 SING N N 17", "MET CE  HE3 SING N N 18",
                             "MET OXT HXT SING N N 19"})));
      data_enum::add(std::make_shared<data>(
          "ASN", 'N',
          data::to_molecule({"ASN N   CA   SING N N 1", "ASN N   H    SING N N 2", "ASN N   H2   SING N N 3",
                             "ASN CA  C    SING N N 4", "ASN CA  CB   SING N N 5", "ASN CA  HA   SING N N 6",
                             "ASN C   O    DOUB N N 7", "ASN C   OXT  SING N N 8", "ASN CB  CG   SING N N 9",
                             "ASN CB  HB2  SING N N 10", "ASN CB  HB3  SING N N 11", "ASN CG  OD1  DOUB N N 12",
                             "ASN CG  ND2  SING N N 13", "ASN ND2 HD21 SING N N 14", "ASN ND2 HD22 SING N N 15",
                             "ASN OXT HXT  SING N N 16"})));
      data_enum::add(std::make_shared<data>(
          "PRO", 'P',
          data::to_molecule({"PRO N   CA  SING N N 1", "PRO N   CD  SING N N 2", "PRO N   H   SING N N 3",
                             "PRO CA  C   SING N N 4", "PRO CA  CB  SING N N 5", "PRO CA  HA  SING N N 6",
                             "PRO C   O   DOUB N N 7", "PRO C   OXT SING N N 8", "PRO CB  CG  SING N N 9",
                             "PRO CB  HB2 SING N N 10", "PRO CB  HB3 SING N N 11", "PRO CG  CD  SING N N 12",
                             "PRO CG  HG2 SING N N 13", "PRO CG  HG3 SING N N 14", "PRO CD  HD2 SING N N 15",
                             "PRO CD  HD3 SING N N 16", "PRO OXT HXT SING N N 17"})));
      data_enum::add(std::make_shared<data>(
          "GLN", 'Q',
          data::to_molecule({"GLN N   CA   SING N N 1", "GLN N   H    SING N N 2", "GLN N   H2   SING N N 3",
                             "GLN CA  C    SING N N 4", "GLN CA  CB   SING N N 5", "GLN CA  HA   SING N N 6",
                             "GLN C   O    DOUB N N 7", "GLN C   OXT  SING N N 8", "GLN CB  CG   SING N N 9",
                             "GLN CB  HB2  SING N N 10", "GLN CB  HB3  SING N N 11", "GLN CG  CD   SING N N 12",
                             "GLN CG  HG2  SING N N 13", "GLN CG  HG3  SING N N 14", "GLN CD  OE1  DOUB N N 15",
                             "GLN CD  NE2  SING N N 16", "GLN NE2 HE21 SING N N 17", "GLN NE2 HE22 SING N N 18",
                             "GLN OXT HXT  SING N N 19"})));
      data_enum::add(std::make_shared<data>(
          "ARG", 'R',
          data::to_molecule({"ARG N   CA   SING N N 1", "ARG N   H    SING N N 2", "ARG N   H2   SING N N 3",
                             "ARG CA  C    SING N N 4", "ARG CA  CB   SING N N 5", "ARG CA  HA   SING N N 6",
                             "ARG C   O    DOUB N N 7", "ARG C   OXT  SING N N 8", "ARG CB  CG   SING N N 9",
                             "ARG CB  HB2  SING N N 10", "ARG CB  HB3  SING N N 11", "ARG CG  CD   SING N N 12",
                             "ARG CG  HG2  SING N N 13", "ARG CG  HG3  SING N N 14", "ARG CD  NE   SING N N 15",
                             "ARG CD  HD2  SING N N 16", "ARG CD  HD3  SING N N 17", "ARG NE  CZ   SING N N 18",
                             "ARG NE  HE   SING N N 19", "ARG CZ  NH1  SING N N 20", "ARG CZ  NH2  DOUB N N 21",
                             "ARG NH1 HH11 SING N N 22", "ARG NH1 HH12 SING N N 23", "ARG NH2 HH21 SING N N 24",
                             "ARG NH2 HH22 SING N N 25", "ARG OXT HXT  SING N N 26"})));
      data_enum::add(std::make_shared<data>(
          "SER", 'S', //
          data::to_molecule({"SER N   CA  SING N N 1", "SER N   H   SING N N 2", "SER N   H2  SING N N 3",
                             "SER CA  C   SING N N 4", "SER CA  CB  SING N N 5", "SER CA  HA  SING N N 6",
                             "SER C   O   DOUB N N 7", "SER C   OXT SING N N 8", "SER CB  OG  SING N N 9",
                             "SER CB  HB2 SING N N 10", "SER CB  HB3 SING N N 11", "SER OG  HG  SING N N 12",
                             "SER OXT HXT SING N N 13"})));
      data_enum::add(std::make_shared<data>(
          "THR", 'T',
          data::to_molecule({"THR N   CA   SING N N 1", "THR N   H    SING N N 2", "THR N   H2   SING N N 3",
                             "THR CA  C    SING N N 4", "THR CA  CB   SING N N 5", "THR CA  HA   SING N N 6",
                             "THR C   O    DOUB N N 7", "THR C   OXT  SING N N 8", "THR CB  OG1  SING N N 9",
                             "THR CB  CG2  SING N N 10", "THR CB  HB   SING N N 11", "THR OG1 HG1  SING N N 12",
                             "THR CG2 HG21 SING N N 13", "THR CG2 HG22 SING N N 14", "THR CG2 HG23 SING N N 15",
                             "THR OXT HXT  SING N N 16"})));
      data_enum::add(std::make_shared<data>(
          "VAL", 'V',
          data::to_molecule({"VAL N   CA   SING N N 1", "VAL N   H    SING N N 2", "VAL N   H2   SING N N 3",
                             "VAL CA  C    SING N N 4", "VAL CA  CB   SING N N 5", "VAL CA  HA   SING N N 6",
                             "VAL C   O    DOUB N N 7", "VAL C   OXT  SING N N 8", "VAL CB  CG1  SING N N 9",
                             "VAL CB  CG2  SING N N 10", "VAL CB  HB   SING N N 11", "VAL CG1 HG11 SING N N 12",
                             "VAL CG1 HG12 SING N N 13", "VAL CG1 HG13 SING N N 14", "VAL CG2 HG21 SING N N 15",
                             "VAL CG2 HG22 SING N N 16", "VAL CG2 HG23 SING N N 17", "VAL OXT HXT  SING N N 18"})));
      data_enum::add(std::make_shared<data>(
          "TRP", 'W',
          data::to_molecule({"TRP N   CA  SING N N 1", "TRP N   H   SING N N 2", "TRP N   H2  SING N N 3",
                             "TRP CA  C   SING N N 4", "TRP CA  CB  SING N N 5", "TRP CA  HA  SING N N 6",
                             "TRP C   O   DOUB N N 7", "TRP C   OXT SING N N 8", "TRP CB  CG  SING N N 9",
                             "TRP CB  HB2 SING N N 10", "TRP CB  HB3 SING N N 11", "TRP CG  CD1 DOUB Y N 12",
                             "TRP CG  CD2 SING Y N 13", "TRP CD1 NE1 SING Y N 14", "TRP CD1 HD1 SING N N 15",
                             "TRP CD2 CE2 DOUB Y N 16", "TRP CD2 CE3 SING Y N 17", "TRP NE1 CE2 SING Y N 18",
                             "TRP NE1 HE1 SING N N 19", "TRP CE2 CZ2 SING Y N 20", "TRP CE3 CZ3 DOUB Y N 21",
                             "TRP CE3 HE3 SING N N 22", "TRP CZ2 CH2 DOUB Y N 23", "TRP CZ2 HZ2 SING N N 24",
                             "TRP CZ3 CH2 SING Y N 25", "TRP CZ3 HZ3 SING N N 26", "TRP CH2 HH2 SING N N 27",
                             "TRP OXT HXT SING N N 28"})));
      data_enum::add(std::make_shared<data>(
          "TYR", 'Y',
          data::to_molecule({"TYR N   CA  SING N N 1", "TYR N   H   SING N N 2", "TYR N   H2  SING N N 3",
                             "TYR CA  C   SING N N 4", "TYR CA  CB  SING N N 5", "TYR CA  HA  SING N N 6",
                             "TYR C   O   DOUB N N 7", "TYR C   OXT SING N N 8", "TYR CB  CG  SING N N 9",
                             "TYR CB  HB2 SING N N 10", "TYR CB  HB3 SING N N 11", "TYR CG  CD1 DOUB Y N 12",
                             "TYR CG  CD2 SING Y N 13", "TYR CD1 CE1 SING Y N 14", "TYR CD1 HD1 SING N N 15",
                             "TYR CD2 CE2 DOUB Y N 16", "TYR CD2 HD2 SING N N 17", "TYR CE1 CZ  DOUB Y N 18",
                             "TYR CE1 HE1 SING N N 19", "TYR CE2 CZ  SING Y N 20", "TYR CE2 HE2 SING N N 21",
                             "TYR CZ  OH  SING N N 22", "TYR OH  HH  SING N N 23", "TYR OXT HXT SING N N 24"})));
      // 5 profile single-letter AA, well only 3
      data_enum::add(std::make_shared<data>(
          "ASX", 'B', "ASN", "ASP",
          data::to_molecule({"ASX N   CA  SING N N 1", "ASX CA  C   SING N N 2", "ASX CA  CB  SING N N 3",
                             "ASX CA  HA  SING N N 4", "ASX C   O   DOUB N N 5", "ASX C   OXT SING N N 6",
                             "ASX CB  CG  SING N N 7", "ASX CB  HB1 SING N N 8", "ASX CB  HB2 SING N N 9",
                             "ASX CG  XD1 DOUB N N 10", "ASX CG  XD2 SING N N 11", "ASX OXT HXT SING N N 12",
                             "ASX N   H   SING N N 13", "ASX N   H2  SING N N 14"})));
      data_enum::add(std::make_shared<data>(
          "GLX", 'Z', "GLN", "GLU",
          data::to_molecule({"GLX N   CA  SING N N 1", "GLX CA  C   SING N N 2", "GLX CA  CB  SING N N 3",
                             "GLX CA  HA  SING N N 4", "GLX C   O   DOUB N N 5", "GLX C   OXT SING N N 6",
                             "GLX CB  CG  SING N N 7", "GLX CB  HB1 SING N N 8", "GLX CB  HB2 SING N N 9",
                             "GLX CG  CD  SING N N 10", "GLX CG  HG1 SING N N 11", "GLX CG  HG2 SING N N 12",
                             "GLX CD  XE1 DOUB N N 13", "GLX CD  XE2 SING N N 14", "GLX N   H   SING N N 15",
                             "GLX N   H2  SING N N 16", "GLX HXT OXT SING N N 17"})));
      // not adding the following with their single letters OUJ b/c these are not or otherwise defined in the CCD
      // PYRROLYSINE single-letter: K or O; three-letter: PYH (obsoleted) or PYL (current);
      // SELENOCYSTEINE single-letter: C or U; three-letter: CSE (obsoleted) or SEC (current);
      // LEU/ILE AMBIGIOUS see next line
      // data_enum::add(std::make_shared<data>("XLE", 'J', "LEU", "ILE", molecule()));

      // unknown AA (special primary single-letter AA)
      data_enum::add(std::make_shared<data>(
          "UNK", 'X', specificity_type::unknown, monomer_type::l_peptide_linking, weight_map({{"ALA", 1.0},
                                                                                              {"CYS", 1.0},
                                                                                              {"ASP", 1.0},
                                                                                              {"GLU", 1.0},
                                                                                              {"PHE", 1.0},
                                                                                              {"GLY", 1.0},
                                                                                              {"HIS", 1.0},
                                                                                              {"ILE", 1.0},
                                                                                              {"LYS", 1.0},
                                                                                              {"LEU", 1.0},
                                                                                              {"MET", 1.0},
                                                                                              {"ASN", 1.0},
                                                                                              {"PRO", 1.0},
                                                                                              {"GLN", 1.0},
                                                                                              {"ARG", 1.0},
                                                                                              {"SER", 1.0},
                                                                                              {"THR", 1.0},
                                                                                              {"VAL", 1.0},
                                                                                              {"TRP", 1.0},
                                                                                              {"TYR", 1.0}}),
          data::to_molecule({"UNK N   CA  SING N N 1", "UNK N   H   SING N N 2", "UNK N   H2  SING N N 3",
                             "UNK CA  C   SING N N 4", "UNK CA  CB  SING N N 5", "UNK CA  HA  SING N N 6",
                             "UNK C   O   DOUB N N 7", "UNK C   OXT SING N N 8", "UNK CB  CG  SING N N 9",
                             "UNK CB  HB1 SING N N 10", "UNK CB  HB2 SING N N 11", "UNK CG  HG1 SING N N 12",
                             "UNK CG  HG2 SING N N 13", "UNK CG  HG3 SING N N 14", "UNK OXT HXT SING N N 15"})));
      // gap (special primary single-letter AA); there seem to be no three letter gap id; used one letter codes: -~.
      data_enum::add(std::make_shared<data>("---", '-', specificity_type::gap, monomer_type::l_peptide_linking,
                                            weight_map({{"---", 1.0}})));

      // add some common non-primary AAs
      // PYRROLYSINE
      data_enum::add(std::make_shared<data>(
          "PYL", "LYS",
          data::to_molecule({"PYL CB2 CG2  SING N N 1", "PYL CB2 HB12 SING N N 2", "PYL CB2 HB22 SING N N 3",
                             "PYL CB2 HB32 SING N N 4", "PYL CG2 CD2  SING N N 5", "PYL CG2 CA2  SING N N 6",
                             "PYL CG2 HG22 SING N N 7", "PYL CD2 CE2  SING N N 8", "PYL CD2 HD32 SING N N 9",
                             "PYL CD2 HD22 SING N N 10", "PYL CE2 N2   DOUB N N 11", "PYL CE2 HE22 SING N N 12",
                             "PYL N2  CA2  SING N N 13", "PYL CA2 C2   SING N N 14", "PYL CA2 HA2  SING N N 15",
                             "PYL C2  NZ   SING N N 16", "PYL C2  O2   DOUB N N 17", "PYL NZ  CE   SING N N 18",
                             "PYL NZ  HZ   SING N N 19", "PYL CE  CD   SING N N 20", "PYL CE  HE3  SING N N 21",
                             "PYL CE  HE2  SING N N 22", "PYL CD  CG   SING N N 23", "PYL CD  HD3  SING N N 24",
                             "PYL CD  HD2  SING N N 25", "PYL CG  CB   SING N N 26", "PYL CG  HG3  SING N N 27",
                             "PYL CG  HG2  SING N N 28", "PYL CB  CA   SING N N 29", "PYL CB  HB3  SING N N 30",
                             "PYL CB  HB2  SING N N 31", "PYL CA  C    SING N N 32", "PYL CA  HA   SING N N 33",
                             "PYL CA  N    SING N N 34", "PYL C   OXT  SING N N 35", "PYL C   O    DOUB N N 36",
                             "PYL OXT HXT  SING N N 37", "PYL N   H    SING N N 38", "PYL N   H2   SING N N 39"})));
      // SELENOCYSTEINE
      data_enum::add(std::make_shared<data>(
          "SEC", "CYS", //
          data::to_molecule({"SEC N   CA  SING N N 1", "SEC N   H   SING N N 2", "SEC N   H2  SING N N 3",
                             "SEC CA  CB  SING N N 4", "SEC CA  C   SING N N 5", "SEC CA  HA  SING N N 6",
                             "SEC CB  SE  SING N N 7", "SEC CB  HB2 SING N N 8", "SEC CB  HB3 SING N N 9",
                             "SEC SE  HE  SING N N 10", "SEC C   O   DOUB N N 11", "SEC C   OXT SING N N 12",
                             "SEC OXT HXT SING N N 13"})));
      // SELENOMETHIONINE
      data_enum::add(std::make_shared<data>(
          "MSE", "MET",
          data::to_molecule({"MSE N   CA  SING N N 1", "MSE N   H   SING N N 2", "MSE N   HN2 SING N N 3",
                             "MSE CA  C   SING N N 4", "MSE CA  CB  SING N N 5", "MSE CA  HA  SING N N 6",
                             "MSE C   O   DOUB N N 7", "MSE C   OXT SING N N 8", "MSE OXT HXT SING N N 9",
                             "MSE CB  CG  SING N N 10", "MSE CB  HB2 SING N N 11", "MSE CB  HB3 SING N N 12",
                             "MSE CG  SE  SING N N 13", "MSE CG  HG2 SING N N 14", "MSE CG  HG3 SING N N 15",
                             "MSE SE  CE  SING N N 16", "MSE CE  HE1 SING N N 17", "MSE CE  HE2 SING N N 18",
                             "MSE CE  HE3 SING N N 19"})));

      // add deoxyribonucleotides
      data_enum::add(std::make_shared<data>(
          "DA", 'A', //
          data::to_molecule({"DA  OP3 P    SING N N 1", "DA  OP3 HOP3 SING N N 2", "DA  P   OP1  DOUB N N 3",
                             "DA  P   OP2  SING N N 4", "DA  P   O5'  SING N N 5", "DA  OP2 HOP2 SING N N 6",
                             "DA  O5' C5'  SING N N 7", "DA  C5' C4'  SING N N 8", "DA  C5' H5'  SING N N 9",
                             "DA  C5' H5'' SING N N 10", "DA  C4' O4'  SING N N 11", "DA  C4' C3'  SING N N 12",
                             "DA  C4' H4'  SING N N 13", "DA  O4' C1'  SING N N 14", "DA  C3' O3'  SING N N 15",
                             "DA  C3' C2'  SING N N 16", "DA  C3' H3'  SING N N 17", "DA  O3' HO3' SING N N 18",
                             "DA  C2' C1'  SING N N 19", "DA  C2' H2'  SING N N 20", "DA  C2' H2'' SING N N 21",
                             "DA  C1' N9   SING N N 22", "DA  C1' H1'  SING N N 23", "DA  N9  C8   SING Y N 24",
                             "DA  N9  C4   SING Y N 25", "DA  C8  N7   DOUB Y N 26", "DA  C8  H8   SING N N 27",
                             "DA  N7  C5   SING Y N 28", "DA  C5  C6   SING Y N 29", "DA  C5  C4   DOUB Y N 30",
                             "DA  C6  N6   SING N N 31", "DA  C6  N1   DOUB Y N 32", "DA  N6  H61  SING N N 33",
                             "DA  N6  H62  SING N N 34", "DA  N1  C2   SING Y N 35", "DA  C2  N3   DOUB Y N 36",
                             "DA  C2  H2   SING N N 37", "DA  N3  C4   SING Y N 38"}),
          monomer_type::dna_linking));
      data_enum::add(std::make_shared<data>(
          "DC", 'C', //
          data::to_molecule({"DC  OP3 P    SING N N 1", "DC  OP3 HOP3 SING N N 2", "DC  P   OP1  DOUB N N 3",
                             "DC  P   OP2  SING N N 4", "DC  P   O5'  SING N N 5", "DC  OP2 HOP2 SING N N 6",
                             "DC  O5' C5'  SING N N 7", "DC  C5' C4'  SING N N 8", "DC  C5' H5'  SING N N 9",
                             "DC  C5' H5'' SING N N 10", "DC  C4' O4'  SING N N 11", "DC  C4' C3'  SING N N 12",
                             "DC  C4' H4'  SING N N 13", "DC  O4' C1'  SING N N 14", "DC  C3' O3'  SING N N 15",
                             "DC  C3' C2'  SING N N 16", "DC  C3' H3'  SING N N 17", "DC  O3' HO3' SING N N 18",
                             "DC  C2' C1'  SING N N 19", "DC  C2' H2'  SING N N 20", "DC  C2' H2'' SING N N 21",
                             "DC  C1' N1   SING N N 22", "DC  C1' H1'  SING N N 23", "DC  N1  C2   SING N N 24",
                             "DC  N1  C6   SING N N 25", "DC  C2  O2   DOUB N N 26", "DC  C2  N3   SING N N 27",
                             "DC  N3  C4   DOUB N N 28", "DC  C4  N4   SING N N 29", "DC  C4  C5   SING N N 30",
                             "DC  N4  H41  SING N N 31", "DC  N4  H42  SING N N 32", "DC  C5  C6   DOUB N N 33",
                             "DC  C5  H5   SING N N 34", "DC  C6  H6   SING N N 35"}),
          monomer_type::dna_linking));
      data_enum::add(std::make_shared<data>(
          "DG", 'G',
          data::to_molecule({"DG  OP3 P    SING N N 1", "DG  OP3 HOP3 SING N N 2", "DG  P   OP1  DOUB N N 3",
                             "DG  P   OP2  SING N N 4", "DG  P   O5'  SING N N 5", "DG  OP2 HOP2 SING N N 6",
                             "DG  O5' C5'  SING N N 7", "DG  C5' C4'  SING N N 8", "DG  C5' H5'  SING N N 9",
                             "DG  C5' H5'' SING N N 10", "DG  C4' O4'  SING N N 11", "DG  C4' C3'  SING N N 12",
                             "DG  C4' H4'  SING N N 13", "DG  O4' C1'  SING N N 14", "DG  C3' O3'  SING N N 15",
                             "DG  C3' C2'  SING N N 16", "DG  C3' H3'  SING N N 17", "DG  O3' HO3' SING N N 18",
                             "DG  C2' C1'  SING N N 19", "DG  C2' H2'  SING N N 20", "DG  C2' H2'' SING N N 21",
                             "DG  C1' N9   SING N N 22", "DG  C1' H1'  SING N N 23", "DG  N9  C8   SING Y N 24",
                             "DG  N9  C4   SING Y N 25", "DG  C8  N7   DOUB Y N 26", "DG  C8  H8   SING N N 27",
                             "DG  N7  C5   SING Y N 28", "DG  C5  C6   SING N N 29", "DG  C5  C4   DOUB Y N 30",
                             "DG  C6  O6   DOUB N N 31", "DG  C6  N1   SING N N 32", "DG  N1  C2   SING N N 33",
                             "DG  N1  H1   SING N N 34", "DG  C2  N2   SING N N 35", "DG  C2  N3   DOUB N N 36",
                             "DG  N2  H21  SING N N 37", "DG  N2  H22  SING N N 38", "DG  N3  C4   SING N N 39"}),
          monomer_type::dna_linking));
      data_enum::add(std::make_shared<data>(
          "DT", 'T', //
          data::to_molecule({"DT  OP3 P    SING N N 1", "DT  OP3 HOP3 SING N N 2", "DT  P   OP1  DOUB N N 3",
                             "DT  P   OP2  SING N N 4", "DT  P   O5'  SING N N 5", "DT  OP2 HOP2 SING N N 6",
                             "DT  O5' C5'  SING N N 7", "DT  C5' C4'  SING N N 8", "DT  C5' H5'  SING N N 9",
                             "DT  C5' H5'' SING N N 10", "DT  C4' O4'  SING N N 11", "DT  C4' C3'  SING N N 12",
                             "DT  C4' H4'  SING N N 13", "DT  O4' C1'  SING N N 14", "DT  C3' O3'  SING N N 15",
                             "DT  C3' C2'  SING N N 16", "DT  C3' H3'  SING N N 17", "DT  O3' HO3' SING N N 18",
                             "DT  C2' C1'  SING N N 19", "DT  C2' H2'  SING N N 20", "DT  C2' H2'' SING N N 21",
                             "DT  C1' N1   SING N N 22", "DT  C1' H1'  SING N N 23", "DT  N1  C2   SING N N 24",
                             "DT  N1  C6   SING N N 25", "DT  C2  O2   DOUB N N 26", "DT  C2  N3   SING N N 27",
                             "DT  N3  C4   SING N N 28", "DT  N3  H3   SING N N 29", "DT  C4  O4   DOUB N N 30",
                             "DT  C4  C5   SING N N 31", "DT  C5  C7   SING N N 32", "DT  C5  C6   DOUB N N 33",
                             "DT  C7  H71  SING N N 34", "DT  C7  H72  SING N N 35", "DT  C7  H73  SING N N 36",
                             "DT  C6  H6   SING N N 37"}),
          monomer_type::dna_linking));
      data_enum::add(std::make_shared<data>(
          "DI", 'I', //
          data::to_molecule({"DI  OP3 P    SING N N 1", "DI  OP3 HOP3 SING N N 2", "DI  P   OP1  DOUB N N 3",
                             "DI  P   OP2  SING N N 4", "DI  P   O5'  SING N N 5", "DI  OP2 HOP2 SING N N 6",
                             "DI  O5' C5'  SING N N 7", "DI  C5' C4'  SING N N 8", "DI  C5' H5'  SING N N 9",
                             "DI  C5' H5'' SING N N 10", "DI  C4' O4'  SING N N 11", "DI  C4' C3'  SING N N 12",
                             "DI  C4' H4'  SING N N 13", "DI  O4' C1'  SING N N 14", "DI  C3' O3'  SING N N 15",
                             "DI  C3' C2'  SING N N 16", "DI  C3' H3'  SING N N 17", "DI  O3' HO3' SING N N 18",
                             "DI  C2' C1'  SING N N 19", "DI  C2' H2'  SING N N 20", "DI  C2' H2'' SING N N 21",
                             "DI  C1' N9   SING N N 22", "DI  C1' H1'  SING N N 23", "DI  N9  C8   SING Y N 24",
                             "DI  N9  C4   SING Y N 25", "DI  C8  N7   DOUB Y N 26", "DI  C8  H8   SING N N 27",
                             "DI  N7  C5   SING Y N 28", "DI  C5  C6   SING N N 29", "DI  C5  C4   DOUB Y N 30",
                             "DI  C6  O6   DOUB N N 31", "DI  C6  N1   SING N N 32", "DI  N1  C2   SING N N 33",
                             "DI  N1  H1   SING N N 34", "DI  C2  N3   DOUB N N 35", "DI  C2  H2   SING N N 36",
                             "DI  N3  C4   SING N N 37"}),
          monomer_type::dna_linking));

      // add ribonucleotides
      data_enum::add(std::make_shared<data>(
          "A", 'A',
          data::to_molecule({"A   OP3 P    SING N N 1", "A   OP3 HOP3 SING N N 2", "A   P   OP1  DOUB N N 3",
                             "A   P   OP2  SING N N 4", "A   P   O5'  SING N N 5", "A   OP2 HOP2 SING N N 6",
                             "A   O5' C5'  SING N N 7", "A   C5' C4'  SING N N 8", "A   C5' H5'  SING N N 9",
                             "A   C5' H5'' SING N N 10", "A   C4' O4'  SING N N 11", "A   C4' C3'  SING N N 12",
                             "A   C4' H4'  SING N N 13", "A   O4' C1'  SING N N 14", "A   C3' O3'  SING N N 15",
                             "A   C3' C2'  SING N N 16", "A   C3' H3'  SING N N 17", "A   O3' HO3' SING N N 18",
                             "A   C2' O2'  SING N N 19", "A   C2' C1'  SING N N 20", "A   C2' H2'  SING N N 21",
                             "A   O2' HO2' SING N N 22", "A   C1' N9   SING N N 23", "A   C1' H1'  SING N N 24",
                             "A   N9  C8   SING Y N 25", "A   N9  C4   SING Y N 26", "A   C8  N7   DOUB Y N 27",
                             "A   C8  H8   SING N N 28", "A   N7  C5   SING Y N 29", "A   C5  C6   SING Y N 30",
                             "A   C5  C4   DOUB Y N 31", "A   C6  N6   SING N N 32", "A   C6  N1   DOUB Y N 33",
                             "A   N6  H61  SING N N 34", "A   N6  H62  SING N N 35", "A   N1  C2   SING Y N 36",
                             "A   C2  N3   DOUB Y N 37", "A   C2  H2   SING N N 38", "A   N3  C4   SING Y N 39"}),
          monomer_type::rna_linking));
      data_enum::add(std::make_shared<data>(
          "C", 'C',
          data::to_molecule({"C   OP3 P    SING N N 1", "C   OP3 HOP3 SING N N 2", "C   P   OP1  DOUB N N 3",
                             "C   P   OP2  SING N N 4", "C   P   O5'  SING N N 5", "C   OP2 HOP2 SING N N 6",
                             "C   O5' C5'  SING N N 7", "C   C5' C4'  SING N N 8", "C   C5' H5'  SING N N 9",
                             "C   C5' H5'' SING N N 10", "C   C4' O4'  SING N N 11", "C   C4' C3'  SING N N 12",
                             "C   C4' H4'  SING N N 13", "C   O4' C1'  SING N N 14", "C   C3' O3'  SING N N 15",
                             "C   C3' C2'  SING N N 16", "C   C3' H3'  SING N N 17", "C   O3' HO3' SING N N 18",
                             "C   C2' O2'  SING N N 19", "C   C2' C1'  SING N N 20", "C   C2' H2'  SING N N 21",
                             "C   O2' HO2' SING N N 22", "C   C1' N1   SING N N 23", "C   C1' H1'  SING N N 24",
                             "C   N1  C2   SING N N 25", "C   N1  C6   SING N N 26", "C   C2  O2   DOUB N N 27",
                             "C   C2  N3   SING N N 28", "C   N3  C4   DOUB N N 29", "C   C4  N4   SING N N 30",
                             "C   C4  C5   SING N N 31", "C   N4  H41  SING N N 32", "C   N4  H42  SING N N 33",
                             "C   C5  C6   DOUB N N 34", "C   C5  H5   SING N N 35", "C   C6  H6   SING N N 36"}),
          monomer_type::rna_linking));
      data_enum::add(std::make_shared<data>(
          "G", 'G',
          data::to_molecule({"G   OP3 P    SING N N 1", "G   OP3 HOP3 SING N N 2", "G   P   OP1  DOUB N N 3",
                             "G   P   OP2  SING N N 4", "G   P   O5'  SING N N 5", "G   OP2 HOP2 SING N N 6",
                             "G   O5' C5'  SING N N 7", "G   C5' C4'  SING N N 8", "G   C5' H5'  SING N N 9",
                             "G   C5' H5'' SING N N 10", "G   C4' O4'  SING N N 11", "G   C4' C3'  SING N N 12",
                             "G   C4' H4'  SING N N 13", "G   O4' C1'  SING N N 14", "G   C3' O3'  SING N N 15",
                             "G   C3' C2'  SING N N 16", "G   C3' H3'  SING N N 17", "G   O3' HO3' SING N N 18",
                             "G   C2' O2'  SING N N 19", "G   C2' C1'  SING N N 20", "G   C2' H2'  SING N N 21",
                             "G   O2' HO2' SING N N 22", "G   C1' N9   SING N N 23", "G   C1' H1'  SING N N 24",
                             "G   N9  C8   SING Y N 25", "G   N9  C4   SING Y N 26", "G   C8  N7   DOUB Y N 27",
                             "G   C8  H8   SING N N 28", "G   N7  C5   SING Y N 29", "G   C5  C6   SING N N 30",
                             "G   C5  C4   DOUB Y N 31", "G   C6  O6   DOUB N N 32", "G   C6  N1   SING N N 33",
                             "G   N1  C2   SING N N 34", "G   N1  H1   SING N N 35", "G   C2  N2   SING N N 36",
                             "G   C2  N3   DOUB N N 37", "G   N2  H21  SING N N 38", "G   N2  H22  SING N N 39",
                             "G   N3  C4   SING N N 40"}),
          monomer_type::rna_linking));
      data_enum::add(std::make_shared<data>(
          "U", 'U',
          data::to_molecule({"U   OP3 P    SING N N 1", "U   OP3 HOP3 SING N N 2", "U   P   OP1  DOUB N N 3",
                             "U   P   OP2  SING N N 4", "U   P   O5'  SING N N 5", "U   OP2 HOP2 SING N N 6",
                             "U   O5' C5'  SING N N 7", "U   C5' C4'  SING N N 8", "U   C5' H5'  SING N N 9",
                             "U   C5' H5'' SING N N 10", "U   C4' O4'  SING N N 11", "U   C4' C3'  SING N N 12",
                             "U   C4' H4'  SING N N 13", "U   O4' C1'  SING N N 14", "U   C3' O3'  SING N N 15",
                             "U   C3' C2'  SING N N 16", "U   C3' H3'  SING N N 17", "U   O3' HO3' SING N N 18",
                             "U   C2' O2'  SING N N 19", "U   C2' C1'  SING N N 20", "U   C2' H2'  SING N N 21",
                             "U   O2' HO2' SING N N 22", "U   C1' N1   SING N N 23", "U   C1' H1'  SING N N 24",
                             "U   N1  C2   SING N N 25", "U   N1  C6   SING N N 26", "U   C2  O2   DOUB N N 27",
                             "U   C2  N3   SING N N 28", "U   N3  C4   SING N N 29", "U   N3  H3   SING N N 30",
                             "U   C4  O4   DOUB N N 31", "U   C4  C5   SING N N 32", "U   C5  C6   DOUB N N 33",
                             "U   C5  H5   SING N N 34", "U   C6  H6   SING N N 35"}),
          monomer_type::rna_linking));
      data_enum::add(std::make_shared<data>(
          "I", 'I',
          data::to_molecule({"I   OP3 P    SING N N 1", "I   OP3 HOP3 SING N N 2", "I   P   OP1  DOUB N N 3",
                             "I   P   OP2  SING N N 4", "I   P   O5'  SING N N 5", "I   OP2 HOP2 SING N N 6",
                             "I   O5' C5'  SING N N 7", "I   C5' C4'  SING N N 8", "I   C5' H5'  SING N N 9",
                             "I   C5' H5'' SING N N 10", "I   C4' O4'  SING N N 11", "I   C4' C3'  SING N N 12",
                             "I   C4' H4'  SING N N 13", "I   O4' C1'  SING N N 14", "I   C3' O3'  SING N N 15",
                             "I   C3' C2'  SING N N 16", "I   C3' H3'  SING N N 17", "I   O3' HO3' SING N N 18",
                             "I   C2' O2'  SING N N 19", "I   C2' C1'  SING N N 20", "I   C2' H2'  SING N N 21",
                             "I   O2' HO2' SING N N 22", "I   C1' N9   SING N N 23", "I   C1' H1'  SING N N 24",
                             "I   N9  C8   SING Y N 25", "I   N9  C4   SING Y N 26", "I   C8  N7   DOUB Y N 27",
                             "I   C8  H8   SING N N 28", "I   N7  C5   SING Y N 29", "I   C5  C6   SING N N 30",
                             "I   C5  C4   DOUB Y N 31", "I   C6  O6   DOUB N N 32", "I   C6  N1   SING N N 33",
                             "I   N1  C2   SING N N 34", "I   N1  H1   SING N N 35", "I   C2  N3   DOUB N N 36",
                             "I   C2  H2   SING N N 37", "I   N3  C4   SING N N 38"}),
          monomer_type::rna_linking));

      return true;
    } // initialize()
    bool cc::_initialized(cc::initialize()); // initialize static variable
  } // namespace che
} // namespace biosim