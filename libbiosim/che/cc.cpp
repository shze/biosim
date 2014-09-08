#include "che/cc.h"

namespace biosim {
  namespace che {
    // construct cc with id; id is assumed unique; throws if not found
    cc::cc(std::string __id) : _impl(find(__id)) {}
    // get identifier
    std::string cc::get_identifier() const { return _impl->_id; }
    // get identifier char, i.e. single letter code
    char cc::get_identifier_char() const { return _impl->_id_char; }
    // get specificity
    cc::specificity_type cc::get_specificity() const { return _impl->_specificity; }
    // return if an unknown (unspecified) symbol; only for unknown, not for gap.
    bool cc::is_unknown() const { return _impl->_specificity == specificity_unknown; }
    // return if gap
    bool cc::is_gap() const { return _impl->_specificity == specificity_gap; }
    // get monomer type
    cc::monomer_type cc::get_monomer_type() const { return _impl->_monomer_type; }
    // get profile weights
    cc::weight_map cc::get_weights() const { return _impl->_weights; }

    // (static) get list of all identifiers; static public interface
    std::list<std::string> cc::get_id_list() {
      std::list<std::string> ids;
      for(auto p : data_enum::get_instances()) {
        ids.push_back(p.first);
      }
      return ids;
    } // get_id_list()

    // ctor for primary unique-single-letter AA; __specificity and __monomer_type have defaults
    cc::data::data(std::string __id, char __id_char, cc::specificity_type __specificity,
                   cc::monomer_type __monomer_type)
        : _id(__id), _id_char(__id_char), _specificity(__specificity), _monomer_type(__monomer_type), _weights() {
      _weights = {{_id, 1.0}};
    } // data ctor
    // ctor for profile derived-single-letter AA
    cc::data::data(std::string __id, std::string __base_id)
        : _id(__id), _id_char(), _specificity(), _monomer_type(), _weights() {
      std::shared_ptr<data> base_ptr(find(__base_id)); // find base id, use it to set the following values
      _id_char = base_ptr->_id_char;
      _specificity = base_ptr->_specificity;
      _monomer_type = base_ptr->_monomer_type;
      _weights = {{__base_id, 1.0}};
    } // data ctor
    // ctor for profile unique-single-letter AA; ctor does not check if base ids have same specificity and
    // monomer values; it assumes they are consistent and uses the values of the first base id.
    cc::data::data(std::string __id, char __id_char, std::string __base_id1, std::string __base_id2)
        : _id(__id), _id_char(__id_char), _specificity(), _monomer_type(), _weights() {
      std::shared_ptr<data> base_ptr1(find(__base_id1));
      std::shared_ptr<data> base_ptr2(find(__base_id2));
      _specificity = base_ptr1->_specificity;
      _monomer_type = base_ptr1->_monomer_type;
      _weights = {{__base_id1, 0.5}, {__base_id2, 0.5}};
    } // data ctor

    // (static) find symbol with id; id is assumed unique; throws if no object was found
    std::shared_ptr<cc::data> cc::find(std::string __id) { return data_enum::get_instances().at(__id).get_object(); }
    // (static) constructs data objects of the enumerated instance set, i.e. all single letter AA, see:
    // http://www.ebi.ac.uk/pdbe-srv/pdbechem/
    // http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bia01as00/full
    // http://www.chem.qmul.ac.uk/iupac/AminoAcid/
    // http://www.rcsb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/sequence.html
    // http://proline.bic.nus.edu.sg/~asif/AminoAcids.ppt
    bool cc::initialize() {
      // 20 primary single-letter AA
      data_enum::add(std::make_shared<data>("ALA", 'A'));
      data_enum::add(std::make_shared<data>("CYS", 'C'));
      data_enum::add(std::make_shared<data>("ASP", 'D'));
      data_enum::add(std::make_shared<data>("GLU", 'E'));
      data_enum::add(std::make_shared<data>("PHE", 'F'));
      data_enum::add(std::make_shared<data>("GLY", 'G'));
      data_enum::add(std::make_shared<data>("HIS", 'H'));
      data_enum::add(std::make_shared<data>("ILE", 'I'));
      data_enum::add(std::make_shared<data>("LYS", 'K'));
      data_enum::add(std::make_shared<data>("LEU", 'L'));
      data_enum::add(std::make_shared<data>("MET", 'M'));
      data_enum::add(std::make_shared<data>("ASN", 'N'));
      data_enum::add(std::make_shared<data>("PRO", 'P'));
      data_enum::add(std::make_shared<data>("GLN", 'Q'));
      data_enum::add(std::make_shared<data>("ARG", 'R'));
      data_enum::add(std::make_shared<data>("SER", 'S'));
      data_enum::add(std::make_shared<data>("THR", 'T'));
      data_enum::add(std::make_shared<data>("VAL", 'V'));
      data_enum::add(std::make_shared<data>("TRP", 'W'));
      data_enum::add(std::make_shared<data>("TYR", 'Y'));
      // 5 profile single-letter AA, well only 3
      data_enum::add(std::make_shared<data>("ASX", 'B', "ASN", "ASP"));
      data_enum::add(std::make_shared<data>("GLX", 'Z', "GLN", "GLU"));
      data_enum::add(std::make_shared<data>("XLE", 'J', "LEU", "ILE"));
      // not adding the following with their single letters OU b/c these are used by other compounds in the CCD
      // PYRROLYSINE single-letter: K or O; three-letter: PYH or PYL;
      // SELENOCYSTEINE single-letter: C or U; three-letter: CSE or SEC;
      // unknown AA (special primary single-letter AA)
      data_enum::add(std::make_shared<data>("UNK", 'X', specificity_unknown)); // l_peptide_linking by default
      // gap (special primary single-letter AA); there seem to be no three letter gap id; used one letter codes: -~.
      data_enum::add(std::make_shared<data>("---", '-', specificity_gap));

      // add some common non-primary AAs
      data_enum::add(std::make_shared<data>("PYH", "LYS")); // PYRROLYSINE
      data_enum::add(std::make_shared<data>("CSE", "CYS")); // SELENOCYSTEINE
      data_enum::add(std::make_shared<data>("MSE", "MET")); // SELENOMETHIONINE

      // add deoxyribonucleotides
      data_enum::add(std::make_shared<data>("DA", 'A', specificity_defined, dna_linking));
      data_enum::add(std::make_shared<data>("DC", 'C', specificity_defined, dna_linking));
      data_enum::add(std::make_shared<data>("DG", 'G', specificity_defined, dna_linking));
      data_enum::add(std::make_shared<data>("DT", 'T', specificity_defined, dna_linking));
      data_enum::add(std::make_shared<data>("DI", 'I', specificity_defined, dna_linking));

      // add ribonucleotides
      data_enum::add(std::make_shared<data>("A", 'A', specificity_defined, rna_linking));
      data_enum::add(std::make_shared<data>("C", 'C', specificity_defined, rna_linking));
      data_enum::add(std::make_shared<data>("G", 'G', specificity_defined, rna_linking));
      data_enum::add(std::make_shared<data>("U", 'U', specificity_defined, rna_linking));
      data_enum::add(std::make_shared<data>("I", 'I', specificity_defined, rna_linking));

      return true;
    } // initialize()

    bool cc::_initialized(cc::initialize()); // initialize static variable
  }
}