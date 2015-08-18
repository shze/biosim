#ifndef che_cc_h
#define che_cc_h

#include <map>
#include <list>
#include <memory>
#include "tools/enumerate.h"
#include "che/molecule.h"

namespace biosim {
  namespace che {
    // exception specific for cc
    struct cc_data_not_found : std::runtime_error {
      explicit cc_data_not_found(std::string const &__s) : std::runtime_error("cc data not found: " + __s) {}
    }; // struct cc_data_not_found

    // stores a chemical component based on the CCD, see: http://www.wwpdb.org/ccd.html;
    // focuses on amino acids, which it defaults to, since the CCD does not contain all one letter codes for DNA and
    // RNA, see http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml, only the following compounds are in the CCD:
    // DNA: DA/A, DC/C, DG/G, DT/T, DI/I, DN/N; RNA: A/A,  C/C,  G/G,  U/U,  I/I,  N/N (id/id_char for each compound).
    // design:
    // oo design containing all information about a cc, not just a char; cc uses a pimpl design, struct data is the
    // enumerated implementation; all cc ctors either find the correct data or throw.
    class cc {
    public:
      // specificity defined as enum for easier comparison; only one specificity allowed, no combinations.
      // primary: all unique-single-letter AA (20 'natural' and 3 ambiguous profile AA);
      // secondary: all profile-but-not-primary and all derived AA
      // (profile AA: AA with _weights.size > 1; derived AA: _id != _weights.first._id)
      enum class specificity_type { primary, secondary, unknown, gap };
      // type of monomer, defined as chemical component type in the CCD; used to differentiate compounds of different
      // subgroups with overlapping identifier chars i.e. one letter codes.
      enum class monomer_type { non_polymer, l_peptide_linking, dna_linking, rna_linking };

      using weight_map = std::map<std::string, double>; // simplify naming

      // ctor from id; id is assumed unique; throws if not found
      explicit cc(std::string const &__id);
      // ctor from id_char and monomer_type, assumed unique; use const ref for args b/c they're no sink args
      explicit cc(char const &__id_char, monomer_type const &__monomer = monomer_type::l_peptide_linking);
      // ctor from weight set
      explicit cc(weight_map __weights);

      // get identifier
      std::string get_identifier() const;
      // get identifier char, i.e. single letter code
      char get_identifier_char() const;
      // get specificity
      specificity_type get_specificity() const;
      // return if an unknown (unspecified) compound; only for unknown, not for gap.
      bool is_unknown() const;
      // return if gap
      bool is_gap() const;
      // get monomer type
      monomer_type get_monomer_type() const;
      // get profile weights
      weight_map get_weights() const;
      // get molecule
      molecule const &get_molecule() const;
      // get determined atoms
      std::set<atom> const &get_determined_atoms() const;

      // get list of all identifiers; static public interface
      static std::list<std::string> get_id_list();
      // create a string of identifier chars; static public interface
      static std::string get_identifier_char_string();

    private:
      // contains the cc data; type that is enumerated
      struct data {
        // ctor for primary unique-single-letter AA
        data(std::string __id, char __id_char, molecule __molecule,
             monomer_type __monomer = monomer_type::l_peptide_linking);
        // ctor for profile unique-single-letter AA; ctor does not check if base ids have same specificity and
        // monomer values; it assumes they are consistent and uses the values of the first base id.
        data(std::string __id, char __id_char, std::string __base_id1, std::string __base_id2, molecule __molecule);
        // ctor for derived AA
        data(std::string __id, std::string __base_id, molecule __molecule);
        // ctor for profile, unknown, and gap AA
        data(std::string __id, char __id_char, cc::specificity_type __specificity, cc::monomer_type __monomer,
             weight_map __weights, molecule __molecule = molecule());

        std::string _id; // identifier (three letter code)
        char _id_char; // identifier char (one letter code)
        specificity_type _specificity; // specificity
        monomer_type _monomer; // type of monomer (chemical component type)
        weight_map _weights; // weights of this data object based on compounds in the enumerated instance set
        molecule _molecule; // the complete atom and bond structure that this cc can have

        // convert an atom-bond definition into a molecule; the atom-bond definition must be in CCD bond format:
        // comp_id atom_id_1 atom_id_2 value_order pdbx_aromatic_flag pdbx_stereo_config pdbx_ordinal
        static molecule to_molecule(std::list<std::string> __atom_bond_definition);
        // normalize the given weight map (so that the sum of all weights is 1)
        static weight_map normalize_weights(weight_map __weights);
      }; // struct data

      // define identifier function for use with enumerate; friend and not part of the class
      friend std::string const &identifier(std::shared_ptr<data const> const &__data_ptr) { return __data_ptr->_id; }

      using data_enum = tools::enumerate<std::shared_ptr<data const>>; // simplify naming
      // find compound with id; id is assumed unique
      static std::shared_ptr<data const> find(std::string const &__id);
      // find compound with id_char and monomer_type; combination is assumed unique
      static std::shared_ptr<data const> find(char const &__id_char, monomer_type const &__monomer);
      // constructs data objects of the enumerated instance set, i.e. all single letter AA, see:
      // http://www.ebi.ac.uk/pdbe-srv/pdbechem/
      // http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bia01as00/full
      // http://www.chem.qmul.ac.uk/iupac/AminoAcid/
      // http://www.rcsb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/sequence.html
      // http://proline.bic.nus.edu.sg/~asif/AminoAcids.ppt
      static bool initialize();
      static bool _initialized; // static variable to initialize data_enum with a minimal set of data objects

      std::shared_ptr<data const> _impl; // ptr to implementation; shared_ptr is acceptable, b/c it points to const data
      std::set<atom> _atoms; // determined atoms in this cc; subset of the atoms in _impl->_molecule
    }; // class cc
  } // namespace che
} // namespace biosim

#endif // che_cc_h
