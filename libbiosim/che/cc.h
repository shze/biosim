#ifndef che_cc_h
#define che_cc_h

#include <map>
#include <list>
#include <memory>
#include "tools/enumerate.h"

namespace biosim {
  namespace che {
    // stores a chemical component based on the CCD, see: http://www.wwpdb.org/ccd.html;
    // focuses on amino acids, which it defaults to, since the CCD does not contain all one letter codes for DNA and
    // RNA, see http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml, only the following symbols are in the CCD:
    // DNA: DA/A, DC/C, DG/G, DT/T, DI/I, DN/N; RNA: A/A,  C/C,  G/G,  U/U,  I/I,  N/N (id/id_char for each symbol).
    // design:
    // oo design containing all information about a cc, not just a char; cc uses a pimpl design, struct data is the
    // enumerated implementation; all cc ctors either find the correct data or throw.
    class cc {
    public:
      // specificity defined as enum for easier comparison; only one specificity allowed, no combinations.
      enum specificity_type { specificity_defined, specificity_unknown, specificity_gap };
      // type of monomer, defined as chemical component type in the CCD; used to differentiate symbols of different
      // subgroups with overlapping identifier chars i.e. one letter codes.
      enum monomer_type { non_polymer, l_peptide_linking, dna_linking, rna_linking };
      using weight_map = std::map<std::string, double>; // simplify naming

      // construct cc with id; id is assumed unique; throws if not found
      cc(std::string __id);
      // get identifier
      std::string get_identifier() const;
      // get identifier char, i.e. single letter code
      char get_identifier_char() const;
      // get specificity
      specificity_type get_specificity() const;
      // return if an unknown (unspecified) symbol; only for unknown, not for gap.
      bool is_unknown() const;
      // return if gap
      bool is_gap() const;
      // get monomer type
      monomer_type get_monomer_type() const;
      // get profile weights
      weight_map get_weights() const;

      // get list of all identifiers; static public interface
      static std::list<std::string> get_id_list();

    private:
      // contains the cc data; type that is enumerated
      struct data {
        // ctor for primary unique-single-letter AA
        data(std::string __id, char __id_char, specificity_type __specificity = specificity_defined,
             monomer_type __monomer_type = l_peptide_linking);
        // ctor for profile derived-single-letter AA
        data(std::string __id, std::string __base_id);
        // ctor for profile unique-single-letter AA; ctor does not check if base ids have same specificity and
        // monomer values; it assumes they are consistent and uses the values of the first base id.
        data(std::string __id, char __id_char, std::string __base_id1, std::string __base_id2);

        std::string _id; // identifier (three letter code)
        char _id_char; // identifier char (one letter code)
        specificity_type _specificity; // specificity
        monomer_type _monomer_type; // type of monomer (chemical component type)
        weight_map _weights; // weights of this data object based on symbols in the enumerated instance set
      }; // struct data

      // define identifier function for use with enumerate; friend and not part of the class
      friend std::string identifier(std::shared_ptr<data const> const &__data_ptr) { return __data_ptr->_id; }

      using data_enum = tools::enumerate<std::shared_ptr<data>>; // simplify naming
      // find symbol with id; id is assumed unique; throws if no object was found
      static std::shared_ptr<data> find(std::string __id);
      static bool _initialized; // static variable to initialize data_enum with a minimal set of data objects
      // constructs data objects of the enumerated instance set, i.e. all single letter AA, see:
      // http://www.ebi.ac.uk/pdbe-srv/pdbechem/
      // http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bia01as00/full
      // http://www.chem.qmul.ac.uk/iupac/AminoAcid/
      // http://www.rcsb.org/pdb/101/static101.do?p=education_discussion/Looking-at-Structures/sequence.html
      // http://proline.bic.nus.edu.sg/~asif/AminoAcids.ppt
      static bool initialize();

      std::shared_ptr<data const> _impl; // ptr to implementation; shared_ptr is acceptable, b/c it points to const data
    }; // class cc
  } // namespace che
} // namespace biosim

#endif // che_cc_h
