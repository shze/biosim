#ifndef che_io_file_pssm_h
#define che_io_file_pssm_h

#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads blast pssm files
      class file_pssm {
      public:
        // read ps from profile data of a given file
        static assembly read(std::string const &__filename);

      private:
        // intermediate data structure to store information about an amino acid position
        struct pssm_amino_acid_data {
          // ctor taking aa id, aa type, information, and pseudocount weight, to simplify setting all members
          pssm_amino_acid_data(size_t __aa_id, char __aa_type, double __information, double __pseudocount_weight)
              : aa_id(__aa_id),
                aa_type(__aa_type),
                score(),
                percent(),
                information(__information),
                pseudocount_weight(__pseudocount_weight) {}

          size_t aa_id; // amino acid id
          char aa_type; // amino acid type (single letter code)
          std::vector<int> score; // score values
          std::vector<size_t> percent; // profile values
          double information; // information content of this position
          double pseudocount_weight; // relative weight of gapless real matches to pseudocounts
        }; // struct pssm_amino_acid_data
        // intermediate data structure to store all pssm information
        struct pssm_data {
          std::vector<pssm_amino_acid_data> aa_data; // amino acid data
        }; // struct pssm_data

        // reads data from a pssm ascii file into an intermediate data structure
        static pssm_data read_pssm_ascii(std::string const &__filename);
      }; // class file_pssm
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_pssm_h
