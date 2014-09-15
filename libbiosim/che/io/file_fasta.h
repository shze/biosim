#ifndef che_io_file_fasta_h
#define che_io_file_fasta_h

#include "che/complex.h"
#include <set>

namespace biosim {
  namespace che {
    namespace io {
      // reads fasta files, including multifasta; see: http://en.wikipedia.org/wiki/FASTA_format;
      // http://www.ebi.ac.uk/help/formats.html#fasta; http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
      // format definition for multifasta (regular grammar, can be parsed by regular expression):
      // multifasta=(;sequence_id\n)+|(>sequence_id\n)(sequence\n)+(>sequence_id\n(sequence\n)+)*
      // sequence_id=.*
      // sequence=[aa_chars]{1..80}|[na_chars]{1..80}
      // aa_chars=abcdefghiklmnprstuvwxyzABCDEFGHIKLMNPRSTUVWXYZ*-
      // na_chars=abcdghkmnrstuvwyABCDGHKMNRSTUVWY-
      class file_fasta {
      public:
        // reads sequences from a given file
        static complex read(std::string const &__filename);

      private:
        // return the set of unique chars in the given string without whitespace
        static std::set<char> get_unique_chars_ignore_whitespace(std::string const &__s);
        // convert a string of id_chars to a ps
        static ps convert_to_ps(std::string __s);
      }; // class file_fasta
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_fasta_h
