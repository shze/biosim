#ifndef che_io_file_fasta_h
#define che_io_file_fasta_h

#include "che/assembly.h"
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
        static assembly read(std::string const &__filename);

      private:
        // convert a string of id_chars to a ps
        static ps convert_to_ps(std::string __s);
        // convert a string of id_chars to an ss
        static ss convert_to_ss(std::string __s);
        // shorten the file storage location to the last two parts
        static std::string shorten_file_storage(std::string __s);
      }; // class file_fasta
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_fasta_h
