#ifndef che_io_file_pdb_h
#define che_io_file_pdb_h

#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads pdb files
      class file_pdb {
      public:
        // read ps and ss from a given file
        static assembly read(std::string const &__filename);
      }; // class file_pdb
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_pdb_h
