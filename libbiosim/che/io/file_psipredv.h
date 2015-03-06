#ifndef che_io_file_psipredv_h
#define che_io_file_psipredv_h

#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads psipred vformat files, see: http://webdocs.cs.ualberta.ca/~lingroup/Software/PSAtip/format.php#ss
      class file_psipredv {
      public:
        // reads ps and ss from a given file
        static assembly read(std::string const &__filename);
      }; // class file_psipredv
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_psipredv_h
