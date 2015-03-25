#ifndef che_io_file_sse_pool_h
#define che_io_file_sse_pool_h

#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads sse pool files (only helix and sheet lines)
      class file_sse_pool {
      public:
        // read ps and ss from a given file
        static assembly read(std::string const &__filename);
        // read ps and ss from a given file and extends the ps to the given reference length
        static assembly read(std::string const &__filename, assembly const &__reference);
      }; // class file_sse_pool
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_sse_pool_h
