#ifndef che_io_file_sse_pool_h
#define che_io_file_sse_pool_h

#include "che/qs.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads sse pool and pdb files (only helix, sheet and seqres lines)
      class file_sse_pool {
      public:
        // creates a quarternary structure from a given file
        static qs read(std::string const &__filename);
      }; // class file_sse_pool
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_sse_pool_h
