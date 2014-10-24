#ifndef che_io_file_qs_h
#define che_io_file_qs_h

#include "che/qs.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads any file type that can be represented as qs, detecting the type based on the file extension
      class file_qs {
      private:
        using reader_function = std::function<qs(std::string const &)>;
        std::map<std::string, reader_function> _readers; // map extension->reader_function

      public:
        // default ctor
        file_qs();

        // creates a qs from a given file
        qs read(std::string const &__filename);
      }; // class file_qs
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_qs_h
