#ifndef che_io_file_assembly_h
#define che_io_file_assembly_h

#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace io {
      // reads any file type that can be represented as an assembly, detecting the type based on the file extension
      class file_assembly {
      public:
        // default ctor
        file_assembly();

        // add reader function to read sse pool files
        void add_sse_pool_reader(assembly const &__reference);

        // reads an assembly from a given file
        assembly read(std::string const &__filename);

      private:
        using reader_function = std::function<assembly(std::string const &)>;
        std::map<std::string, reader_function> _readers; // map extension->reader_function
      }; // class file_assembly
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_assembly_h
