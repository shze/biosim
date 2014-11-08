#include "che/io/file_qs.h"
#include "che/io/file_psipredv.h"
#include "che/io/file_sse_pool.h"
#include <boost/filesystem.hpp>
#include "tools/file.h"

namespace biosim {
  namespace che {
    namespace io {
      // default ctor
      file_qs::file_qs() : _readers() {
        _readers[".psipred_ss"] = &file_psipredv::read;
        _readers[".psipred_ss2"] = &file_psipredv::read;
        _readers[".rdbProf"] = &file_psipredv::read;
        _readers[".jufo9d_ss"] = &file_psipredv::read;
        // add pdb, sequence length defined, no overlapping sse; cast b/c compiler does not know which read() to use
        _readers[".pdb"] = static_cast<qs (*)(std::string const &)>(&file_sse_pool::read);
        // don't add sse_pool, sequence length is not defined (missing coil intervals), and overlap possible
      } // ctor

      // add reader function to read sse pool files
      void file_qs::add_sse_pool_reader(qs const &__reference) {
        _readers[".pool"] = [&](std::string const &__filename) -> qs {
          DEBUG << "Extend length of ss sequence to " << __reference.get_ss("A").get_sequence().size()
                << " to match reference";
          return file_sse_pool::read(__filename, __reference);
        }; // lambda
      } // add_sse_pool_reader()

      // creates a qs from a given file
      qs file_qs::read(std::string const &__filename) {
        boost::filesystem::path filename(__filename); // get file extension
        std::map<std::string, reader_function>::const_iterator itr(_readers.find(filename.extension().string()));
        if(itr == _readers.end()) {
          throw tools::unknown_file_format(filename.filename().string());
        } // if

        DEBUG << "Found qs file format reader for extension " << filename.extension();
        return itr->second.operator()(__filename); // return result based on function in map
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
