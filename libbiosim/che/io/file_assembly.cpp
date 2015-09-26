#include "che/io/file_assembly.h"
#include "che/io/file_fasta.h"
#include "che/io/file_pdb.h"
#include "che/io/file_psipredv.h"
#include "che/io/file_sse_pool.h"
#include <boost/filesystem.hpp>
#include "tools/file.h"

namespace biosim {
  namespace che {
    namespace io {
      // default ctor
      file_assembly::file_assembly() : _readers() {
        _readers[".fasta"] = &file_fasta::read;
        _readers[".psipred_ss"] = &file_psipredv::read;
        _readers[".psipred_ss2"] = &file_psipredv::read;
        _readers[".rdbProf"] = &file_psipredv::read;
        _readers[".jufo9d_ss"] = &file_psipredv::read;
        _readers[".pdb"] = &file_pdb::read;
        // don't add sse_pool, sequence length is not defined (missing coil intervals), and overlap possible
      } // ctor

      // add reader function to read sse pool files
      void file_assembly::add_sse_pool_reader(assembly const &__reference) {
        _readers[".pool"] = [&](std::string const &__filename) -> assembly {
          DEBUG << "Extend length of ss sequence to "
                << __reference.get_first_structure("A").get_ss().get_sequence().size() << " to match reference";
          return file_sse_pool::read(__filename, __reference);
        }; // lambda
      } // add_sse_pool_reader()

      // reads an assembly from a given file
      assembly file_assembly::read(std::string const &__filename) {
        boost::filesystem::path filename(__filename); // get file extension
        std::map<std::string, reader_function>::const_iterator itr(_readers.find(filename.extension().string()));
        if(itr == _readers.end()) {
          throw tools::unknown_file_format(filename.filename().string());
        } // if

        DEBUG << "Found assembly file format reader for extension " << filename.extension();
        return itr->second.operator()(__filename); // return result based on function in map
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
