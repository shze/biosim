#include "tools/file.h"
#include <fstream>
#include <boost/algorithm/string.hpp>

namespace biosim {
  namespace tools {
    // (static) read full content of file of given filename, return a list of strings for each line
    std::list<std::string> file::read_to_string_list(std::string const &__filename) {
      std::ifstream filestream(__filename.c_str()); // create stream, automatically closed when leaving scope/function
      // being unable to open a file could simply result in an empty string list, but it's an important difference for
      // the user if the file could not be opened or if it was empty; distinguish by throwing an exception for the first
      if(!filestream.is_open()) {
        throw std::ios_base::failure("ios_base failure: could not open file " + __filename);
      }

      // read the file line by line, collect content
      std::string buffer;
      std::list<std::string> lines;
      while(filestream.good()) {
        if(std::getline(filestream, buffer).good()) {
          boost::trim_right(buffer);
          lines.push_back(buffer);
        }
      }

      return lines;
    } // read_to_string_list()
  } // namespace tools
} // namespace biosim
