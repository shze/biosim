#ifndef tools_file_h
#define tools_file_h

#include <string>
#include <list>

namespace biosim {
  namespace tools {
    // file io specific functions
    class file {
    public:
      // read full content of file of given filename, return a list of strings for each line
      static std::list<std::string> read_to_string_list(std::string const &__filename);
    }; // class file
  } // namespace tools
} // namespace biosim

#endif // tools_file_h
