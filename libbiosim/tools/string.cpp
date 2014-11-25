#include "tools/string.h"

namespace biosim {
  namespace tools {
    // classify function that always classifies chars as not in this class, i.e. returns always false
    int nothing(int ch) { return 0; }

    // returns the unique char set without the ignored chars from a given string
    std::set<char> get_unique_char_set_ignore(std::string const &__s, char_function const &__ignore) {
      std::set<char> unique;
      for(auto const &c : __s) {
        if(!__ignore(c)) {
          unique.insert(c);
        } // if
      } // for

      return unique;
    } // get_unique_char_set_ignore()
    // returns the unique char set from a given string
    std::set<char> get_unique_char_set(std::string const &__s) { return get_unique_char_set_ignore(__s, &nothing); }
    // returns the unique char string without the ignored chars from a given string
    std::string get_unique_char_string_ignore(std::string const &__s, char_function const &__ignore) {
      std::set<char> char_set(get_unique_char_set_ignore(__s, __ignore));
      return std::string(char_set.begin(), char_set.end());
    } // get_unique_char_string_ignore()
    // returns the unique char string from a given string
    std::string get_unique_char_string(std::string const &__s) { return get_unique_char_string_ignore(__s, &nothing); }
  } // namespace tools
} // namespace biosim
