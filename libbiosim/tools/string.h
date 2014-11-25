#ifndef tools_string_h
#define tools_string_h

#include <string>
#include <functional>
#include <set>

namespace biosim {
  namespace tools {
    // function type taking a char and returning if the char is of a certain class
    using char_function = std::function<int(int)>;
    // classify function that always classifies chars as not in this class, i.e. returns always false
    int nothing(int ch);

    // returns the unique char set without the ignored chars from a given string
    std::set<char> get_unique_char_set_ignore(std::string const &__s, char_function const &__ignore);
    // returns the unique char set from a given string
    std::set<char> get_unique_char_set(std::string const &__s);
    // returns the unique char string without the ignored chars from a given string
    std::string get_unique_char_string_ignore(std::string const &__s, char_function const &__ignore);
    // returns the unique char string from a given string
    std::string get_unique_char_string(std::string const &__s);
  } // namespace tools
} // namespace biosim

#endif // tools_string_h
