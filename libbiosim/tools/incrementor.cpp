#include "tools/incrementor.h"
#include <stdexcept>
#include "tools/log.h"

namespace biosim {
  namespace tools {
    // ctor taking alphabet; default ctor
    incrementor::incrementor(std::string __alphabet)
        : _alphabet(__alphabet), _first_char_is_neutral(__alphabet[0] == '0') {}
    // incrementing function
    std::string incrementor::next(std::string __s, size_t increment) const {
      int pos = __s.length() - 1; // start at last char; use int because we're decreasing it, and it can go < 0

      while(increment) {
        size_t current_letter_pos = _alphabet.find(__s[pos]);
        if(current_letter_pos == std::string::npos) {
          throw std::range_error("string contains letter outside of alphabet range");
        } else if(current_letter_pos + increment >= _alphabet.size()) { // if we have carryover
          __s[pos] = _alphabet[(current_letter_pos + increment) % _alphabet.size()];
          increment = (current_letter_pos + increment) / _alphabet.size();
        } else { // no carryover
          __s[pos] = _alphabet[current_letter_pos + increment];
          increment = 0;
        }
        pos--;
        if(increment && pos < 0) { // if we run out of string, add at the beginning
          __s.insert(__s.begin(), _first_char_is_neutral ? _alphabet[increment] : _alphabet[increment - 1]);
          increment = 0; // technically we also need to do pos++, but this is always the last step, so skip this
        } // if
      } // while

      return __s;
    } // next()

    // (static) get uppercase alphabet
    std::string const &incrementor::get_uppercase_alphabet() {
      static std::string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
      return letters;
    } // get_uppercase_alphabet()
  } // namespace tools
} // namespace biosim
