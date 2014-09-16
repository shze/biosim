#ifndef tools_incrementor_h
#define tools_incrementor_h

#include <string>

namespace biosim {
  namespace tools {
    // increments a given string based on the ordered alphabet given to ctor
    class incrementor {
    public:
      // ctor taking alphabet; default ctor
      incrementor(std::string __alphabet = get_uppercase_alphabet());
      // incrementing function
      std::string next(std::string __s, size_t increment = 1) const;

      // get uppercase alphabet
      static std::string const &get_uppercase_alphabet();

    private:
      std::string _alphabet; // the ordered alphabet for incrementing
      // neutral element is not used for carryover, e.g. alphabet = 01: next(0) = 10, not 00;
      // non-neutral elements are used e.g. alphabet = AB: next(B) = AA, not BA
      bool _first_char_is_neutral;
    }; // class incrementor
  } // namespace tools
} // namespace biosim

#endif // tools_incrementor_h
