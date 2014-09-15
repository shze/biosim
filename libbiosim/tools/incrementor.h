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
    }; // class incrementor
  } // namespace tools
} // namespace biosim

#endif // tools_incrementor_h
