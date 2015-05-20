#ifndef tools_incrementor_h
#define tools_incrementor_h

#include <stdexcept>
#include <vector>
#include <algorithm>

namespace biosim {
  namespace tools {
    // increments the elements of a given container based on the alphabet(s) given to ctor
    template <class T>
    class incrementor {
    public:
      // ctor taking alphabets, at least one is required; default ctor
      incrementor(std::vector<std::vector<typename T::value_type>> const &__alphabets = {get_uppercase_alphabet()})
          : _alphabets(__alphabets), _fixed_length(__alphabets.size() > 1) {
        if(__alphabets.empty()) {
          throw std::invalid_argument("at least one alphabet is needed");
        } // if
      } // ctor
      // incrementing function
      T next(T __container, size_t increment = 1) const {
        int pos(__container.size() - 1); // start at last element; use int because we're decreasing, and it can go < 0
        if(_fixed_length && _alphabets.size() != __container.size()) {
          throw std::out_of_range("input length differs from alphabet count for fixed length input");
        } // if

        while(increment) {
          std::vector<typename T::value_type> alphabet(_alphabets.size() == 1 ? _alphabets[0] : _alphabets[pos]);
          bool first_char_is_neutral(alphabet[0] == '0' || alphabet[0] == 0);

          size_t current_letter_pos(std::find(alphabet.begin(), alphabet.end(), __container[pos]) - alphabet.begin());
          if(current_letter_pos == alphabet.size()) {
            throw std::out_of_range("input contains letter outside of alphabet range");
          } else if(current_letter_pos + increment >= alphabet.size()) { // if we have carryover
            __container[pos] = alphabet[(current_letter_pos + increment) % alphabet.size()];
            increment = (current_letter_pos + increment) / alphabet.size();
          } else { // no carryover
            __container[pos] = alphabet[current_letter_pos + increment];
            increment = 0;
          } // else
          pos--;
          if(increment && pos < 0) { // if we run out of elements, add at the beginning if allowed
            if(_fixed_length) {
              throw std::overflow_error("adding increment causes overflow for fixed length input");
            } // if
            __container.insert(__container.begin(),
                               first_char_is_neutral ? alphabet[increment] : alphabet[increment - 1]);
            increment = 0; // technically we also need to do pos++, but this is always the last step, so skip this
          } // if
        } // while

        return __container;
      } // next()

      // get uppercase alphabet
      static std::vector<char> get_uppercase_alphabet() {
        static std::string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        return std::vector<char>(letters.begin(), letters.end());
      } // get_uppercase_alphabet()
    private:
      std::vector<std::vector<typename T::value_type>> _alphabets; // the ordered alphabet for incrementing
      bool _fixed_length; // if the increment function does not extend the length of a given argument
    }; // class incrementor
  } // namespace tools
} // namespace biosim

#endif // tools_incrementor_h
