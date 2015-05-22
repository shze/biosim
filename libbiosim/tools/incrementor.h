#ifndef tools_incrementor_h
#define tools_incrementor_h

#include <stdexcept>
#include <vector>
#include <algorithm>
#include "tools/log.h"

namespace biosim {
  namespace tools {
    // increments the elements of a given container based on the alphabet(s) given to ctor
    template <class T>
    class incrementor {
    public:
      // ctor taking alphabets, at least one is required; default ctor
      explicit incrementor(std::vector<std::vector<typename T::value_type>> const &__alphabets = {
                               get_uppercase_alphabet()})
          : _alphabets(__alphabets), _fixed_length(__alphabets.size() > 1) {
        if(__alphabets.empty()) {
          throw std::invalid_argument("at least one alphabet is needed");
        } // if
      } // ctor
      // incrementing function
      T next(T __container, size_t increment = 1) const {
        if(_fixed_length && _alphabets.size() != __container.size()) { // check container size if fixed length container
          throw std::out_of_range("input length differs from alphabet count for fixed length input");
        } // if
        for(size_t pos(0); pos < __container.size(); ++pos) { // check all positions have valid letters
          std::vector<typename T::value_type> alphabet(_fixed_length ? _alphabets[pos] : _alphabets[0]);
          if(std::find(alphabet.begin(), alphabet.end(), __container[pos]) == alphabet.end()) {
            throw std::out_of_range("input contains letter outside of alphabet range");
          } // if
        } // for
        int pos(__container.size() - 1); // start at last element; use int because we're decreasing, and it can go < 0

        while(increment) {
          std::vector<typename T::value_type> alphabet(_fixed_length ? _alphabets[pos] : _alphabets[0]);
          size_t current_letter_pos(std::find(alphabet.begin(), alphabet.end(), __container[pos]) - alphabet.begin());
          if(current_letter_pos + increment >= alphabet.size()) { // if we have carryover
            DEBUG << "Increasing element at position " << pos << ", keeping carryover for next element";
            __container[pos] = alphabet[(current_letter_pos + increment) % alphabet.size()];
            increment = (current_letter_pos + increment) / alphabet.size();
          } else { // no carryover
            DEBUG << "Increasing element at position " << pos << ", done incrementing";
            __container[pos] = alphabet[current_letter_pos + increment];
            increment = 0;
          } // else
          pos--;
          if(increment && pos < 0) {
            if(_fixed_length) {
              throw std::overflow_error("adding increment causes overflow for fixed length input");
            } // if
            DEBUG << "Inserting new element at the front of the container";
            bool first_char_is_neutral(alphabet.size() > 0 && (alphabet[0] == '0' || alphabet[0] == 0));
            __container.insert(__container.begin(), first_char_is_neutral ? alphabet[1] : alphabet[0]);
            pos++;
            increment--; // decrease by 1, b/c we inserted the first element (b/c sometimes we have no neutral element)
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
