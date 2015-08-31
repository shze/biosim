#ifndef tools_incrementor_h
#define tools_incrementor_h

#include <stdexcept>
#include <vector>
#include <algorithm>

namespace biosim {
  namespace tools {
    // increments the elements of a container based on the alphabets; containers and alphabets always have fixed length
    template <class T>
    class incrementor {
    public:
      // ctor taking one alphabet for each position; no alphabets are allowed (for subtensor), empty alphabets are not
      explicit incrementor(std::vector<std::vector<typename T::value_type>> const &__alphabets = {
                               get_uppercase_alphabet()})
          : _alphabets(__alphabets), _overflow(false) {
        for(auto alphabet : _alphabets) {
          if(alphabet.empty()) {
            throw std::invalid_argument("cannot increment with an empty alphabet");
          } // if
        } // for
      } // ctor

      // incrementing function; does not check at the beginning if all positions have allowed letters, only step by step
      T next(T __container, size_t increment = 1) {
        if(_alphabets.size() != __container.size()) { // an alphabet for each container position is required
          throw std::out_of_range("input length differs from alphabet count for fixed length input");
        } // if
        int pos(__container.size() - 1); // start with last element; use int because we're decreasing, and it can go < 0
        _overflow = (increment && pos < 0); // reset overflow, it could have been set before calling next()

        while(increment && !_overflow) { // as long as we still need to increment and we don't have an overflow
          std::vector<typename T::value_type> const &a(_alphabets[pos]); // get the alphabet
          typename std::vector<typename T::value_type>::const_iterator current_letter_itr(
              std::find(a.begin(), a.end(), __container[pos]));
          if(current_letter_itr == a.end()) { // check if current position has an allowed letter
            throw std::out_of_range("input contains letter not in alphabet");
          } // if
          size_t current_letter_pos(current_letter_itr - a.begin());
          __container[pos] = a[(current_letter_pos + increment) % a.size()]; // update current letter
          increment = (current_letter_pos + increment) / a.size(); // update increment
          --pos;
          _overflow = (increment && pos < 0);
        } // while

        return __container;
      } // next()

      // return if an overflow occured on the last call to next()
      bool overflow() const { return _overflow; }

      // get uppercase alphabet
      static std::vector<char> get_uppercase_alphabet() {
        static std::string letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        return std::vector<char>(letters.begin(), letters.end());
      } // get_uppercase_alphabet()

    private:
      std::vector<std::vector<typename T::value_type>> _alphabets; // the ordered alphabet for incrementing
      bool _overflow; // overflow state
    }; // class incrementor
  } // namespace tools
} // namespace biosim

#endif // tools_incrementor_h
