#ifndef tools_mapper_h
#define tools_mapper_h

#include <stdexcept>
#include <iterator>
#include "tools/format.h"

namespace biosim {
  namespace tools {
    // map size_t to alphabets; assumes first letter of an alphabet denotes 0; fills all digits (one per alphabet)
    template <class T>
    class mapper {
    public:
      using alphabet = std::vector<typename T::value_type>; // simplify naming of alphabet
      using alphabet_container = std::vector<alphabet>; // simplify naming of alphabet_container

      // ctor from alphabets and optional big_endian bool
      explicit mapper(alphabet_container __alphabets, bool __big_endian = false)
          : _alphabets(__alphabets),
            _big_endian(__big_endian),
            _digit_factor(calculate_digit_factor()), // needs _alphabets and _big_endian
            _max(calculate_max_container()) {} // needs _alphabets and _digit_factor
      // get alphabets
      alphabet_container const &get_alphabets() const { return _alphabets; }
      // return if big endian
      bool big_endian() const { return _big_endian; }
      // get minimum representable number
      size_t get_min() const { return 0; }
      // get maximum representable number
      size_t get_max() const { return _max; }
      // decode a container to its corresponding number
      size_t decode(T const &__container) const {
        if(__container.size() != _alphabets.size()) {
          throw std::invalid_argument("container size differs from alphabet size");
        } // if

        size_t value(0);
        for(size_t dim(0); dim < __container.size(); ++dim) {
          typename alphabet::const_iterator current_letter_itr(
              std::find(_alphabets[dim].begin(), _alphabets[dim].end(), __container[dim]));
          if(current_letter_itr == _alphabets[dim].end()) {
            throw std::out_of_range("input contains letter not in alphabet");
          } // if
          value += std::distance(_alphabets[dim].begin(), current_letter_itr) * _digit_factor[dim];
        } // for
        return value;
      } // decode()
      // encode a number into its corresponding container
      T encode(size_t __value) const {
        T container(_alphabets.size(), typename T::value_type()); // create container with full length
        __value = std::min(__value, get_max()); // ensure value is at most get_max()

        int step(_big_endian ? -1 : 1);
        int start_dim(_big_endian ? _alphabets.size() - 1 : 0), end_dim(_big_endian ? -1 : _alphabets.size());
        for(int dim(start_dim); dim != end_dim; dim += step) { // loop until we have no dimensions left to fill
          alphabet const &current_alphabet(_alphabets[dim]);
          size_t current_alphabet_size = current_alphabet.size();
          size_t next_value(__value / current_alphabet_size);
          container[dim] = current_alphabet[__value - next_value * current_alphabet_size];
          __value = next_value;
        } // for

        return container;
      } // encode()

    private:
      // calculate the factor for each digit; needs _alphabets and _big_endian to be set
      std::vector<size_t> calculate_digit_factor() const {
        std::vector<size_t> digit_factor(_alphabets.size());
        int step(_big_endian ? -1 : 1); // from lowest to highest digit
        int start_dim(_big_endian ? _alphabets.size() - 1 : 0), end_dim(_big_endian ? -1 : _alphabets.size());
        for(int dim(start_dim); dim != end_dim; dim += step) {
          digit_factor[dim] = dim == start_dim ? 1 : (digit_factor[dim - step] * _alphabets[dim - step].size());
        } // for
        return digit_factor;
      } // calculate_digit_factor()
      // calculate maximum representable number; needs _alphabets and _digit_factor (through decode) to be set
      size_t calculate_max_container() const {
        for(auto const &alphabet : _alphabets) {
          if(alphabet.empty()) { // make sure a given alphabet is not empty (otherwise mapping will be tricky)
            throw std::invalid_argument("cannot map to an empty alphabet");
          } // if
        } // for

        T max_container(_alphabets.size(), typename T::value_type());
        for(size_t dim(0); dim < _alphabets.size(); ++dim) {
          max_container[dim] = _alphabets[dim].back();
        } // for
        return decode(max_container);
      } // calculate_max()

      alphabet_container _alphabets; // ordered alphabets for mapping
      bool _big_endian; // if mapping is big endian ordered
      std::vector<size_t> _digit_factor; // factor for each digit
      size_t _max; // maximum representable number
    }; // class mapper

    template <class T>
    inline std::ostream &operator<<(std::ostream &__out, mapper<T> const &__m) {
      std::stringstream s;
      typename mapper<T>::alphabet_container const &alphabets(__m.get_alphabets());
      if(!alphabets.empty()) {
        std::for_each(alphabets.begin(), --alphabets.end(),
                      [&](typename mapper<T>::alphabet const &__a) { s << "{" << tools::to_string(__a) << "}, "; });
        s << "{" << tools::to_string(alphabets.back()) << "}";
      } // if

      __out << "mapper: alphabets={" << s.str() << "}, big_endian=" << __m.big_endian() << ", min=" << __m.get_min()
            << ", max=" << __m.get_max();
      return __out;
    } // operator<<()
  } // namespace tools
} // namespace biosim

#endif // tools_mapper_h
