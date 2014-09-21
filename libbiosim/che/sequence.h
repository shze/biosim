#ifndef che_sequence_h
#define che_sequence_h

#include <iostream>
#include <vector>
#include "che/cc.h"

namespace biosim {
  namespace che {
    // sequence
    template <typename T>
    class sequence : public std::vector<T> {};

    // output operator for sequence
    template <typename T>
    inline std::ostream &operator<<(std::ostream &__out, sequence<T> const &__seq) {
      for(auto const &component : __seq) { // component is of type T
        __out << component.get_identifier_char(); // T is required to have T::get_identifier_char()
      } // for
      return __out;
    } // operator<<()

    // primary structure
    using ps = sequence<cc>;
  } // namespace che
} // namespace biosim

#endif // che_sequence_h
