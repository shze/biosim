#ifndef che_ps_h
#define che_ps_h

#include <iostream>
#include <vector>
#include "che/cc.h"

namespace biosim {
  namespace che {
    // primary structure
    class ps : public std::vector<cc> {};

    // output operator for ps
    inline std::ostream &operator<<(std::ostream &__out, ps const &__ps) {
      for(auto const &component : __ps) { // component is of type cc
        __out << component.get_identifier_char();
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_ps_h
