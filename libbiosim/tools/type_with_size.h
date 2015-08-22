#ifndef tools_type_with_size_h
#define tools_type_with_size_h

#include <boost/math/special_functions/next.hpp>
#include "tools/log.h"

namespace biosim {
  namespace tools {
    // template for compile time conversion from size to type; currently only implemented for uint;
    // default is void, which cannot be used, because it is an incomplete type.
    template <size_t size>
    struct type_with_size {
      using uint = void;
    }; // struct type_with_size

    // specialized template for data types for specified size
    template <>
    struct type_with_size<4> {
      using uint = uint32_t;
    }; // struct type_with_size

    // specialized template for data types for specified size
    template <>
    struct type_with_size<8> {
      using uint = uint64_t;
    }; // struct type_with_size
  } // namespace tools
} // namespace biosim

#endif // tools_type_with_size_h
