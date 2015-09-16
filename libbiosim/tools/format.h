#ifndef tools_format_h
#define tools_format_h

#include <algorithm>
#include <sstream>

namespace biosim {
  namespace tools {
    // format a vector as string; should be replaced by a more general solution like:
    // http://louisdx.github.io/cxx-prettyprint/
    template <class T>
    std::string to_string(std::vector<T> const &__pos, std::string __separator = ", ") {
      std::stringstream s;
      if(!__pos.empty()) {
        std::for_each(__pos.begin(), --__pos.end(), [&](T const &__v) { s << __v << __separator; });
        s << __pos.back();
      } // if
      return s.str();
    } // to_string()
  } // namespace tools
} // namespace biosim

#endif // tools_format_h
