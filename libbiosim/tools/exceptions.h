#ifndef tools_exceptions_h
#define tools_exceptions_h

#include <boost/exception/all.hpp>
#include <boost/iterator/iterator_concepts.hpp>

namespace biosim {
  namespace tools {
    using exception_desc = boost::error_info<struct tag_exception_desc, std::string>;
  } // namespace tools
} // namespace biosim

#endif // tools_exceptions_h
