#ifndef che_atom_h
#define che_atom_h

#include <string>
#include "math/point.h"
#include <boost/optional.hpp>

namespace biosim {
  namespace che {
    // stores an atom
    class atom {
    public:
      // ctor from identifier and position; default ctor
      explicit atom(std::string __identifier = std::string(), math::point __position = math::point());
      // get identifier
      std::string const &get_identifier() const;
      // returns true if this identifier is less than the rhs identifier
      bool operator<(atom const &__rhs) const;

    private:
      std::string _id; // identifier
      boost::optional<math::point> _pos; // optional position
    }; // class atom
  } // namespace che
} // namespace biosim

#endif // che_atom_h
