#include "che/atom.h"

namespace biosim {
  namespace che {
    // ctor from identifier and position; default ctor
    atom::atom(std::string __identifier, math::point __position) : _id(__identifier), _pos(__position) {}
    // get identifier
    std::string const &atom::get_identifier() const { return _id; }
    // get position
    boost::optional<math::point> const &atom::get_position() const { return _pos; }
    // returns true if this identifier is less than the rhs identifier
    bool atom::operator<(atom const &__rhs) const { return this->get_identifier() < __rhs.get_identifier(); }
  } // namespace che
} // namespace biosim