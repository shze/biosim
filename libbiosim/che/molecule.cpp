#include "che/molecule.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    molecule::molecule(std::string const &__storage, std::string const &__identifier)
        : _storage(__storage), _id(__identifier), _ps(), _ts() {}
    // get storage
    std::string const &molecule::get_storage() const { return _storage; }
    // get identifier
    std::string const &molecule::get_identifier() const { return _id; }
    // get ps
    ps const &molecule::get_ps() const { return _ps; }
    // set ps
    void molecule::set_ps(ps __ps) { _ps = __ps; }
  } // namespace che
} // namespace biosim