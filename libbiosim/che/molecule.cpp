#include "che/molecule.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    molecule::molecule(std::string const &__storage, std::string const &__identifier)
        : _storage(__storage), _id(__identifier), _ps(), _ts() {}
    // set ps
    void molecule::set_ps(ps __ps) { _ps = __ps; }
  } // namespace che
} // namespace biosim