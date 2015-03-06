#include "che/molecule.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    molecule::molecule(std::string __storage, std::string __identifier)
        : _storage(__storage), _id(__identifier), _ps() {}
    // ctor from storage, identifier, and ps
    molecule::molecule(std::string __storage, std::string __identifier, ps __ps)
        : _storage(__storage), _id(__identifier), _ps(__ps) {}
    // get storage
    std::string const &molecule::get_storage() const { return _storage; }
    // get identifier
    std::string const &molecule::get_identifier() const { return _id; }
    // get length
    size_t molecule::get_length() const { return _ps.size(); }
    // get ps
    ps const &molecule::get_ps() const { return _ps; }
    // set ps
    void molecule::set_ps(ps __ps) { _ps = __ps; }
  } // namespace che
} // namespace biosim