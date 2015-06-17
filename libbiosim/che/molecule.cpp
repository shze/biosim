#include "che/molecule.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    molecule::molecule(std::string __storage, std::string __identifier)
        : _storage(__storage), _id(__identifier), _ps(), _ss() {}
    // ctor from storage, identifier, and ps
    molecule::molecule(std::string __storage, std::string __identifier, ps __ps)
        : _storage(__storage), _id(__identifier), _ps(__ps), _ss() {}
    // ctor from storage, identifier, ps, and ss
    molecule::molecule(std::string __storage, std::string __identifier, ps __ps, ss __ss)
        : _storage(__storage), _id(__identifier), _ps(__ps), _ss(__ss) {
      if(_ps.size() != _ss.get_length()) {
        throw std::invalid_argument("length of primary and secondary structure differs");
      } // if
    } // ctor
    // get storage
    std::string const &molecule::get_storage() const { return _storage; }
    // get identifier
    std::string const &molecule::get_identifier() const { return _id; }
    // get length
    size_t molecule::get_length() const { return _ps.size(); }
    // get ps
    ps const &molecule::get_ps() const { return _ps; }
    // get ss
    ss const &molecule::get_ss() const { return _ss; }
  } // namespace che
} // namespace biosim