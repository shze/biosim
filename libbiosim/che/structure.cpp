#include "che/structure.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    structure::structure(std::string __storage, std::string __identifier)
        : _storage(__storage), _id(__identifier), _ps(), _ss() {}
    // ctor from storage, identifier, and ps
    structure::structure(std::string __storage, std::string __identifier, ps __ps)
        : _storage(__storage), _id(__identifier), _ps(__ps), _ss() {}
    // ctor from storage, identifier, ps, and ss
    structure::structure(std::string __storage, std::string __identifier, ps __ps, ss __ss)
        : _storage(__storage), _id(__identifier), _ps(__ps), _ss(__ss) {
      if(_ps.size() != _ss.get_length()) {
        throw std::invalid_argument("length of primary and secondary structure differs");
      } // if
    } // ctor
    // get storage
    std::string const &structure::get_storage() const { return _storage; }
    // get identifier
    std::string const &structure::get_identifier() const { return _id; }
    // get length
    size_t structure::get_length() const { return _ps.size(); }
    // get ps
    ps const &structure::get_ps() const { return _ps; }
    // get ss
    ss const &structure::get_ss() const { return _ss; }
  } // namespace che
} // namespace biosim