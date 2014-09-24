#include "che/ts.h"

namespace biosim {
  namespace che {
    // ctor from storage and identifier (also default ctor)
    ts::ts(std::string __storage, std::string __identifier) : _storage(__storage), _id(__identifier), _ps() {}
    // ctor from storage, identifier, and ps
    ts::ts(std::string __storage, std::string __identifier, ps __ps)
        : _storage(__storage), _id(__identifier), _ps(__ps) {}
    // get storage
    std::string const &ts::get_storage() const { return _storage; }
    // get identifier
    std::string const &ts::get_identifier() const { return _id; }
    // get length
    size_t ts::get_length() const { return _ps.size(); }
    // get ps
    ps const &ts::get_ps() const { return _ps; }
    // set ps
    void ts::set_ps(ps __ps) { _ps = __ps; }
  } // namespace che
} // namespace biosim