#include "che/structure.h"

namespace biosim {
  namespace che {
    // ctor from ps (also default ctor)
    structure::structure(ps __ps) : _ps(__ps), _storage(), _id(), _ss() {}
    // ctor from storage and identifier 
    structure::structure(std::string __storage, std::string __identifier)
        : _ps(), _storage(__storage), _id(__identifier), _ss() {}
    // ctor from storage, identifier, and ps
    structure::structure(std::string __storage, std::string __identifier, ps __ps)
        : _ps(__ps), _storage(__storage), _id(__identifier), _ss() {}
    // ctor from storage, identifier, ps, and ss
    structure::structure(std::string __storage, std::string __identifier, ps __ps, ss __ss)
        : _ps(__ps), _storage(__storage), _id(__identifier), _ss(__ss) {
      if(__ps.size() != __ss.get_length()) {
        throw std::invalid_argument("length of primary and secondary structure differs");
      } // if
    } // ctor
    // get ps
    ps const &structure::get_ps() const { return _ps; }
    // get length
    size_t structure::get_length() const { return get_ps().size(); }
    // get storage
    std::string const &structure::get_storage() const {
      static std::string storage_default;
      return _storage ? *_storage : storage_default;
    } // get_storage()
    // get identifier
    std::string const &structure::get_identifier() const { 
      static std::string id_default;
      return _id ? *_id : id_default;
    } // get_identifier()
    // get ss
    ss const &structure::get_ss() const {
      static ss ss_default;
      return _ss ? *_ss : ss_default;
    } // get_ss()
  } // namespace che
} // namespace biosim