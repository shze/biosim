#include "che/assembly.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // default ctor
    assembly::assembly() : _molecules(), _ss() {}
    // ctor from molecule and ss
    assembly::assembly(molecule __m, ss __ss) : _molecules(), _ss() { add(__m, __ss); }
    // get list of all chain_ids
    std::list<std::string> assembly::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto m : _molecules) {
        chain_ids.push_back(m.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // if the molecule with the given chain_id exists
    bool assembly::has_molecule(std::string const &__chain_id) const {
      return _molecules.find(__chain_id) != _molecules.end();
    }
    // if the ss with the given chain_id exists
    bool assembly::has_ss(std::string const &__chain_id) const { return _ss.find(__chain_id) != _ss.end(); }
    // get molecule with given chain_id, throws out of range if no ts with this chain_id exists
    molecule const &assembly::get_molecule(std::string const &__chain_id) const { return _molecules.at(__chain_id); }
    // get ss with given chain_id, throws out of range if no ts with this chain_id exists
    ss const &assembly::get_ss(std::string const &__chain_id) const { return _ss.at(__chain_id); }
    // add molecule to assembly, assign next available chain_id, and returns the chain_id
    std::string assembly::add(molecule __m) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor<std::string> inc;
      while(_molecules.find(chain_id) != _molecules.end()) {
        chain_id = inc.next(chain_id);
      } // while

      _molecules.insert(std::pair<std::string, molecule>(chain_id, __m));
      return chain_id;
    } // add()
    // add corresponding molecule and ss to assembly; assign next available chain_id, and returns the chain_id
    std::string assembly::add(molecule __m, ss __ss) {
      std::string chain_id = add(__m);
      _ss.insert(std::pair<std::string, ss>(chain_id, __ss));
      return chain_id;
    } // add()
    // set a molecule and ss to have a specific chain_id
    void assembly::set(std::string const &__chain_id, molecule __m, ss __ss) {
      _molecules[__chain_id] = __m;
      // for ss use erase() and insert(), because operator[] requires a default ctor which ss does not have
      _ss.erase(__chain_id);
      _ss.insert(std::pair<std::string, ss>(__chain_id, __ss));
    } // set()
  } // namespace che
} // namespace biosim