#include "che/assembly.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // default ctor
    assembly::assembly() : _molecules() {}
    // ctor from molecule
    assembly::assembly(molecule __m) : _molecules() { add(__m); }
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
    // get molecule with given chain_id, throws out of range if no ts with this chain_id exists
    molecule const &assembly::get_molecule(std::string const &__chain_id) const { return _molecules.at(__chain_id); }
    // add molecule to assembly, assign next available chain_id, and returns the chain_id
    std::string assembly::add(molecule __m) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor<std::string> inc;
      while(_molecules.find(chain_id) != _molecules.end()) {
        chain_id = inc.next(chain_id);
      } // while

      set(chain_id, __m);
      return chain_id;
    } // add()
    // set a molecule and ss to have a specific chain_id
    void assembly::set(std::string const &__chain_id, molecule __m) { _molecules[__chain_id] = __m; }
  } // namespace che
} // namespace biosim