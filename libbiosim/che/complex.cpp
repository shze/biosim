#include "che/complex.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // get list of all chain_ids
    std::list<std::string> complex::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto p : _molecules) {
        chain_ids.push_back(p.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // get molecule with given chain_id, throws out of range if no molecule with this chain_id exists
    molecule const &complex::get(std::string __chain_id) const { return _molecules.at(__chain_id); }
    // add molecule to complex, assign next available chain_id, and returns the chain_id
    std::string complex::add(molecule __molecule) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor inc;
      while(_molecules.find(chain_id) != _molecules.end()) {
        chain_id = inc.next(chain_id);
      }

      _molecules[chain_id] = __molecule;
      return chain_id;
    }
  } // namespace che
} // namespace biosim