#include "che/complex.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // add molecule to complex, assign next available chain_id
    void complex::add(molecule __molecule) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor inc;
      while(_molecules.find(chain_id) != _molecules.end()) {
        chain_id = inc.next(chain_id);
      }

      _molecules[chain_id] = __molecule;
    }
  } // namespace che
} // namespace biosim