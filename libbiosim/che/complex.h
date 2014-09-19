#ifndef che_complex_h
#define che_complex_h

#include "che/molecule.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // complex of multiple molecules with additional information
    class complex {
    public:
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // get molecule with given chain_id, throws out of range if no molecule with this chain_id exists
      molecule const &get(std::string __chain_id) const;
      // add molecule to complex, assign next available chain_id, and returns the chain_id
      std::string add(molecule __molecule);

    private:
      std::map<std::string, molecule> _molecules; // maps chain_id -> molecule
      ss _ss; // secondary structure
    }; // class complex
  } // namespace che
} // namespace biosim

#endif // che_complex_h
