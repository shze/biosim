#ifndef che_complex_h
#define che_complex_h

#include "che/molecule.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // complex of multiple molecules with additional information
    class complex {
    public:
      // add molecule to complex, assign next available chain_id
      void add(molecule __molecule);

    private:
      std::map<std::string, molecule> _molecules; // maps chain_id -> molecule
      ss _ss; // secondary structure
    }; // class complex
  } // namespace che
} // namespace biosim

#endif // che_complex_h
