#include "che/alignment.h"

namespace biosim {
  namespace che {
    // ctor from a molecule
    alignment::alignment(molecule __molecule) : _molecules(1, __molecule) {}
    // ctor from a list of molecules
    alignment::alignment(std::list<molecule> __molecules) : _molecules(__molecules) {}
    // get aligned molecules
    std::list<molecule> const &alignment::get_molecules() { return _molecules; }
  } // namespace che
} // namespace biosim