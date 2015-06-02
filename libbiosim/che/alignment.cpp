#include "che/alignment.h"

namespace biosim {
  namespace che {
    // ctor from a molecule
    alignment::alignment(molecule __molecule) : _molecules(1, __molecule) {}
    // ctor from a list of molecules
    alignment::alignment(std::list<molecule> __molecules) : _molecules(__molecules) {
      if(!_molecules.empty()) {
        for(auto m : _molecules) {
          if(m.get_length() != _molecules.begin()->get_length()) {
            throw std::invalid_argument("not all aligned molecules have the same length");
          } // if
        } // for
      } // if
    } // ctor
    // get alignment length
    size_t alignment::get_length() const { return _molecules.empty() ? 0 : _molecules.begin()->get_length(); }
    // get alignment depth, i.e. number of molecules
    size_t alignment::get_depth() const { return _molecules.size(); }
    // get aligned molecules
    std::list<molecule> const &alignment::get_molecules() const { return _molecules; }
  } // namespace che
} // namespace biosim