#ifndef che_alignment_h
#define che_alignment_h

#include "che/molecule.h"
#include <list>

namespace biosim {
  namespace che {
    // alignment of molecules, consensus/profile;
    // design: unified alignment class for sequence and structure alignment, cleaner than replicating everything.
    class alignment {
    public:
      // ctor from a molecule
      explicit alignment(molecule __molecule);
      // ctor from a list of molecules
      explicit alignment(std::list<molecule> __molecules);
      // get alignment length
      size_t get_length() const;
      // get alignment depth, i.e. number of molecules
      size_t get_depth() const;
      // get aligned molecules
      std::list<molecule> const &get_molecules() const;

    private:
      std::list<molecule> _molecules; // list of aligned molecules
    }; // class alignment
  } // namespace che
} // namespace

#endif // che_alignment_h
