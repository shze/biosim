#ifndef che_alignment_h
#define che_alignment_h

#include "che/molecule.h"

namespace biosim {
  namespace che {
    // alignment of molecules, consensus/profile;
    // design: unified alignment class for sequence and structure alignment, cleaner than replicating everything.
    class alignment {
    public:
      // ctor from a molecule
      explicit alignment(molecule __molecule);
      // ctor from a list of molecules
      explicit alignment(std::vector<molecule> __molecules);
      // get alignment length
      size_t get_length() const;
      // get alignment depth, i.e. number of molecules
      size_t get_depth() const;
      // get aligned molecules
      std::vector<molecule> const &get_molecules() const;
      // get cc at specific position and depth
      cc const &get_cc(size_t __position, size_t __depth) const;

    private:
      std::vector<molecule> _molecules; // list of aligned molecules
    }; // class alignment

    // output operator for alignment
    inline std::ostream &operator<<(std::ostream &__out, alignment const &__a) {
      __out << "Alignment (depth=" << __a.get_depth() << ", size=" << __a.get_length() << ")\n";
      for(auto const &m : __a.get_molecules()) {
        __out << m;
      } // for
      return __out;
    } // operator<<()

    // container scoring an alignment together with its score; design: simple container
    class scored_alignment {
    public:
      // ctor from alignment and score
      scored_alignment(alignment __alignment, double __score);
      // get the alignment
      alignment const &get_alignment() const;
      // get the score
      double get_score() const;

    private:
      alignment _alignment; // alignment
      double _score; // score
    }; // class scored_alignment

    // output operator for scored_alignment
    inline std::ostream &operator<<(std::ostream &__out, scored_alignment const &__a) {
      __out << "Alignment (depth=" << __a.get_alignment().get_depth() << ", size=" << __a.get_alignment().get_length()
            << ", score=" << __a.get_score() << ")\n";
      for(auto const &m : __a.get_alignment().get_molecules()) {
        __out << m;
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace

#endif // che_alignment_h
