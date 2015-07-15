#include "che/alignment.h"

namespace biosim {
  namespace che {
    // ctor from a molecule
    alignment::alignment(molecule __molecule) : _molecules(1, __molecule) {}
    // ctor from a list of molecules
    alignment::alignment(std::vector<molecule> __molecules) : _molecules(__molecules) {
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
    std::vector<molecule> const &alignment::get_molecules() const { return _molecules; }
    // get cc at specific position and depth
    cc const &alignment::get_cc(size_t __position, size_t __depth) const {
      return _molecules.at(__depth).get_ps().at(__position);
    } // get_cc()

    // ctor from alignment and score
    scored_alignment::scored_alignment(alignment __alignment, double __score)
        : _alignment(__alignment), _score(__score) {}
    // get the alignment
    alignment const &scored_alignment::get_alignment() const { return _alignment; }
    // get the score
    double scored_alignment::get_score() const { return _score; }
  } // namespace che
} // namespace biosim