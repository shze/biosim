#include "che/alignment.h"

namespace biosim {
  namespace che {
    // ctor from a molecule (takes begin and end from given molecule)
    alignment::alignment(molecule __molecule)
        : _molecules(1, __molecule), _begins(1, 0), _ends(1, __molecule.get_length()) {}
    // ctor from a list of molecules, begins, and ends
    alignment::alignment(std::vector<molecule> __molecules, std::vector<size_t> __begins, std::vector<size_t> __ends)
        : _molecules(__molecules), _begins(__begins), _ends(__ends) {
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
    // get begins (1-based) of the alignment in the complete molecules
    std::vector<size_t> const &alignment::get_begins() const { return _begins; }
    // get ends (1-based) of the alignment in the complete molecules
    std::vector<size_t> const &alignment::get_ends() const { return _ends; }
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