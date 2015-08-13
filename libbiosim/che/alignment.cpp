#include "che/alignment.h"

namespace biosim {
  namespace che {
    // ctor from a structure (takes begin and end from given structure)
    alignment::alignment(structure __s) : _structures(1, __s), _begins(1, 0), _ends(1, __s.get_length()) {}
    // ctor from a list of structures, begins, and ends
    alignment::alignment(std::vector<structure> __structures, std::vector<size_t> __begins, std::vector<size_t> __ends)
        : _structures(__structures), _begins(__begins), _ends(__ends) {
      if(!_structures.empty()) {
        for(auto s : _structures) {
          if(s.get_length() != _structures.begin()->get_length()) {
            throw std::invalid_argument("not all aligned structures have the same length");
          } // if
        } // for
      } // if
    } // ctor
    // get alignment length
    size_t alignment::get_length() const { return _structures.empty() ? 0 : _structures.begin()->get_length(); }
    // get alignment depth, i.e. number of structures
    size_t alignment::get_depth() const { return _structures.size(); }
    // get aligned structures
    std::vector<structure> const &alignment::get_structures() const { return _structures; }
    // get begins (1-based) of the alignment in the complete structures
    std::vector<size_t> const &alignment::get_begins() const { return _begins; }
    // get ends (1-based) of the alignment in the complete structures
    std::vector<size_t> const &alignment::get_ends() const { return _ends; }
    // get cc at specific position and depth
    cc const &alignment::get_cc(size_t __position, size_t __depth) const {
      return _structures.at(__depth).get_ps().at(__position);
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