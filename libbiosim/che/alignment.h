#ifndef che_alignment_h
#define che_alignment_h

#include "che/structure.h"

namespace biosim {
  namespace che {
    // alignment of structures, consensus/profile;
    // design: unified alignment class for sequence and structure alignment, cleaner than replicating everything.
    class alignment {
    public:
      // ctor from a structure (takes begin and end from given structure)
      explicit alignment(structure __s);
      // ctor from a list of structures, begins, and ends
      alignment(std::vector<structure> __structures, std::vector<size_t> __begins, std::vector<size_t> __ends);
      // get alignment length
      size_t get_length() const;
      // get alignment depth, i.e. number of structures
      size_t get_depth() const;
      // get aligned structures
      std::vector<structure> const &get_structures() const;
      // get begins (1-based) of the alignment in the complete structures
      std::vector<size_t> const &get_begins() const;
      // get ends (1-based) of the alignment in the complete structures
      std::vector<size_t> const &get_ends() const;
      // get cc at specific position and depth
      cc const &get_cc(size_t __position, size_t __depth) const;

    private:
      std::vector<structure> _structures; // vector of aligned structures
      std::vector<size_t> _begins; // vector of begins
      std::vector<size_t> _ends; // vector of ends
    }; // class alignment

    // output operator for alignment
    inline std::ostream &operator<<(std::ostream &__out, alignment const &__a) {
      __out << "Alignment (depth=" << __a.get_depth() << ", size=" << __a.get_length() << ")\n";
      for(size_t pos(0); pos < __a.get_depth(); ++pos) { // use pos to access structures, begins, and ends
        che::structure const &p(__a.get_structures()[pos]);
        std::string id(p.get_identifier());
        id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
        __out << p.get_storage() << ": " << id << " length=" << p.get_length()
              << " subsequence=" << __a.get_begins()[pos] << ".." << __a.get_ends()[pos] << "\n" << p.get_ps() << "\n";
        if(p.get_ss().defined()) {
          __out << p.get_ss();
        } // if
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
            << ", score=" << __a.get_score()
            << ", normalized_score=" << (__a.get_score() / __a.get_alignment().get_length()) << ")\n";
      for(size_t pos(0); pos < __a.get_alignment().get_depth(); ++pos) {
        che::structure p(__a.get_alignment().get_structures()[pos]);
        std::string id(p.get_identifier());
        id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
        __out << p.get_storage() << ": " << id << " length=" << p.get_length()
              << " subsequence=" << __a.get_alignment().get_begins()[pos] << ".." << __a.get_alignment().get_ends()[pos]
              << "\n" << p.get_ps() << "\n";
        if(p.get_ss().defined()) {
          __out << p.get_ss();
        } // if
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace

#endif // che_alignment_h
