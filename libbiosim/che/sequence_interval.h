#ifndef che_sequence_interval_h
#define che_sequence_interval_h

#include "math/interval.h"
#include <set>
#include "che/cchb_dssp.h"
#include "che/sequence.h"

namespace biosim {
  namespace che {
    // defines an interval of a sequence with a closed interval denoting the consecutive identical elements
    template <class T>
    class sequence_interval : public math::interval<size_t> { // derived from closed interval
    public:
      // ctor taking interval describing the dimension and a T describing the type; no default ctor
      sequence_interval(size_t __min, size_t __max, T __type) : interval(__min, __max), _type(__type) {}

      // get the T describing the type
      T const &get_type() const { return _type; }

      // test if this interval and interval2 overlap, depending on type; overwrites interval::overlaps()
      bool overlaps(sequence_interval<T> const &__interval2, bool __all_types_overlap = true) const {
        return (__all_types_overlap || get_type() == __interval2.get_type()) && __interval2.interval::overlaps(*this);
      } // overlaps()
      // return the intervals from the iterator interval begin to end that overlap with the given interval; considers
      // all intervals as overlapping if __all_types_overlap=true, otherwise only intervals of same type.
      // for generic iterator: http://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator
      template <class I>
      std::set<sequence_interval<T>> get_overlapping_subset(I __begin, I __end, bool __all_types_overlap = true) {
        std::set<che::sequence_interval<T>> overlapping_set;
        for(; __begin != __end; ++__begin) {
          if(__begin->overlaps(*this, __all_types_overlap)) {
            overlapping_set.insert(*__begin);
          } // if
        } // for
        return overlapping_set;
      } // overlapping_subset()
      // check if the given interval overlaps with any intervals in the sequence_interval iterator range begin to end
      template <class I>
      bool overlaps(I __begin, I __end, bool __all_types_overlap = true) {
        return !get_overlapping_subset(__begin, __end, __all_types_overlap).empty();
      } // overlaps()

      // convert a set of sequence intervals into a sequence
      static che::sequence<T> to_sequence(std::set<che::sequence_interval<T>> const &__intervals,
                                          size_t const &__length) {
        // find unknown; used to fill ss_symbols
        T unknown(T::specificity_type::unknown);

        che::sequence<T> seq;
        for(auto const &this_interval : __intervals) {
          DEBUG << "Extend ss_symbols with interval " << this_interval;
          // intervals has 0-based positions, extend_symbols() takes length
          // extend until one before min(), i.e. position 0..min()-1, i.e. length of min()
          seq.resize(this_interval.get_min(), unknown);
          // extend until max(), i.e. positions min()..max(), i.e. length of max()+1
          seq.resize(this_interval.get_max() + 1, this_interval.get_type());
        } // for

        // fills up until the end defined by __length, if __length > ss_symbols.size()
        if(__length > seq.size()) {
          DEBUG << "Extend ss_symbols to length " << __length;
          seq.resize(__length, unknown);
        } // if

        return seq;
      } // to_sequence

      // convert a sequence into a set of sequence intervals
      static std::set<che::sequence_interval<T>> to_sequence_intervals(che::sequence<T> const &__sequence) {
        // return an empty set right away and avoid checking for empty sets in the logic below
        if(__sequence.size() == 0) {
          return std::set<che::sequence_interval<T>>();
        } // if

        // create intervals for all types
        std::list<sequence_interval<T>> sequence_interval_list; // first create list to extend the last intervals
        for(size_t pos(0); pos < __sequence.size(); ++pos) {
          if(sequence_interval_list.empty()) {
            sequence_interval<T> new_sequence_interval(pos, pos, __sequence[pos]);
            DEBUG << "sequence_interval_list empty, adding new sequence_interval " << new_sequence_interval;
            sequence_interval_list.push_back(new_sequence_interval);
            continue; // to avoid the else and intendation
          } // if

          sequence_interval<T> &last_sequence_interval(sequence_interval_list.back()); // reference is important!
          if(last_sequence_interval.get_type().get_identifier() == __sequence[pos].get_identifier()) {
            last_sequence_interval.set_max(pos);
          } // if
          else {
            // print this once to tell how much we extended when we find something different, instead of printing for
            // every extend in the if part above
            DEBUG << "identifier was equal, extended sequence_interval " << last_sequence_interval;

            sequence_interval<T> new_sequence_interval(pos, pos, __sequence[pos]);
            DEBUG << "identifier is different, adding new sequence_interval " << new_sequence_interval;
            sequence_interval_list.push_back(new_sequence_interval);
          } // else
        } // for
        // print this once instead of for every single extend
        DEBUG << "identifier was equal, extended sequence_interval " << sequence_interval_list.back();

        // insert all sequence_intervals except ones with unknown type into sequence_interval_set
        std::set<che::sequence_interval<T>> sequence_interval_set;
        T unknown(T::specificity_type::unknown);
        for(auto const &this_sequence_interval : sequence_interval_list) {
          if(this_sequence_interval.get_type().get_identifier() != unknown.get_identifier()) { // do NOT insert unknown
            sequence_interval_set.insert(this_sequence_interval);
          } // if
        } // for

        return sequence_interval_set;
      } // to_sequence_intervals()

    private:
      T _type; // type of all elements from [begin .. end] (closed interval)
    }; // class sequence_interval

    // output operator for sequence_interval
    template <class T>
    inline std::ostream &operator<<(std::ostream &__out, sequence_interval<T> const &__sequence_interval) {
      __out << "[" << __sequence_interval.get_min() << ", " << __sequence_interval.get_max()
            << "; epsilon=" << sequence_interval<T>::get_epsilon() << "]->"
            << __sequence_interval.get_type().get_identifier_char();
      return __out;
    } // operator<<()

    using cchb_dssp_interval = sequence_interval<cchb_dssp>;
  } // namespace che
} // namespace biosim

#endif // che_sequence_interval_h
