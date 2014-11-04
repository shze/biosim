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
    private:
      T _type; // type of all elements from [begin .. end] (closed interval)

    public:
      // ctor taking interval describing the dimension and a T describing the type; no default ctor
      sequence_interval(size_t __min, size_t __max, T __type) : interval(__min, __max), _type(__type) {}

      // get the T describing the type
      T const &get_type() const { return _type; }

      // returns if interval1 is less than i.e. before (starts further left, at a smaller number) interval2
      static bool less_min_max(sequence_interval<T> const &__interval1, sequence_interval<T> const &__interval2) {
        return math::interval<size_t>::less_min_max(__interval1, __interval2);
      } // less_min_max()
      // returns if the length of interval1 is less than the length of interval2
      static bool less_length(sequence_interval<T> const &__interval1, sequence_interval<T> const &__interval2) {
        return math::interval<size_t>::less_length(__interval1, __interval2);
      } // less_length()
      // returns if interval1 ends before before interval2 (ends further left, at a smaller number)
      static bool less_max(sequence_interval<T> const &__interval1, sequence_interval<T> const &__interval2) {
        return math::interval<size_t>::less_max(__interval1, __interval2);
      } // less_max()

      // returns if both intervals have same min and max positions
      static bool equal_min_max(sequence_interval<T> const &__interval1, sequence_interval<T> const &__interval2) {
        return math::interval<size_t>::equal_min_max(__interval1, __interval2);
      } // equal_min_max()
      // returns if both intervals have same length
      static bool equal_length(sequence_interval<T> const &__interval1, sequence_interval<T> const &__interval2) {
        return math::interval<size_t>::equal_length(__interval1, __interval2);
      } // equal_length()

      // test if this interval and interval2 overlap, depending on type; overwrites interval::overlaps()
      bool overlaps(sequence_interval<T> const &__interval2, bool __all_types_overlap = true) const {
        return __interval2.interval::overlaps(*this) && (__all_types_overlap || get_type() == __interval2.get_type());
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

      // convert a set of sequence_intervals into a sequence
      static che::sequence<T> to_sequence(std::set<che::sequence_interval<T>> const &__intervals,
                                          size_t const &__length) {
        // find unknown; used to fill ss_symbols
        T unknown(T::specificity_unknown);

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
