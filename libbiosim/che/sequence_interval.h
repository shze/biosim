#ifndef che_sequence_interval_h
#define che_sequence_interval_h

#include "math/interval.h"
#include <set>
#include "che/cchb_dssp.h"

namespace biosim {
  namespace che {
    // defines an interval of a sequence with a closed interval denoting the consecutive identical elements
    template <typename T>
    class sequence_interval : public math::interval<size_t> { // derived from closed interval
    private:
      T _type; // type of all elements from [begin .. end] (closed interval)

    public:
      // ctor taking interval describing the dimension and a T describing the type; no default ctor
      sequence_interval(size_t __min, size_t __max, T __type) : interval(__min, __max), _type(__type) {}

      // get the T describing the type
      T const &get_type() const { return _type; }

      // test if this interval and interval2 overlap, depending on type; overwrites interval::overlaps()
      bool overlaps(sequence_interval<T> const &__interval2, bool __all_types_overlap = false) const {
        return __interval2.interval::overlaps(*this) && (__all_types_overlap || get_type() == __interval2.get_type());
      } // overlaps()
      // return the intervals from the iterator interval begin to end that overlap with the given interval; considers
      // all intervals as overlapping if __all_types_overlap=true, otherwise only intervals of same type.
      // for generic iterator: http://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator
      template <typename I>
      std::set<sequence_interval<T>> get_overlapping_subset(I __begin, I __end, bool __all_types_overlap = false) {
        std::set<che::sequence_interval<T>> overlapping_set;
        for(; __begin != __end; ++__begin) {
          if(__begin->overlaps(*this, __all_types_overlap)) {
            overlapping_set.insert(*__begin);
          }
        } // for
        return overlapping_set;
      } // overlapping_subset()
      // check if the given interval overlaps with any intervals in the sequence_interval iterator range begin to end
      template <typename I>
      bool overlaps(I __begin, I __end, bool __all_types_overlap = false) {
        return !get_overlapping_subset(__begin, __end, __all_types_overlap).empty();
      } // overlaps()
    }; // class sequence_interval

    // output operator for sequence_interval
    template <typename T>
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
