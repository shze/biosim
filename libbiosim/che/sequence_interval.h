#ifndef che_sequence_interval_h
#define che_sequence_interval_h

#include "math/interval.h"
#include <set>
#include "che/cchb.h"

namespace biosim {
  namespace che {
    // defined a interval of a sequence with a number of consecutive identical elements
    template <typename T>
    class sequence_interval {
    private:
      math::interval<size_t> _dimension; // begin and end (closed interval, includes begin and end)
      T _type; // type of all elements from [begin .. end]

    public:
      // ctor taking interval describing the dimension and a T describing the type; no default ctor
      sequence_interval(math::interval<size_t> __dimension, T __type) : _dimension(__dimension), _type(__type) {}

      // get the interval describing the dimension
      math::interval<size_t> const &get_dimension() const { return _dimension; }
      // get the T describing the type
      T const &get_type() const { return _type; }

      // return the intervals from the iterator interval begin to end that overlap with the given interval; considers
      // all intervals as overlapping if __all_types_overlap=true, otherwise only intervals of same type.
      // for generic iterator: http://stackoverflow.com/questions/5054087/declare-a-function-accepting-generic-iterator
      template <typename I>
      static std::set<sequence_interval<T>> get_overlapping_subset(sequence_interval const &__interval, I __begin,
                                                                   I __end, bool __all_types_overlap = false) {
        std::set<che::sequence_interval<T>> overlapping_set;
        for(; __begin != __end; ++__begin) {
          if(__begin->get_dimension().overlaps(__interval.get_dimension())) {
            // intervals overlap, if their dimensions overlap and either __all_types_overlap or they have the same type
            if(__all_types_overlap || __begin->get_type() == __interval.get_type()) {
              overlapping_set.insert(*__begin);
            } // if
          } // if
        } // for
        return overlapping_set;
      } // overlapping_subset()
      // check if the given interval overlaps with any intervals in the sequence_interval iterator range begin to end
      template <typename I>
      static bool overlaps(sequence_interval const &__interval, I __begin, I __end, bool __all_types_overlap = false) {
        return !get_overlapping_subset(__interval, __begin, __end, __all_types_overlap).empty();
      } // overlap()

      // compares two sequence intervals by comparing dimension, with the default interval min/max, and type
      static bool less_dimension_min_max(sequence_interval const &__interval1, sequence_interval const &__interval2) {
        return math::interval<size_t>::less_min_max(__interval1.get_dimension(), __interval2.get_dimension()) ||
               (!(math::interval<size_t>::less_min_max(__interval2.get_dimension(), __interval1.get_dimension())) &&
                __interval1.get_type() < __interval2.get_type());
      } // less_dimension_min_max()
      // compares two sequence intervals by comparing max value of the dimension, and type
      static bool less_dimension_max(sequence_interval const &__interval1, sequence_interval const &__interval2) {
        return math::interval<size_t>::less_max(__interval1.get_dimension(), __interval2.get_dimension()) ||
               (!(math::interval<size_t>::less_max(__interval2.get_dimension(), __interval1.get_dimension())) &&
                __interval1.get_type() < __interval2.get_type());
      } // less_dimension_max()

      // define ordering
      bool operator<(sequence_interval const &__rhs) const {
        return less_dimension_min_max(*this, __rhs);
      } // operator<()
    }; // class sequence_interval

    // output operator for sequence_interval
    template <typename T>
    inline std::ostream &operator<<(std::ostream &__out, sequence_interval<T> const &__sequence_interval) {
      __out << __sequence_interval.get_dimension() << "->" << *__sequence_interval.get_type();
      return __out;
    } // operator<<()

    using cchb_interval = sequence_interval<cchb>;
  } // namespace che
} // namespace biosim

#endif // che_sequence_interval_h
