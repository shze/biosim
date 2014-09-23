#ifndef math_interval_h
#define math_interval_h

#include <limits>
#include <list>
#include "tools/log.h"

namespace biosim {
  namespace math {
    // holds an interval (range); see http://stackoverflow.com/questions/1089192/c-stl-range-container
    // implements an closed interval [a, b] = {x | a <= x <= b}, see http://en.wikipedia.org/wiki/Interval_(mathematics)
    template <typename T>
    class interval {
    private:
      T _min, _max; // min, max of the interval

    public:
      // default ctor
      interval()
          : _min(std::numeric_limits<T>::is_signed ? -std::numeric_limits<T>::max() : std::numeric_limits<T>::min()),
            _max(std::numeric_limits<T>::max()) {}
      // ctor from min and max value; use this to construct an interval from a single value
      interval(T __min, T __max) : _min(__min), _max(__max) {
        if(__min > __max) { // if min > max, swap them
          std::swap(_min, _max);
        }
      } // ctor

      // get epsilon for this T (static b/c it is only type dependent, and interval independent)
      static T get_epsilon() {
        // for size_t epsilon was 0, but the next number is 1, otherwise to the numeric_limits epsilon value
        return std::numeric_limits<T>::is_integer ? 1 : std::numeric_limits<T>::epsilon();
      } // get_epsilon()

      // get min of interval
      T const &get_min() const { return _min; }
      // get max of interval
      T const &get_max() const { return _max; }
      // set min of interval
      void set_min(T __min) {
        _min = __min;
        if(__min > _max) {
          _max = __min;
        } // if
      } // set_min()
      // set max of interval
      void set_max(T __max) {
        _max = __max;
        if(__max < _min) {
          _min = __max;
        } // if
      } // set_max()

      // get the length of the interval
      T get_length() const { return _max - _min; }

      // test if this interval and interval2 are continuous (1->2; overlapping or adjacent with a gap <= epsilon)
      bool is_continuous(interval<T> const &__interval2) const {
        // continuous if _max is in [__interval2._min-eps, __interval2._max), and _interval2._min is in (_min, _max+eps]
        return (_max >= __interval2._min - get_epsilon() && _max < __interval2._max) &&
               (__interval2._min > _min && __interval2._min <= _max + get_epsilon());
      } // is_continuous()
      // test if this interval and interval2 overlap (if this.max is in interval2 or interval2.max is in this interval)
      bool overlaps(interval<T> const &__interval2) const {
        return (_max >= __interval2._min && _max <= __interval2._max) ||
               (__interval2._max >= _min && __interval2._max <= _max);
      } // overlaps()

      // returns if interval1 is less than i.e. before (starts further left, at a smaller number) interval2
      static bool less_min_max(interval<T> const &__interval1, interval<T> const &__interval2) {
        return __interval1._min < __interval2._min ||
               (__interval1._min == __interval2._min && __interval1._max < __interval2._max);
      } // less_min_max()
      // returns if the length of interval1 is less than the length of interval2
      static bool less_length(interval<T> const &__interval1, interval<T> const &__interval2) {
        return __interval1.get_length() < __interval2.get_length();
      } // less_length()
      // returns if interval1 ends before before interval2 (ends further left, at a smaller number)
      static bool less_max(interval<T> const &__interval1, interval<T> const &__interval2) {
        return __interval1._max < __interval2._max;
      } // less_max()
      // returns true if this interval is less than the rhs interval
      bool operator<(interval<T> const &__rhs) const {
        bool const result(less_min_max(*this, __rhs));
        DEBUG << "operator<(): " << *this << " < " << __rhs << ", result=" << result;
        return result;
      } // operator<()

      // returns if both intervals have same min and max positions
      static bool equal_min_max(interval<T> const &__interval1, interval<T> const &__interval2) {
        return __interval1._min == __interval2._min && __interval1._max == __interval2._max;
      } // equal_min_max()
      // returns if both intervals have same length
      static bool equal_length(interval<T> const &__interval1, interval<T> const &__interval2) {
        return __interval1.get_length() == __interval2.get_length();
      } // equal_length()
      // returns true if this interval is equal to the rhs interval
      bool operator==(interval<T> const &__rhs) const {
        bool const result(equal_min_max(*this, __rhs));
        DEBUG << "operator==(): " << *this << " == " << __rhs << ", result=" << result;
        return result;
      } // operator==()

      // tries to merge two intervals; returns merged interval if merging is possible, or an empty list if not
      static std::list<interval<T>> merge(interval<T> const &__interval1, interval<T> const &__interval2) {
        std::list<interval<T>> merged_intervals; // to collect result intervals

        if(__interval1.is_continuous(__interval2) || __interval2.is_continuous(__interval1) ||
           __interval1.overlaps(__interval2)) {
          // create a merged interval
          interval merged_interval(std::min(__interval1._min, __interval2._min),
                                   std::max(__interval1._max, __interval2._max));
          merged_intervals.push_back(merged_interval);
        } // if
        // if intervals are not continuous, return an empty list

        return merged_intervals;
      } // merge()
      // tries to merge a list of intervals
      static std::list<interval<T>> merge(std::list<interval<T>> __intervals) {
        std::list<interval<T>> merged_intervals; // create empty merged_intervals
        if(__intervals.empty()) {
          return merged_intervals;
        } // if

        __intervals.sort(); // sort input

        // get first interval, remove from input, add to merged_intervals
        DEBUG << "merge(list): adding " << __intervals.front() << " (first)";
        merged_intervals.push_back(__intervals.front()); // __intervals contains at least of element
        __intervals.pop_front();

        // while not input empty // this loop moves itr2 over the input
        while(!__intervals.empty()) {
          // try to merge last interval of merged_intervals with first in __intervals (they are sorted)
          std::list<interval<T>> merged_sublist(merge(merged_intervals.back(), __intervals.front()));

          if(merged_sublist.size() == 1) { // if merge was successful, replace the last merged interval with this
            DEBUG << "merge(list): replacing " << merged_intervals.back() << " with " << merged_sublist.front();
            merged_intervals.pop_back(); // remove last
            merged_intervals.push_back(merged_sublist.front()); // add new merged interval
          } else { // if( merged_sublist.size() == 0) // if not successful, add new interval to merged_intervals
            DEBUG << "merge(list): adding " << __intervals.front();
            merged_intervals.push_back(__intervals.front()); // add first unprocessed to merged intervals
          }
          __intervals.pop_front(); // remove first from unprocessed intervals (common between both cases)
        }

        return merged_intervals;
      } // merge()
    }; // class interval

    // output operator for interval; inline to avoid multiple definitions error when .h is included multiple times
    template <typename T>
    inline std::ostream &operator<<(std::ostream &__out, interval<T> const &__interval) {
      __out << "[" << __interval.get_min() << ", " << __interval.get_max() << "; epsilon=" << interval<T>::get_epsilon()
            << "]";
      return __out;
    } // operator<<()

    // less for intervals; to sort intervals in containers differently than the default defined in operator<()
    // take a function pointer, so that all less functions implemented in interval<T>::less_* can be used
    template <typename T>
    struct interval_less : std::binary_function<interval<T>, interval<T>, bool> {
      using less_function = bool (*)(interval<T> const &, interval<T> const &); // make reading easier

      less_function _less_function; // member to store the function to call for sorting

      // ctor; no default value as it would be confusing; std::less works for interval as it uses the operator<()
      interval_less(less_function __less_function) : _less_function(__less_function) {}

      // function operator
      bool operator()(interval<T> const &__interval1, interval<T> const &__interval2) {
        return _less_function(__interval1, __interval2);
      } // operator()
    }; // struct interval_less
  } // namespace math
} // namespace biosim

#endif // math_interval_h
