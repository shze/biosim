#ifndef math_interval_scheduler_maximize_h
#define math_interval_scheduler_maximize_h

#include "math/algo/interval_scheduler.h"
#include "tools/less.h"
#include <functional>

namespace biosim {
  namespace math {
    namespace algo {
      // interval scheduling algorithm maximizing the number of intervals by selecting the earliest ending intervals
      template <class T>
      class interval_scheduler_maximize : public interval_scheduler<T> {
      public:
        // method to select a non-overlapping subset of intervals from the input set
        std::set<T> schedule(std::set<T> __intervals) {
          // sort pool by ending, use multiset as multiple intervals could end at the same position
          using interval_max_less_multiset = std::multiset<T, tools::less<T, less_max<T>>>;
          interval_max_less_multiset end_sorted_intervals(__intervals.begin(), __intervals.end());

          // create an empty set that will collect all sequence intervals that are converted into the result sequence
          std::set<T> maximized_intervals;

          // while not empty: take first, remove all overlapping
          DEBUG << "end_sorted_intervals.size=" << end_sorted_intervals.size()
                << "; maximized_intervals.size=" << maximized_intervals.size() << " (begin)";
          while(!end_sorted_intervals.empty()) {
            // there might be multipe intervals ending at the same position as the first one; this converter is a simple
            // method for converting non-overlapping intervals in a sequence, and for this scenario this is a non-issue;
            // so if multiple intervals end at the same position, assume the first one is good enough, ignore the rest.
            typename interval_max_less_multiset::const_iterator begin_itr(end_sorted_intervals.begin());
            T new_interval(*begin_itr); // save sequence_interval somewhere outside the set
            end_sorted_intervals.erase(begin_itr); // remove this sequence_interval from the set
            maximized_intervals.insert(new_interval);
            DEBUG << "end_sorted_intervals.size=" << end_sorted_intervals.size()
                  << "; maximized_intervals.size=" << maximized_intervals.size() << " (interval " << new_interval
                  << ": end_sorted_intervals->max_count_intervals)";

            // remove all intervals that are still in end_sorted_intervals and overlap with new_interval; use itrs here
            // and not a range-based for loop to only remove the specific interval, not all interval with this end
            for(typename interval_max_less_multiset::const_iterator other_interval_itr(end_sorted_intervals.begin()),
                itr_end(end_sorted_intervals.end());
                other_interval_itr != itr_end;) {
              if(new_interval.overlaps(*other_interval_itr)) {
                // print before erasing the itr
                DEBUG << "end_sorted_intervals.size=" << end_sorted_intervals.size()
                      << "; maximized_intervals.size=" << maximized_intervals.size() << " (overlapping interval "
                      << *other_interval_itr << ": end_sorted_intervals->removed)";
                // original other_interval_itr is invalid after erase(), assign next element itr that erase() returns
                other_interval_itr = end_sorted_intervals.erase(other_interval_itr);
              } // if
              else {
                ++other_interval_itr;
              } // else
            } // for
          } // while
          DEBUG << "end_sorted_intervals.size=" << end_sorted_intervals.size()
                << "; maximized_intervals.size=" << maximized_intervals.size() << " (finished)";

          return maximized_intervals;
        } // schedule()
      }; // class interval_scheduler_maximize
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_interval_scheduler_maximize_h
