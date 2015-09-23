#ifndef math_average_h
#define math_average_h

#include <cstddef>

namespace biosim {
  namespace math {
    namespace algo {
      // cumulative moving average, see https://en.wikipedia.org/wiki/Moving_average
      class cma {
      public:
        // default ctor
        cma();
        // observe a given value
        void observe(double __value);
        // get current number of observed values
        size_t get_n() const;
        // get current average
        double get_average() const;

      private:
        double _average; // current average
        size_t _n; // current number of observed values
      }; // class cma
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_average_h
