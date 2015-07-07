#include "math/algo/factorial.h"
#include <stdexcept>

namespace biosim {
  namespace math {
    namespace algo {
      // calculates n factorial (n!); even with long double, results are imprecise for n > 25, inf for n > 1754 on
      // x86_64; higher values can only be calculated using bigger datatypes like BigInt, etc and better algorithms, see
      // http://www.luschny.de/math/factorial/FastFactorialFunctions.htm
      long double factorial(unsigned short __n) {
        long double fac(1.0);
        for(; __n > 1; --__n) {
          fac *= __n;
        } // for
        return fac;
      } // factorial()
      // calculates the binomial coefficent; the product algorithm allows complex n (here only double), integer k;
      // generalization allowing complex k is possible, see https://en.wikipedia.org/wiki/Binomial_coefficient
      long double binomial_coefficent(double __n, short __k) {
        if(__k < 0) {
          return 0;
        }
        if(__k == 0) {
          return 1;
        }
        // k > 0
        long double result(__n - __k + 1); // initialize with value for i = 1 b/c k > 0
        for(size_t i(2); i <= __k; ++i) {
          result = result * (__n - __k + i) / i;
        } // for
        return result;
      } // binomial_coefficent
    } // namespace algo
  } // namespace math
} // namespace biosim