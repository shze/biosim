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
      // calculates the binomial coefficent
      long double binomial_coefficent(unsigned short __n, unsigned short __k) {
        if(__n < __k) {
          throw std::invalid_argument("binomial_coefficent cannot be calculated for n < k");
        } // if
        return factorial(__n) / factorial(__k) / factorial(__n - __k);
      } // binomial_coefficent()
    } // namespace algo
  } // namespace math
} // namespace biosim