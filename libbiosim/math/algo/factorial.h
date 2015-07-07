#ifndef math_factorial_h
#define math_factorial_h

namespace biosim {
  namespace math {
    namespace algo {
      // calculates n factorial (n!); even with long double, results are imprecise for n > 25, inf for n > 1754 on
      // x86_64; higher values can only be calculated using bigger datatypes like BigInt, etc and better algorithms, see
      // http://www.luschny.de/math/factorial/FastFactorialFunctions.htm
      long double factorial(unsigned short __n);
      // calculates the binomial coefficent
      long double binomial_coefficent(unsigned short __n, unsigned short __k);
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_factorial_h
