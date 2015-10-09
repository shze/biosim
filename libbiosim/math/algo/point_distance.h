#ifndef math_algo_distance_h
#define math_algo_distance_h

#include "math/point.h"
#include <algorithm>

namespace biosim {
  namespace math {
    namespace algo {
      // calculates the power of the minkowski distance between the two given points
      template <class T, size_t N>
      double minkowski_distance_power(std::array<T, N> const &__p1, std::array<T, N> const &__p2,
                                      double const &__power) {
        double d = 0;
        for(size_t pos(0); pos < N; ++pos) {
          d += std::pow(std::abs(__p1[pos] - __p2[pos]), __power);
        } // for
        return d;
      } // minkowski_distance_power()
      // calculates the minkowski distance between the two given points
      template <class T, size_t N>
      double minkowski_distance(std::array<T, N> const &__p1, std::array<T, N> const &__p2, double const &__power) {
        return std::pow(minkowski_distance_power(__p1, __p2, __power), 1.0 / __power);
      } // minkowski_distance()
      // calculates the square of the euclidean distance between the two given points
      template <class T, size_t N>
      double euclidean_distance_power(std::array<T, N> const &__p1, std::array<T, N> const &__p2) {
        return minkowski_distance_power(__p1, __p2, 2.0);
      } // euclidean_distance_power()
      // calculates the euclidean distance between the two given points
      template <class T, size_t N>
      double euclidean_distance(std::array<T, N> const &__p1, std::array<T, N> const &__p2) {
        return minkowski_distance(__p1, __p2, 2.0);
      } // euclidean_distance()
      // calculates the manhattan distance between the two given points
      template <class T, size_t N>
      double manhattan_distance(std::array<T, N> const &__p1, std::array<T, N> const &__p2) {
        return minkowski_distance(__p1, __p2, 1.0);
      } // manhattan_distance()
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_algo_distance_h
