#ifndef math_algo_dp_h
#define math_algo_dp_h

#include "math/tensor.h"

namespace biosim {
  namespace math {
    namespace algo {
      // dynamic programming algorithm class; fills the given tensor using the scoring function
      template <class T>
      class dp {
      public:
        // scoring function, taking a tensor and position and returning a score
        using score_function = std::function<T(math::tensor<T> const &, std::vector<size_t> const &)>;

        // calculates the scores for all positions of the tensor __input using the scoring function and returns it
        math::tensor<T> calculate(math::tensor<T> __input, score_function const &__score) {
          std::vector<size_t> pos(__input.get_rank(), 0);
          std::vector<std::vector<size_t>> alphabets;
          for(size_t dim(0); dim < __input.get_rank(); ++dim) {
            std::vector<size_t> alphabet(__input.get_size(dim));
            size_t element(0);
            std::generate(alphabet.begin(), alphabet.end(), [&] { return element++; });
            alphabets.push_back(alphabet);
          } // for
          tools::incrementor<std::vector<size_t>> inc(alphabets);
          while(!inc.overflow()) {
            DEBUG << "Calculating score for tensor(" << math::tensor<T>::to_string(pos) << ")";
            __input(pos) = __score(__input, pos);
            pos = inc.next(pos);
          } // while

          return __input;
        } // calculate()
      }; // class dp
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_algo_dp_h
