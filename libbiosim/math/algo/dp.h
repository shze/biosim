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
        using score_function = std::function<T(tensor<T> const &, std::vector<size_t> const &)>;

        // calculates the scores for all positions of the tensor __input using the scoring function and returns it
        tensor<T> calculate(tensor<T> __input, score_function const &__score) {
          tools::mapper<std::vector<size_t>> m(__input.get_mapper_alphabets());
          std::vector<size_t> pos;
          for(size_t i(m.get_min()); i <= m.get_max(); ++i) {
            pos = m.encode(i);
            DEBUG << "Calculating score for tensor(" << tools::to_string(pos) << ")";
            __input(pos) = __score(__input, pos);
          } // while

          return __input;
        } // calculate()
      }; // class dp
    } // namespace algo
  } // namespace math
} // namespace biosim

#endif // math_algo_dp_h
