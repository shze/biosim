#ifndef che_score_cm_cc_function_h
#define che_score_cm_cc_function_h

#include "math/function.h"
#include "che/cc.h"

namespace biosim {
  namespace che {
    namespace score {
      // base class for cc comparison functions
      class cc_cm_function : public math::cm_function<che::cc> {
      public:
        // returns the minimum score
        virtual double get_min_score() const = 0;
        // returns the score of comparing two unknown cc
        double get_unknown_score() const { return compare(che::cc('X'), che::cc('X')); }
      }; // class cc_cm_function
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_cm_cc_function_h
