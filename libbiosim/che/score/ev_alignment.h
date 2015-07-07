#ifndef che_score_ev_alignment
#define che_score_ev_alignment

#include "math/function.h"
#include "che/alignment.h"

namespace biosim {
  namespace che {
    namespace score {
      // class to evaluate the probability of alignments
      class ev_alignment : public math::ev_function<che::alignment> {
      public:
        // returns identifier
        std::string get_identifier() const;
        // evaluate the probability of the given alignment
        double evaluate(che::alignment const &__al) const;
      }; // class ev_alignment
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_ev_alignment
