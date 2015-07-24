#ifndef che_score_ev_alignment
#define che_score_ev_alignment

#include "math/function.h"
#include "che/alignment.h"
#include "che/score/cm_cc_blosum.h"

namespace biosim {
  namespace che {
    namespace score {
      // class to evaluate the probability of alignments
      class ev_alignment : public math::ev_function<che::alignment> {
      public:
        // ctor taking a cc score offset; default ctor
        explicit ev_alignment(double __cc_score_offset = 0.0);
        // returns identifier
        std::string get_identifier() const;
        // evaluate the probability of the given alignment
        double evaluate(che::alignment const &__al) const;

      private:
        cm_cc_blosum _cc_score_f; // cc score function
        double _cc_score_offset; // cc score offset
        double _gap_score; // gap score
      }; // class ev_alignment
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_ev_alignment
