#ifndef che_score_ev_alignment
#define che_score_ev_alignment

#include "che/score/cm_cc_function.h"
#include "che/alignment.h"
#include "che/score/cm_cc_blosum.h"

namespace biosim {
  namespace che {
    namespace score {
      // class to evaluate the probability of alignments
      class ev_alignment : public math::ev_function<che::alignment> {
      public:
        using cc_cm_function_ptr = std::shared_ptr<cc_cm_function>; // simplify naming

        // default ctor
        ev_alignment();
        // ctor taking a cc compare function, cc compare function offset, and gap score
        ev_alignment(cc_cm_function_ptr __cm_f, double __cm_f_offset, double __gap_score);
        // returns identifier
        std::string get_identifier() const;
        // evaluate the probability of the given alignment
        double evaluate(che::alignment const &__al) const;

      private:
        cc_cm_function_ptr _cm_f; // cc compare function
        double _cm_f_offset; // cc compare function offset
        double _gap_score; // gap score
      }; // class ev_alignment
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_ev_alignment
