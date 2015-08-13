#ifndef che_score_cm_cc_identity_h
#define che_score_cm_cc_identity_h

#include "che/score/cm_cc_function.h"
#include "che/cc.h"

namespace biosim {
  namespace che {
    namespace score {
      // comparison function for comparing two cc based on identity
      class cm_cc_identity : public cc_cm_function {
      public:
        // default ctor
        cm_cc_identity();
        // returns identifier
        std::string get_identifier() const;
        // returns the minimum score
        double get_min_score() const;
        // compares the given two instances of cc
        double compare(che::cc const &__cc1, che::cc const &__cc2) const;

      private:
        std::string _id; // identifier

        // get score for matching cc
        static double get_match_score();
        // get score for mismatching cc
        static double get_mismatch_score();
      }; // class cm_cc_identity
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_cm_cc_identity_h
