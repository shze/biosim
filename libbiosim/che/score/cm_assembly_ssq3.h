#ifndef che_score_cm_assembly_ssq3_h
#define che_score_cm_assembly_ssq3_h

#include "math/function.h"
#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace score {
      // computes the q3 distance [0..1] between the two secondary structures of two assemblies
      class cm_assembly_ssq3 : public math::cm_function<che::assembly> {
      public:
        // returns identifier
        std::string get_identifier() const;
        // compares the given two instances of assemblies
        double compare(che::assembly const &__first, che::assembly const &__second) const;
        // compares the given two instances of ss
        double compare(che::ss const &__first, che::ss const &__second) const;
      }; // class cm_assembly_ssq3
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_cm_assembly_ssq3_h
