#ifndef score_cm_assembly_ssq3_h
#define score_cm_assembly_ssq3_h

#include "score/function.h"
#include "che/assembly.h"

namespace biosim {
  namespace score {
    // computes the q3 distance [0..1] between the two secondary structures of two assemblies
    class cm_assembly_ssq3 : public cm_function<che::assembly> {
    public:
      // returns identifier
      std::string get_identifier() const;
      // compares the given two instances of assemblies
      double compare(che::assembly const &__first, che::assembly const &__second) const;
      // compares the given two instances of ss
      double compare(che::ss const &__first, che::ss const &__second) const;
    }; // class cm_assembly_ssq3
  } // namespace score
} // namespace biosim

#endif // score_cm_assembly_ssq3_h
