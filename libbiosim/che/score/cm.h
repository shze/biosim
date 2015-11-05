#ifndef che_score_cm_h
#define che_score_cm_h

#include "math/function.h"
#include "tools/enumerate.h"
#include "che/assembly.h"

namespace biosim {
  namespace che {
    namespace score {
      // enumerates all compare functions for assembly; added manually
      class assembly_compares {
      public:
        using assembly_compare_ptr = std::shared_ptr<math::cm_function<assembly>>; // simplify naming
        using assembly_enum = tools::enumerate<assembly_compare_ptr>; // simplify naming

        // adds new compare function ptrs; this allows cm_function ptrs to be added from outside this class
        static void add(assembly_compare_ptr __obj);
        // returns const ref to static instance_container
        static assembly_enum::instance_container get_instances();

      private:
        // constructs all compare function ptrs
        static bool initialize();
        // static variable to initialize assembly_enum
        static bool _initialized;
      }; // class assembly_compares
    } // namespace score
  } // namespace che
} // namespace biosim

#endif // che_score_cm_h
