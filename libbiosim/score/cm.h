#ifndef score_cm_h
#define score_cm_h

#include <memory>
#include "score/function.h"
#include "tools/enumerate.h"
#include "che/assembly.h"

namespace biosim {
  namespace score {
    // define identifier function for enumerate as template for all cm_function ptrs
    template <class T>
    std::string identifier(std::shared_ptr<cm_function<T>> const &__cm_function_ptr) {
      return __cm_function_ptr->get_identifier();
    } // identifier()

    // enumerates all compare functions for assembly; added manually
    class assembly_compares {
    private:
      using assembly_compare_ptr = std::shared_ptr<score::cm_function<che::assembly>>; // simplify naming
      using assembly_enum = tools::enumerate<assembly_compare_ptr>; // simplify naming

      // constructs all compare function ptrs
      static bool initialize();
      // static variable to initialize assembly_enum
      static bool _initialized;

    public:
      // adds new compare function ptrs; this allows cm_function ptrs to be added from outside this class
      static void add(assembly_compare_ptr __obj);
      // returns const ref to static instance_container
      static assembly_enum::instance_container get_instances();
    }; // class assembly_compares
  } // namespace score
} // namespace biosim

#endif // score_cm_h
