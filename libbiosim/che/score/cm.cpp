#include "che/score/cm.h"
#include "che/score/cm_assembly_ssq3.h"

namespace biosim {
  namespace che {
    namespace score {
      // (static) adds new compare function ptrs; this allows cm_function ptrs to be added from outside this class
      void assembly_compares::add(assembly_compare_ptr __obj) { assembly_enum::add(__obj); }
      // (static) returns const ref to static instance_container
      assembly_compares::assembly_enum::instance_container assembly_compares::get_instances() {
        return assembly_enum::get_instances();
      } // get_instances()

      // (static) constructs all compare function ptrs
      bool assembly_compares::initialize() {
        add(assembly_compare_ptr(new cm_assembly_ssq3));
        return true;
      } // initialize()

      bool assembly_compares::_initialized(assembly_compares::initialize()); // initialize static variable
    } // namespace score
  } // namespace che
} // namespace biosim