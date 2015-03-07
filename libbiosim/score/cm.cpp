#include "score/cm.h"
#include "score/cm_assembly_ssq3.h"

namespace biosim {
  namespace score {
    // (static) constructs all compare function ptrs
    bool assembly_compares::initialize() {
      add(assembly_compare_ptr(new score::cm_assembly_ssq3));

      return true;
    } // initialize()

    bool assembly_compares::_initialized(assembly_compares::initialize()); // initialize static variable

    // (static) adds new compare function ptrs; this allows cm_function ptrs to be added from outside this class
    void assembly_compares::add(assembly_compare_ptr __obj) { assembly_enum::add(__obj); }
    // (static) returns const ref to static instance_container
    assembly_compares::assembly_enum::instance_container assembly_compares::get_instances() {
      return assembly_enum::get_instances();
    } // get_instances()
  } // namespace score
} // namespace biosim