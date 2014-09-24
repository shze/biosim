#include "che/ss.h"

namespace biosim {
  namespace che {
    // ctor from secondary structure sequence
    ss::ss(std::vector<cchb> __sequence) : _sequence(__sequence), _sses() {}
  } // namespace che
} // namespace biosim