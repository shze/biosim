#include "che/structure_ensemble.h"

namespace biosim {
  namespace che {
    // ctor from structure and occupancy
    structure_sample::structure_sample(structure __sample, double __occupancy)
        : sample(__sample), occupancy(__occupancy) {}

    // ctor from single structure
    structure_ensemble::structure_ensemble(structure const &__sample) : _samples(1, structure_sample(__sample, 1.0)) {}
    // get samples
    std::list<structure_sample> const &structure_ensemble::get_samples() const { return _samples; }
    // (static) normalize the occupancies
    std::list<structure_sample> structure_ensemble::normalize_occupancies(std::list<structure_sample> __samples) {
      double sum(0.0);
      for(auto const &p : __samples) {
        sum += p.occupancy;
      } // for
      for(auto &p : __samples) {
        p.occupancy /= sum;
      } // for
      return __samples;
    } // normalize_occupancies()
  } // namespace che
} // namespace biosim