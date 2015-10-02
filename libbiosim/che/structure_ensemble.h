#ifndef che_structure_ensemble_h
#define che_structure_ensemble_h

#include "che/structure.h"
#include <iterator>

namespace biosim {
  namespace che {
    // structure_sample is one element of a structure_ensemble; essentially a pair, but with specifically named members
    struct structure_sample {
      // ctor from structure and occupancy
      structure_sample(structure __sample, double __occupancy);

      structure sample; // structure
      double occupancy; // occupancy of this structure, 0..1
    }; // struct structure_sample

    // ensemble of structures; no default ctor to avoid an empty ensemble; class to ensure normalized occupancies
    class structure_ensemble {
    public:
      // ctor from single structure
      explicit structure_ensemble(structure const &__sample);
      // ctor from begin and end iterators of container<structure>
      template <class I>
      structure_ensemble(I __begin, I __end) {
        if(std::distance(__begin, __end) == 0) {
          throw std::invalid_argument("cannot create an empty structure_ensemble");
        } // if

        std::for_each(__begin, __end, [&](che::structure __sample) { _samples.emplace_back(__sample, 1.0); });
        _samples = normalize_occupancies(_samples);
      } // ctor
      // get samples
      std::list<structure_sample> const &get_samples() const;

    private:
      // normalize the occupancies
      static std::list<structure_sample> normalize_occupancies(std::list<structure_sample> __samples);

      std::list<structure_sample> _samples; // samples of the ensemble
    }; // class structure_ensemble
  } // namespace che
} // namespace biosim

#endif // che_structure_ensemble_h
