#ifndef che_algo_aligner_dp_h
#define che_algo_aligner_dp_h

#include "che/algo/aligner.h"
#include "che/score/ev_alignment.h"
#include "math/tensor.h"

namespace biosim {
  namespace che {
    namespace algo {
      // dynamic programming aligner for pairwise and multiple sequence alignment
      class aligner_dp : public psa_aligner, msa_aligner {
      public:
        // ctor from alignment score function; also default ctor
        explicit aligner_dp(score::ev_alignment __score_f = score::ev_alignment());
        // aligns two alignments
        std::list<scored_alignment> align_pair(alignment const &__al1, alignment const &__al2) const;
        // aligns multiple alignments
        std::list<scored_alignment> align_multiple(std::vector<alignment> const &__alignments) const;

      private:
        // find the best scoring alignments by backtracking through the filled score tensor
        std::list<scored_alignment> backtrack(math::tensor<double> const &__scores,
                                              std::vector<alignment> const &__alignments) const;

        score::ev_alignment _score_f; // alignment score
      }; // class aligner_dp
    } // namespace algo
  } // namespace che
} // namespace biosim

#endif // che_algo_aligner_dp_h
