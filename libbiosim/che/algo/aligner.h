#ifndef che_algo_aligner_h
#define che_algo_aligner_h

#include "che/alignment.h"

namespace biosim {
  namespace che {
    namespace algo {
      // interface for aligner doing pairwise sequence/structure alignment (psa)
      class psa_aligner {
      public:
        // aligns two alignments
        virtual std::list<scored_alignment> align_pair(alignment const &__al1, alignment const &__al2) const = 0;
        // aligns an alignment to a number of alignments
        virtual std::list<scored_alignment> align_pairwise(alignment const &__al1,
                                                           std::vector<alignment> const &__alignments) const {
          std::list<scored_alignment> alignments;
          for(auto al2 : __alignments) {
            alignments.splice(alignments.end(), align_pair(__al1, al2));
          } // for

          return alignments;
        } // align_pairwise()
      }; // class psa_aligner

      // interface for aligner doing multiple sequence/structure alignment (msa)
      class msa_aligner {
      public:
        // aligns multiple alignments
        virtual std::list<scored_alignment> align_multiple(std::vector<alignment> const &__alignments) const = 0;
      }; // class msa_aligner
    } // namespace algo
  } // namespace che
} // namespace biosim

#endif // che_algo_aligner_h
