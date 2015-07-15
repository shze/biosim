#include "che/algo/aligner_dp.h"
#include "math/algo/dp.h"

namespace biosim {
  namespace che {
    namespace algo {
      // returns whether an underflow occurs when substracting vector __pos - __minus
      bool tensor_pos_underflow(std::vector<size_t> const &__pos, std::vector<size_t> const &__minus) {
        for(size_t dim(0); dim < __pos.size() && dim < __minus.size(); ++dim) {
          if(__pos[dim] < __minus[dim]) {
            return true;
          } // if
        } // for
        return false;
      } // tensor_pos_underflow()
      // returns the vector substraction __pos - __minus
      std::vector<size_t> tensor_pos_subtract(std::vector<size_t> const &__pos, std::vector<size_t> const &__minus) {
        std::vector<size_t> result(__pos.size(), 0);
        for(size_t dim(0); dim < __pos.size() && dim < __minus.size(); ++dim) {
          result[dim] = __pos[dim] - __minus[dim];
        } // for
        return result;
      } // tensor_pos_subtract()

      // internal class passed to the dp algorithm and used to fill the tensor based on alignments and alignment score
      class tensor_score_function {
      public:
        // ctor from alignments to be aligned and an alignment score function
        tensor_score_function(std::vector<alignment> __alignments, score::ev_alignment __score_f)
            : _alignments(__alignments), _score_f(__score_f) {}
        // operator called by the dp algorithm to calculate the score for the given position in the given tensor
        double operator()(math::tensor<double> const &__input, std::vector<size_t> const &__pos) {
          if(_alignments.size() != __pos.size()) {
            throw std::invalid_argument("need to have as many input alignments as dimensions in tensor position");
          } // if

          // create incrementor with __pos.size length and digits 0=gap, 1=match
          tools::incrementor<std::vector<size_t>> inc(std::vector<std::vector<size_t>>(__pos.size(), {0, 1}));
          std::vector<size_t> direction(__pos.size(), 0); // create starting position, all zero
          direction = inc.next(direction); // do one increment to avoid the all-gap case

          bool done(false);
          double score(-std::numeric_limits<double>::max());
          while(!done) {
            if(!tensor_pos_underflow(__pos, direction)) { // calculate only if the needed position is within the tensor
              std::vector<size_t> previous_pos(tensor_pos_subtract(__pos, direction));
              double previous_score(__input(previous_pos));

              std::vector<che::molecule> molecule_parts; // create a small alignment only with the necessary parts
              for(size_t dim(0); dim < direction.size(); ++dim) {
                che::ps seq;
                seq.push_back(direction[dim] == 1
                                  ? _alignments[dim].get_cc(__pos[dim] - 1, 0)
                                  : che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
                che::molecule m("", "", seq);
                molecule_parts.emplace_back(m);
              } // for
              double this_score(_score_f.evaluate(che::alignment(molecule_parts)));

              score = std::max(score, previous_score + this_score);
            } // if

            try {
              direction = inc.next(direction);
            } catch(std::overflow_error &e) {
              done = true;
            } // catch
          } // while

          return score == -std::numeric_limits<double>::max() ? 0.0 : score;
        } // operator()

      private:
        std::vector<alignment> _alignments; // alignments to be aligned
        score::ev_alignment _score_f; // alignment score function
      }; // class dp_score

      // ctor from alignment score; also default ctor
      aligner_dp::aligner_dp(score::ev_alignment __score_f) : _score_f(__score_f) {}
      // aligns two alignments
      std::list<scored_alignment> aligner_dp::align_pair(alignment const &__al1, alignment const &__al2) const {
        return align_multiple({__al1, __al2});
      } // align_pair()
      // aligns multiple alignments
      std::list<scored_alignment> aligner_dp::align_multiple(std::vector<alignment> const &__alignments) const {
        math::algo::dp<double> dp_algorithm;
        tensor_score_function tensor_f(__alignments, _score_f);
        std::vector<size_t> tensor_size;
        for(auto a : __alignments) {
          tensor_size.push_back(a.get_length() + 1);
        } // for

        // the two key steps: fill the tensor using dp; return the best alignments found in backtracking
        math::tensor<double> filled_tensor(dp_algorithm.calculate(math::tensor<double>(tensor_size), tensor_f));
        return backtrack(filled_tensor, __alignments);
      } // align_multiple()

      // find the best scoring alignments by backtracking through the filled score tensor
      std::list<scored_alignment> aligner_dp::backtrack(math::tensor<double> const &__scores,
                                                        std::vector<alignment> const &__alignments) const {
        return std::list<scored_alignment>();
      } // backtrack()
    } // namespace algo
  } // namespace che
} // namespace biosim
