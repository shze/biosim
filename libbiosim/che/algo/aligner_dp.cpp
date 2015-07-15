#include "che/algo/aligner_dp.h"
#include "math/algo/dp.h"
#include <boost/algorithm/string.hpp>

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
          while(!done) { // go over all directions
            if(!tensor_pos_underflow(__pos, direction)) { // calculate only if the needed position is within the tensor
              std::vector<size_t> previous_pos(tensor_pos_subtract(__pos, direction));
              double previous_score(__input(previous_pos));

              std::vector<che::molecule> molecule_parts; // create a small alignment only with the necessary parts
              for(size_t dim(0); dim < direction.size(); ++dim) {
                che::ps seq;
                seq.push_back(direction[dim] == 1
                                  ? _alignments[dim].get_cc(__pos[dim] - 1, 0)
                                  : che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
                molecule_parts.emplace_back("", "", seq);
              } // for
              double step_score(_score_f.evaluate(che::alignment(molecule_parts)));

              score = std::max(score, previous_score + step_score);
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

      // internal data structure for building the alignment
      struct alignment_data {
        // ctor from size
        explicit alignment_data(size_t __size) : _molecule_data() {
          _molecule_data.insert(_molecule_data.begin(), __size, std::list<che::cc>());
        } // ctor

        std::vector<std::list<che::cc>> _molecule_data; // a list<cc> for each molecule
      }; // struct alignment_data

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
        std::list<std::pair<std::vector<size_t>, alignment_data>> work; // list of work items
        std::list<alignment_data> best_alignment_data; // list of completed alignments

        // create first work item
        std::vector<size_t> last_pos; // fill in last position
        for(size_t dim(0); dim < __scores.get_rank(); ++dim) {
          last_pos.push_back(__scores.get_size(dim) - 1);
        } // for
        work.emplace_back(std::make_pair(last_pos, alignment_data(__alignments.size())));

        // create incrementor with __pos.size length and digits 0=gap, 1=match, for use in while
        tools::incrementor<std::vector<size_t>> inc(std::vector<std::vector<size_t>>(__scores.get_rank(), {0, 1}));
        // backtracking
        while(!work.empty()) {
          std::vector<size_t> current_pos(work.front().first);
          alignment_data current_alignment_data(work.front().second);
          work.pop_front();

          std::vector<size_t> direction(__scores.get_rank(), 0); // create starting position, all zero
          direction = inc.next(direction); // do one increment to avoid the all-gap case

          bool done(false);
          while(!done) { // go over all directions
            if(!tensor_pos_underflow(current_pos, direction)) { // calculate only if needed position is in tensor
              std::vector<size_t> previous_pos(tensor_pos_subtract(current_pos, direction));
              double previous_score(__scores(previous_pos));

              std::vector<che::molecule> molecule_parts; // create a small alignment only with the necessary parts
              for(size_t dim(0); dim < direction.size(); ++dim) {
                che::ps seq;
                seq.push_back(direction[dim] == 1
                                  ? __alignments[dim].get_cc(current_pos[dim] - 1, 0)
                                  : che::cc(che::cc::specificity_type::gap, che::cc::monomer_type::l_peptide_linking));
                molecule_parts.emplace_back("", "", seq);
              } // for
              double step_score(_score_f.evaluate(che::alignment(molecule_parts)));

              if(previous_score + step_score == __scores(current_pos)) {
                alignment_data copy(current_alignment_data); // make copy to not affect other directions
                for(size_t dim(0); dim < copy._molecule_data.size(); ++dim) {
                  copy._molecule_data[dim].push_front(molecule_parts[dim].get_ps()[0]); // extend alignment
                } // for

                if(previous_pos == std::vector<size_t>(__scores.get_rank(), 0)) {
                  best_alignment_data.push_back(copy);
                } // if
                else {
                  work.push_back(std::make_pair(previous_pos, copy));
                } // else
              } // if
            } // if

            try {
              direction = inc.next(direction);
            } catch(std::overflow_error &e) {
              done = true;
            } // catch
          } // while
        } // while

        // convert alignment_data into scored_alignment
        std::list<scored_alignment> alignments;
        for(alignment_data a : best_alignment_data) {
          std::vector<molecule> molecules;

          std::vector<std::list<che::cc>>::const_iterator m_itr(a._molecule_data.begin()),
              m_itr_end(a._molecule_data.end());
          std::vector<alignment>::const_iterator a_itr(__alignments.begin()), a_itr_end(__alignments.end());
          for(; m_itr != m_itr_end && a_itr != a_itr_end; ++m_itr, ++a_itr) {
            std::vector<std::string> storages, identifiers; // to store all the parts of storage and identifier
            for(che::molecule m : a_itr->get_molecules()) { // for each molecule, use part of storage and identifier
              std::vector<std::string> parts;
              boost::algorithm::split(parts, m.get_storage(), boost::is_any_of("/"));
              storages.emplace_back(parts.size() > 1 ? parts[parts.size() - 2] + "/" + parts[parts.size() - 1] : "~");
              boost::algorithm::split(parts, m.get_identifier(), boost::is_any_of(" "));
              identifiers.emplace_back(parts.size() > 0 ? parts[0] : "sequence");
            }
            std::string storage(boost::algorithm::join(storages, "+"));
            std::string identifier(boost::algorithm::join(identifiers, "+"));

            ps seq; // create the sequence
            seq.insert(seq.begin(), m_itr->begin(), m_itr->end()); // insert the aligned sequence
            molecules.emplace_back(storage, identifier, seq); // create the final molecule from all parts
          }

          alignments.emplace_back(alignment(molecules), __scores(last_pos));
        } // for

        return alignments;
      } // backtrack()
    } // namespace algo
  } // namespace che
} // namespace biosim
