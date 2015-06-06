#include "che/score/cm_assembly_ssq3.h"

namespace biosim {
  namespace che {
    namespace score {
      // returns identifier
      std::string cm_assembly_ssq3::get_identifier() const { return "ssq3"; }
      // compares the given two instances of assemblies
      double cm_assembly_ssq3::compare(che::assembly const &__first, che::assembly const &__second) const {
        std::list<std::string> chain_ids1(__first.get_chain_id_list()), chain_ids2(__second.get_chain_id_list());
        std::set<std::string> unique_chain_ids(chain_ids1.begin(), chain_ids1.end());
        unique_chain_ids.insert(chain_ids2.begin(), chain_ids2.end());

        size_t chain_count(0);
        double q3_sum(0.0);
        for(auto chain_id : unique_chain_ids) {
          if(!__first.has_ss(chain_id) || !__second.has_ss(chain_id)) {
            DEBUG << get_identifier() << ": chain id '" << chain_id << "' does not exist in both structures, ignoring.";
            continue;
          }

          q3_sum += compare(__first.get_ss(chain_id), __second.get_ss(chain_id));
          ++chain_count;
        }

        return q3_sum / chain_count;
      } // compare()
      // compares the given two instances of ss
      double cm_assembly_ssq3::compare(che::ss const &__first, che::ss const &__second) const {
        che::sequence<che::cchb_dssp> const seq1(__first.get_sequence()), seq2(__second.get_sequence());

        if(seq1.size() != seq2.size()) {
          throw math::compare_error(get_identifier() + ": Sequence length differ, sequence1.length=" +
                                    std::to_string(seq1.size()) + ", sequence2.length=" + std::to_string(seq2.size()));
        } // if

        size_t c_correct(0), h_correct(0), e_correct(0); // initialize
        for(size_t pos(0); pos < seq1.size(); ++pos) { // works for both sequences b/c of same length
          char const sequence1_ss(seq1[pos].get_identifier_char());
          char const sequence2_ss(seq2[pos].get_identifier_char());
          if(sequence1_ss == sequence2_ss) { // actually it's only important that they are the same, not which one
            if(sequence1_ss == 'C') {
              ++c_correct;
            } // if
            else if(sequence1_ss == 'H') {
              ++h_correct;
            } // else if
            else if(sequence1_ss == 'E') {
              ++e_correct;
            } // else if
          } // if
        } // for

        DEBUG << get_identifier() << ": c_correct=" << c_correct << " h_correct=" << h_correct
              << " e_correct=" << e_correct << " seq_len=" << seq1.size();

        return ((double)c_correct + h_correct + e_correct) / seq1.size();
      } // compare()
    } // namespace score
  } // namespace che
} // namespace biosim
