#include "che/score/ev_alignment.h"
#include "tools/log.h"
#include "math/algo/factorial.h"

namespace biosim {
  namespace che {
    namespace score {
      // default ctor
      ev_alignment::ev_alignment()
          : _cm_f(new cm_cc_blosum(true)), _cm_f_offset(0.0), _gap_score(_cm_f->get_min_score()) {}
      // ctor taking a cc compare function, cc score offset, and gap score
      ev_alignment::ev_alignment(ev_alignment::cc_cm_function_ptr __cm_f, double __cm_f_offset, double __gap_score)
          : _cm_f(__cm_f), _cm_f_offset(__cm_f_offset), _gap_score(__gap_score) {}
      // returns identifier
      std::string ev_alignment::get_identifier() const { return "alignment"; }
      // evaluate the probability of the given alignment
      double ev_alignment::evaluate(const alignment &__al) const {
        double total_score(0.0);
        double weight(1 / math::algo::binomial_coefficent(__al.get_depth(), 2));
        DEBUG << "Evaluate alignment: length=" << __al.get_length() << "; depth=" << __al.get_depth()
              << "; weight=" << weight;
        for(size_t pos(0); pos < __al.get_length(); ++pos) {
          for(int depth1(0); depth1 < __al.get_depth() - 1; ++depth1) {
            for(int depth2(depth1 + 1); depth2 < __al.get_depth(); ++depth2) {
              cc const &cc1(__al.get_cc(pos, depth1)), &cc2(__al.get_cc(pos, depth2));
              double score(cc1.is_gap() || cc2.is_gap() ? _gap_score : _cm_f->compare(cc1, cc2));
              total_score += weight * score;
              DEBUG << "Evaluate alignment: pos=" << pos << "; depth_pair=(" << depth1 << ", " << depth2
                    << "); score=" << score << "; total_score=" << total_score;
            } // for
          } // for
        } // for

        return total_score + _cm_f_offset * __al.get_length();
      } // evaluate()
    } // namespace score
  } // namespace che
} // namespace biosim
