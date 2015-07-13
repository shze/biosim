#include "che/score/ev_alignment.h"
#include "tools/log.h"
#include "che/score/cm_cc_blosum.h"
#include "math/algo/factorial.h"
#include <limits>

namespace biosim {
  namespace che {
    namespace score {
      // returns identifier
      std::string ev_alignment::get_identifier() const { return "alignment"; }
      // evaluate the probability of the given alignment
      double ev_alignment::evaluate(const alignment &__al) const {
        cm_cc_blosum s(true);
        double min_score(std::numeric_limits<double>::max()); // get minimal score, use as gap score
        for(size_t i(0); i < s.get_dbl_bitscore_matrix().size1(); ++i) {
          for(size_t j(0); j <= i; ++j) {
            min_score = std::min(min_score, s.get_dbl_bitscore_matrix()(i, j));
          } // for
        } // for

        double total_score(0.0);
        double weight(1 / math::algo::binomial_coefficent(__al.get_depth(), 2));
        DEBUG << "Evaluate alignment: length=" << __al.get_length() << "; depth=" << __al.get_depth()
              << "; weight=" << weight;
        for(size_t pos(0); pos < __al.get_length(); ++pos) {
          for(int depth1(0); depth1 < __al.get_depth() - 1; ++depth1) {
            for(int depth2(depth1 + 1); depth2 < __al.get_depth(); ++depth2) {
              cc const &cc1(__al.get_cc(pos, depth1)), &cc2(__al.get_cc(pos, depth2));
              double score(cc1.is_gap() || cc2.is_gap() ? min_score : s.compare(cc1, cc2));
              total_score += weight * score;
              DEBUG << "Evaluate alignment: pos=" << pos << "; depth_pair=(" << depth1 << ", " << depth2
                    << "); score=" << score << "; total_score=" << total_score;
            } // for
          } // for
        } // for

        return total_score;
      } // evaluate()
    } // namespace score
  } // namespace che
} // namespace biosim
