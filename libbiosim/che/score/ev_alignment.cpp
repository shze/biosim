#include "che/score/ev_alignment.h"
#include "tools/log.h"
#include "che/score/cm_cc_blosum.h"
#include "math/algo/factorial.h"

namespace biosim {
  namespace che {
    namespace score {
      // returns identifier
      std::string ev_alignment::get_identifier() const { return "alignment"; }
      // evaluate the probability of the given alignment
      double ev_alignment::evaluate(const alignment &__al) const {
        double total_score(0.0);
        if(__al.get_depth() < 2) { // binomial_coefficent needs alignment depth >= 2
          return total_score;
        } // if

        cm_cc_blosum s(true);
        double weight(1 / math::algo::binomial_coefficent(__al.get_depth(), 2));
        LOG << "Evaluate alignment: length=" << __al.get_length() << "; depth=" << __al.get_depth()
            << "; weight=" << weight;
        for(size_t pos(0); pos < __al.get_length(); ++pos) {
          for(int depth1(0); depth1 < __al.get_depth() - 1; ++depth1) {
            for(int depth2(depth1 + 1); depth2 < __al.get_depth(); ++depth2) {
              double score(weight * s.compare(__al.get_cc(pos, depth1), __al.get_cc(pos, depth2)));
              total_score += score;
              LOG << "Evaluate alignment: pos=" << pos << "; depth_pair=(" << depth1 << ", " << depth2
                  << "); score=" << score << "; total_score=" << total_score;
            } // for
          } // for
        } // for

        return total_score;
      } // evaluate()
    } // namespace score
  } // namespace che
} // namespace biosim
