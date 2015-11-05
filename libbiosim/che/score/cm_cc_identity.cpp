#include "che/score/cm_cc_identity.h"
#include "tools/log.h"

namespace biosim {
  namespace che {
    namespace score {
      // default ctor
      cm_cc_identity::cm_cc_identity() : _id("identity") {}
      // returns identifier
      std::string cm_cc_identity::get_identifier() const { return std::string(_id); }
      // returns the minimum score
      double cm_cc_identity::get_min_score() const { return std::min(get_match_score(), get_mismatch_score()); }
      // compares the given two instances of cc
      double cm_cc_identity::compare(cc const &__cc1, cc const &__cc2) const {
        cc::weight_map map1(__cc1.get_weights()), map2(__cc2.get_weights()); // use maps to handle profile cc
        double total_score(0.0);
        for(auto p1 : map1) { // iterate through each combination
          for(auto p2 : map2) {
            cc base_cc1(p1.first), base_cc2(p2.first);
            bool is_match(base_cc1.get_identifier_char() == base_cc2.get_identifier_char());
            double score((is_match ? get_match_score() : get_mismatch_score()) * p1.second * p2.second);
            total_score += score;
            DEBUG << "Identity: compare cc1=" << __cc1.get_identifier_char() << " (s=" << (int)__cc1.get_specificity()
                  << ") to cc2=" << __cc2.get_identifier_char() << " (s=" << (int)__cc2.get_specificity()
                  << "), map_cc1=" << p1.first << " (" << base_cc1.get_identifier_char() << "; w=" << p1.second
                  << ") to map_cc2=" << p2.first << " (" << base_cc2.get_identifier_char() << "; w=" << p2.second
                  << "): score=" << score;
          } // for
        } // for

        return total_score;
      } // compare()

      // (static) get score for matching cc
      double cm_cc_identity::get_match_score() { return 1.0; }
      // (static) get score for mismatching cc
      double cm_cc_identity::get_mismatch_score() { return -1.0; }
    } // namespace score
  } // namespace che
} // namespace biosim
