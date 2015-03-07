#include <boost/test/unit_test.hpp>
#include <boost/concept_check.hpp>

#include "score/cm.h" // header to test

using namespace biosim;

class cm_assembly_fun : public score::cm_function<che::assembly> {
public:
  std::string get_identifier() const { return "fun"; }
  double compare(che::assembly const &, che::assembly const &) const { return -1.0; }
};

BOOST_AUTO_TEST_SUITE(suite_cm)

BOOST_AUTO_TEST_CASE(cm_instances) { BOOST_CHECK(score::assembly_compares::get_instances().size() == 1); }

BOOST_AUTO_TEST_CASE(cm_add) {
  size_t instance_count(score::assembly_compares::get_instances().size());
  score::assembly_compares::add(std::shared_ptr<score::cm_function<che::assembly>>(new cm_assembly_fun));
  BOOST_CHECK(score::assembly_compares::get_instances().size() == instance_count + 1);
}

BOOST_AUTO_TEST_SUITE_END()
