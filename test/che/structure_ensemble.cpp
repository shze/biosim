#include <boost/test/unit_test.hpp>

#include "che/structure_ensemble.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_structure_ensemble)

BOOST_AUTO_TEST_CASE(structure_sample) {
  std::string storage = "some_storage";
  std::string id = "some_id";
  double occ = 0.7;
  che::structure_sample s(che::structure(storage, id), occ);
  BOOST_CHECK(s.sample.get_storage() == storage);
  BOOST_CHECK(s.sample.get_identifier() == id);
  BOOST_CHECK(s.occupancy == occ);
}

BOOST_AUTO_TEST_CASE(structure_ensemble) {
  che::structure_ensemble e{che::structure()};
  BOOST_CHECK(e.get_samples().size() == 1);
  BOOST_CHECK(e.get_samples().front().occupancy == 1.0);

  std::vector<che::structure> v(2);
  e = che::structure_ensemble(v.begin(), v.end());
  BOOST_CHECK(e.get_samples().size() == 2);
  BOOST_CHECK(e.get_samples().front().occupancy == 0.5);
}

BOOST_AUTO_TEST_SUITE_END()
