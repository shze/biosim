#include <boost/test/unit_test.hpp>

#include "math/floating_point.h" // header to test
#include <map>
#include "tools/log.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_floating_point)

BOOST_AUTO_TEST_CASE(floating_point_ctor) {
  math::floating_point<float> f;
  BOOST_CHECK(f.value == 0.0);
  f = math::floating_point<float>(1.2);
  BOOST_CHECK(f.value == float(1.2));

  math::floating_point<double> d;
  BOOST_CHECK(d.value == 0.0);
  d = math::floating_point<double>(1.2);
  BOOST_CHECK(d.value == 1.2);
}

BOOST_AUTO_TEST_CASE(floating_point_cmp) {
  std::map<int, bool> fp_distance_equality{
      {0, true}, {1, true}, {2, true}, {3, true}, {4, true}, {5, false}, {10, false}};
  std::list<float> float_numbers{0, 1, 1e4, 1e30};
  std::list<double> double_numbers{0, 1, 1e4, 1e30};

  BOOST_CHECK(math::floating_point<float>::get_max_ulp_tolerance() == 4);

  BOOST_CHECK(math::almost_equal((float)0.0, (float)0.0) == true);

  for(auto p : fp_distance_equality) {
    for(auto n : float_numbers) {
      BOOST_CHECK(math::almost_equal(n, boost::math::float_advance(n, p.first)) == p.second);
      BOOST_CHECK(math::almost_equal(n, boost::math::float_advance(n, -p.first)) == p.second);
    } // for
  } // for

  math::floating_point<float> f0, f1;
  f1.sign = 1;
  BOOST_CHECK(math::almost_equal(f0, f1) == true);

  f0.value = 0;
  f0.exponent = (1 << math::floating_point<float>::get_exponent_bit_count()) - 1;
  f0.mantissa = 1;
  f1 = f0;
  f1.sign = 1;
  math::floating_point<float> f2(f1);
  f2.mantissa = 10;
  BOOST_CHECK(math::almost_equal(f0, f0) == false);
  BOOST_CHECK(math::almost_equal(f0, f1) == false);
  BOOST_CHECK(math::almost_equal(f0, f2) == false);
  BOOST_CHECK(math::almost_equal(f1, f2) == false);
  BOOST_CHECK(math::almost_equal(f0, math::floating_point<float>(1.0)) == false);

  f0.value = 0;
  f0.exponent = (1 << math::floating_point<float>::get_exponent_bit_count()) - 1;
  f0.mantissa = 0;
  f1 = f0;
  f1.sign = 1;
  f2 = f0;
  --f2.representation;
  math::floating_point<float> f3(f2);
  f3.sign = 1;
  BOOST_CHECK(math::almost_equal(f0, f0) == true);
  BOOST_CHECK(math::almost_equal(f0, f2) == true);
  BOOST_CHECK(math::almost_equal(f1, f1) == true);
  BOOST_CHECK(math::almost_equal(f1, f3) == true);
  BOOST_CHECK(math::almost_equal(f0, f1) == false);
  BOOST_CHECK(math::almost_equal(f0, f3) == false);
  BOOST_CHECK(math::almost_equal(f1, f2) == false);

  BOOST_CHECK(math::floating_point<double>::get_max_ulp_tolerance() == 4);

  BOOST_CHECK(math::almost_equal((double)0.0, (double)0.0) == true);

  for(auto p : fp_distance_equality) {
    for(auto n : double_numbers) {
      BOOST_CHECK(math::almost_equal(n, boost::math::float_advance(n, p.first)) == p.second);
      BOOST_CHECK(math::almost_equal(n, boost::math::float_advance(n, -p.first)) == p.second);
    } // for
  } // for

  math::floating_point<double> d0, d1;
  d1.sign = 1;
  BOOST_CHECK(math::almost_equal(d0, d1) == true);

  d0.value = 0;
  d0.exponent = (1 << math::floating_point<double>::get_exponent_bit_count()) - 1;
  d0.mantissa = 1;
  d1 = d0;
  d1.sign = 1;
  math::floating_point<double> d2(d1);
  d2.mantissa = 10;
  BOOST_CHECK(math::almost_equal(d0, d0) == false);
  BOOST_CHECK(math::almost_equal(d0, d1) == false);
  BOOST_CHECK(math::almost_equal(d0, d2) == false);
  BOOST_CHECK(math::almost_equal(d1, d2) == false);
  BOOST_CHECK(math::almost_equal(d0, math::floating_point<double>(1.0)) == false);

  d0.value = 0;
  d0.exponent = (1 << math::floating_point<double>::get_exponent_bit_count()) - 1;
  d0.mantissa = 0;
  d1 = d0;
  d1.sign = 1;
  d2 = d0;
  --d2.representation;
  math::floating_point<double> d3(d2);
  d3.sign = 1;
  BOOST_CHECK(math::almost_equal(d0, d0) == true);
  BOOST_CHECK(math::almost_equal(d0, d2) == true);
  BOOST_CHECK(math::almost_equal(d1, d1) == true);
  BOOST_CHECK(math::almost_equal(d1, d3) == true);
  BOOST_CHECK(math::almost_equal(d0, d1) == false);
  BOOST_CHECK(math::almost_equal(d0, d3) == false);
  BOOST_CHECK(math::almost_equal(d1, d2) == false);
}

BOOST_AUTO_TEST_SUITE_END()
