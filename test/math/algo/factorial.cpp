#include <boost/test/unit_test.hpp>

#include "math/algo/factorial.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_factorial)

BOOST_AUTO_TEST_CASE(factorial) {
  BOOST_CHECK(std::isinf(math::algo::factorial(-1)));
  BOOST_CHECK(math::algo::factorial(0) == 1.0);
  BOOST_CHECK(math::algo::factorial(1) == 1.0);
  BOOST_CHECK(math::algo::factorial(2) == 2.0);
  BOOST_CHECK(math::algo::factorial(25) == 15511210043330985984000000.0L);
  BOOST_CHECK_CLOSE(math::algo::factorial(1754), 1.97926189010501005367e+4930L, 1e-20);
  BOOST_CHECK(std::isinf(math::algo::factorial(1755)));
}

BOOST_AUTO_TEST_CASE(binomial_coefficient) {
  BOOST_CHECK(math::algo::binomial_coefficent(0, 0) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(1, 0) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(1, 1) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(2, 0) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(2, 1) == 2);
  BOOST_CHECK(math::algo::binomial_coefficent(2, 2) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(3, 0) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(3, 1) == 3);
  BOOST_CHECK(math::algo::binomial_coefficent(3, 2) == 3);
  BOOST_CHECK(math::algo::binomial_coefficent(3, 3) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(4, 0) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(4, 1) == 4);
  BOOST_CHECK(math::algo::binomial_coefficent(4, 2) == 6);
  BOOST_CHECK(math::algo::binomial_coefficent(4, 3) == 4);
  BOOST_CHECK(math::algo::binomial_coefficent(4, 4) == 1);

  BOOST_CHECK(math::algo::binomial_coefficent(0, 1) == 0);
  BOOST_CHECK(math::algo::binomial_coefficent(0, 2) == 0);
  BOOST_CHECK(math::algo::binomial_coefficent(49, 6) == 13983816);
  BOOST_CHECK(math::algo::binomial_coefficent(2.5, 2) == 1.875);
  BOOST_CHECK(math::algo::binomial_coefficent(-1, 2) == 1);
  BOOST_CHECK(math::algo::binomial_coefficent(-1, 3) == -1);
}

BOOST_AUTO_TEST_SUITE_END()
