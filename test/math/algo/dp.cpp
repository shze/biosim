#include <boost/test/unit_test.hpp>

#include "math/algo/dp.h" // header to test

using namespace biosim;

size_t fibonacci_function(math::tensor<size_t> const &__input, std::vector<size_t> const &__pos) {
  if(__pos.size() != 1) {
    throw std::invalid_argument("fibonacci numbers cannot be calculated for multidimensional positions");
  }
  return __pos[0] < 2 ? 1 : __input({(__pos[0] - 1)}) + __input({(__pos[0] - 2)});
}

BOOST_AUTO_TEST_SUITE(suite_dp)

BOOST_AUTO_TEST_CASE(dp_fib) {
  math::algo::dp<size_t> fib;
  math::tensor<size_t> fibonacci_numbers(fib.calculate(math::tensor<size_t>({10}), &fibonacci_function));
  BOOST_CHECK(fibonacci_numbers({9}) == 55);
  BOOST_REQUIRE_THROW(fib.calculate(math::tensor<size_t>({10, 3}), &fibonacci_function), std::invalid_argument);
}

BOOST_AUTO_TEST_SUITE_END()
