#include <boost/test/unit_test.hpp>

#include "math/tensor.h" // header to test
#include "tools/log.h"
#include "tools/incrementor.h"

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_tensor)

BOOST_AUTO_TEST_CASE(tensor_rank0) {
  double d(3.14);
  math::tensor<double> dbl_tensor({});
  BOOST_CHECK(dbl_tensor.get_rank() == 0);
  BOOST_REQUIRE_THROW(dbl_tensor.get_size(0), std::out_of_range);
  dbl_tensor({}) = d;
  BOOST_CHECK(dbl_tensor({}) == d);
  BOOST_REQUIRE_THROW(dbl_tensor({5}), std::out_of_range);

  float f(2.17);
  math::tensor<float> flt_tensor({});
  BOOST_CHECK(flt_tensor.get_rank() == 0);
  BOOST_REQUIRE_THROW(flt_tensor.get_size(0), std::out_of_range);
  flt_tensor({}) = f;
  BOOST_CHECK(flt_tensor({}) == f);
  BOOST_REQUIRE_THROW(flt_tensor({5}), std::out_of_range);

  int i(-10);
  math::tensor<int> int_tensor({});
  BOOST_CHECK(int_tensor.get_rank() == 0);
  BOOST_REQUIRE_THROW(int_tensor.get_size(0), std::out_of_range);
  int_tensor({}) = i;
  BOOST_CHECK(int_tensor({}) == i);
  BOOST_REQUIRE_THROW(int_tensor({5}), std::out_of_range);

  size_t s(25);
  math::tensor<size_t> szt_tensor({});
  BOOST_CHECK(szt_tensor.get_rank() == 0);
  BOOST_REQUIRE_THROW(szt_tensor.get_size(0), std::out_of_range);
  szt_tensor({}) = s;
  BOOST_CHECK(szt_tensor({}) == s);
  BOOST_REQUIRE_THROW(szt_tensor({5}), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(tensor_rank1) {
  math::tensor<double> t({2});
  BOOST_CHECK(t.get_rank() == 1);
  BOOST_CHECK(t.get_size(0) == 2);
  t({0}) = 10;
  BOOST_CHECK(t({0}) == 10);
  t({1}) = 11;
  BOOST_CHECK(t({1}) == 11);
  BOOST_REQUIRE_THROW(t({2}), std::out_of_range);
  BOOST_REQUIRE_THROW(t({1, 5}), std::out_of_range);
  BOOST_REQUIRE_THROW(t({}), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(tensor_rank2) {
  math::tensor<size_t> t({2, 3});
  BOOST_CHECK(t.get_rank() == 2);
  BOOST_CHECK(t.get_size(0) == 2);
  BOOST_CHECK(t.get_size(1) == 3);
  BOOST_REQUIRE_THROW(t({0, 3}), std::out_of_range);

  size_t input(10);
  for(size_t pos0(0); pos0 < t.get_size(0); ++pos0) {
    for(size_t pos1(0); pos1 < t.get_size(1); ++pos1, ++input) {
      t({pos0, pos1}) = input;
      BOOST_CHECK(t({pos0, pos1}) == input);
    } // for
  } // for
}

BOOST_AUTO_TEST_CASE(tensor_rank3) {
  math::tensor<size_t> t({2, 3, 4});
  BOOST_CHECK(t.get_rank() == 3);
  BOOST_CHECK(t.get_size(0) == 2);
  BOOST_CHECK(t.get_size(1) == 3);
  BOOST_CHECK(t.get_size(2) == 4);

  size_t input(10);
  for(size_t pos0(0); pos0 < t.get_size(0); ++pos0) {
    for(size_t pos1(0); pos1 < t.get_size(1); ++pos1) {
      for(size_t pos2(0); pos2 < t.get_size(2); ++pos2, ++input) {
        t({pos0, pos1, pos2}) = input;
        BOOST_CHECK(t({pos0, pos1, pos2}) == input);
      } // for
    } // for
  } // for
}

BOOST_AUTO_TEST_CASE(tensor_rank5) {
  math::tensor<size_t> t({2, 2, 2, 2, 2});
  tools::incrementor<std::vector<size_t>> inc({{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}});
  std::vector<size_t> pos(5, 0);
  size_t input(10);
  bool done(false);
  while(!done) {
    t(pos) = input;
    BOOST_CHECK(t(pos) == input);

    ++input;
    try {
      pos = inc.next(pos);
    } catch(std::overflow_error &e) {
      done = true;
    }
  } // while
}

BOOST_AUTO_TEST_CASE(tensor_sub) {
  math::tensor<size_t> t({3, 3, 3});
  size_t input(10);
  for(size_t pos0(0); pos0 < t.get_size(0); ++pos0) {
    for(size_t pos1(0); pos1 < t.get_size(1); ++pos1) {
      for(size_t pos2(0); pos2 < t.get_size(2); ++pos2, ++input) {
        t({pos0, pos1, pos2}) = input;
        BOOST_CHECK(t({pos0, pos1, pos2}) == input);
      } // for
    } // for
  } // for

  math::tensor<size_t> subt(t.sub({}, {0, 0, 0}));
  BOOST_CHECK(subt.get_rank() == 0);
  BOOST_CHECK(subt({}) == 10);
  BOOST_REQUIRE_THROW(t.sub({}, {0, 1, 3}), std::out_of_range);
  BOOST_REQUIRE_THROW(t.sub({}, {0, 1, 2, 4}), std::out_of_range);

  subt = t.sub({0}, {1, 2});
  BOOST_CHECK(subt.get_rank() == 1);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt({0}) == 15);
  BOOST_CHECK(subt({1}) == 24);
  BOOST_CHECK(subt({2}) == 33);
  subt = t.sub({1}, {1, 2});
  BOOST_CHECK(subt.get_rank() == 1);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt({0}) == 21);
  BOOST_CHECK(subt({1}) == 24);
  BOOST_CHECK(subt({2}) == 27);
  subt = t.sub({2}, {1, 2});
  BOOST_CHECK(subt.get_rank() == 1);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt({0}) == 25);
  BOOST_CHECK(subt({1}) == 26);
  BOOST_CHECK(subt({2}) == 27);

  subt = t.sub({0, 1}, {0});
  BOOST_CHECK(subt.get_rank() == 2);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt.get_size(1) == 3);
  BOOST_CHECK(subt({0, 0}) == 10);
  BOOST_CHECK(subt({0, 1}) == 13);
  BOOST_CHECK(subt({0, 2}) == 16);
  BOOST_CHECK(subt({1, 0}) == 19);
  BOOST_CHECK(subt({2, 0}) == 28);
  BOOST_CHECK(subt({2, 2}) == 34);
  subt = t.sub({0, 2}, {0});
  BOOST_CHECK(subt.get_rank() == 2);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt.get_size(1) == 3);
  BOOST_CHECK(subt({0, 0}) == 10);
  BOOST_CHECK(subt({0, 1}) == 11);
  BOOST_CHECK(subt({0, 2}) == 12);
  BOOST_CHECK(subt({1, 0}) == 19);
  BOOST_CHECK(subt({2, 0}) == 28);
  BOOST_CHECK(subt({2, 2}) == 30);
  subt = t.sub({1, 2}, {0});
  BOOST_CHECK(subt.get_rank() == 2);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt.get_size(1) == 3);
  BOOST_CHECK(subt({0, 0}) == 10);
  BOOST_CHECK(subt({0, 1}) == 11);
  BOOST_CHECK(subt({0, 2}) == 12);
  BOOST_CHECK(subt({1, 0}) == 13);
  BOOST_CHECK(subt({2, 0}) == 16);
  BOOST_CHECK(subt({2, 2}) == 18);
  subt = t.sub({0, 1, 2}, {});
  BOOST_CHECK(subt.get_rank() == 3);
  BOOST_CHECK(subt.get_size(0) == 3);
  BOOST_CHECK(subt.get_size(1) == 3);
  BOOST_CHECK(subt.get_size(2) == 3);
  for(size_t pos0(0); pos0 < t.get_size(0); ++pos0) {
    for(size_t pos1(0); pos1 < t.get_size(1); ++pos1) {
      for(size_t pos2(0); pos2 < t.get_size(2); ++pos2, ++input) {
        BOOST_CHECK(t({pos0, pos1, pos2}) == subt({pos0, pos1, pos2}));
      } // for
    } // for
  } // for
}

BOOST_AUTO_TEST_SUITE_END()
