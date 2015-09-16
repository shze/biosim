#include <boost/test/unit_test.hpp>

#include "tools/mapper.h" // header to test
#include <map>

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_mapper)

BOOST_AUTO_TEST_CASE(mapper_ctor) {
  tools::mapper<std::string> m({});
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 0);
  BOOST_CHECK(m.encode(0) == std::string());
  BOOST_CHECK(m.encode(1) == std::string());

  BOOST_REQUIRE_THROW(tools::mapper<std::string>({{}}), std::invalid_argument);

  m = tools::mapper<std::string>({{'0', '1'}});
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 1);

  m = tools::mapper<std::string>({{'0', '1'}, {'0', '1'}});
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 3);

  m = tools::mapper<std::string>({{'0', '1'}, {'0', '1'}, {'0', '1'}});
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 7);
}

BOOST_AUTO_TEST_CASE(mapper_string) {
  std::map<size_t, std::string> expected_big_endian{{0, "000"},
                                                    {1, "001"},
                                                    {2, "002"},
                                                    {3, "010"},
                                                    {4, "011"},
                                                    {5, "012"},
                                                    {6, "100"},
                                                    {7, "101"},
                                                    {8, "102"},
                                                    {9, "110"},
                                                    {10, "111"},
                                                    {11, "112"},
                                                    {12, "112"}};
  tools::mapper<std::string> m({{'0', '1'}, {'0', '1'}, {'0', '1', '2'}}, true);
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 11);
  for(size_t index(m.get_min()); index <= m.get_max(); ++index) {
    BOOST_CHECK(m.encode(index) == expected_big_endian[index]);
    BOOST_CHECK(m.decode(m.encode(index)) == index);
  } // for

  std::map<size_t, std::string> expected_little_endian{{0, "000"},
                                                       {1, "100"},
                                                       {2, "010"},
                                                       {3, "110"},
                                                       {4, "001"},
                                                       {5, "101"},
                                                       {6, "011"},
                                                       {7, "111"},
                                                       {8, "002"},
                                                       {9, "102"},
                                                       {10, "012"},
                                                       {11, "112"},
                                                       {12, "112"}};
  m = tools::mapper<std::string>({{'0', '1'}, {'0', '1'}, {'0', '1', '2'}}, false);
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 11);
  for(size_t index(m.get_min()); index <= m.get_max(); ++index) {
    BOOST_CHECK(m.encode(index) == expected_little_endian[index]);
    BOOST_CHECK(m.decode(m.encode(index)) == index);
  } // for
}

BOOST_AUTO_TEST_CASE(mapper_vector) {
  std::map<size_t, std::vector<size_t>> expected_big_endian{{0, {0, 0, 0}},
                                                            {1, {0, 0, 1}},
                                                            {2, {0, 0, 2}},
                                                            {3, {0, 1, 0}},
                                                            {4, {0, 1, 1}},
                                                            {5, {0, 1, 2}},
                                                            {6, {1, 0, 0}},
                                                            {7, {1, 0, 1}},
                                                            {8, {1, 0, 2}},
                                                            {9, {1, 1, 0}},
                                                            {10, {1, 1, 1}},
                                                            {11, {1, 1, 2}},
                                                            {12, {1, 1, 2}}};
  tools::mapper<std::vector<size_t>> m({{0, 1}, {0, 1}, {0, 1, 2}}, true);
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 11);
  for(size_t index(m.get_min()); index <= m.get_max(); ++index) {
    BOOST_CHECK(m.encode(index) == expected_big_endian[index]);
    BOOST_CHECK(m.decode(m.encode(index)) == index);
  } // for

  std::map<size_t, std::vector<size_t>> expected_little_endian{{0, {0, 0, 0}},
                                                               {1, {1, 0, 0}},
                                                               {2, {0, 1, 0}},
                                                               {3, {1, 1, 0}},
                                                               {4, {0, 0, 1}},
                                                               {5, {1, 0, 1}},
                                                               {6, {0, 1, 1}},
                                                               {7, {1, 1, 1}},
                                                               {8, {0, 0, 2}},
                                                               {9, {1, 0, 2}},
                                                               {10, {0, 1, 2}},
                                                               {11, {1, 1, 2}},
                                                               {12, {1, 1, 2}}};
  m = tools::mapper<std::vector<size_t>>({{0, 1}, {0, 1}, {0, 1, 2}}, false);
  BOOST_CHECK(m.get_min() == 0);
  BOOST_CHECK(m.get_max() == 11);
  for(size_t index(m.get_min()); index <= m.get_max(); ++index) {
    BOOST_CHECK(m.encode(index) == expected_little_endian[index]);
    BOOST_CHECK(m.decode(m.encode(index)) == index);
  } // for
}

BOOST_AUTO_TEST_SUITE_END()
