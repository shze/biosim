#include <boost/test/unit_test.hpp>

#include "tools/file.h" // header to test
#include <boost/algorithm/string.hpp>

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file)

BOOST_AUTO_TEST_CASE(file_read) {
  BOOST_REQUIRE_THROW(tools::file::read_to_string_list("file_does_not_exist"), std::ios_base::failure);

  std::string const expected_content =
      "MHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPYTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDN\n"
      "LKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
  std::string const file_name = "../test/data/1a00_a_no_id.fasta"; // assumes unit_tester is starter from build folder
  std::string const actual_content(boost::algorithm::join(tools::file::read_to_string_list(file_name), "\n"));
  BOOST_CHECK(actual_content == expected_content);
}

BOOST_AUTO_TEST_SUITE_END()
