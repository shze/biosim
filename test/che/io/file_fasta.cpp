#include <boost/test/unit_test.hpp>

#include "che/io/file_fasta.h" // header to test

using namespace biosim;

BOOST_AUTO_TEST_SUITE(suite_file_fasta)

BOOST_AUTO_TEST_CASE(file_fasta_read) {
  BOOST_CHECK(che::io::file_fasta::read("../test/data/T0666.psipred_ss2").get_chain_id_list().size() == 0);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/1a00_a.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/1a00_a_no_id.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/1a00_b.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P01236.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P01236-multifasta2.fasta").get_chain_id_list().size() == 2);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P01236-multifasta.fasta").get_chain_id_list().size() == 2);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P01236-wikipedia.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P01942.fasta").get_chain_id_list().size() == 1);
  BOOST_CHECK(che::io::file_fasta::read("../test/data/P02088.fasta").get_chain_id_list().size() == 1);

  che::assembly a(che::io::file_fasta::read("../test/data/1UBI_A_DSSP.fasta"));
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.has_molecule("A"));
  BOOST_CHECK(a.has_molecule("B"));
  BOOST_CHECK(a.has_ss("A"));
  BOOST_CHECK(a.has_ss("B"));
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_identifier() == "MET"); // first dssp sequence is matched to cc seq
  BOOST_CHECK(a.get_molecule("B").get_ps()[0].get_identifier() == "UNK"); // second dssp sequence has no match

  a = che::io::file_fasta::read("../test/data/1UBI_A_DSSP2.fasta");
  BOOST_CHECK(a.get_chain_id_list().size() == 1);
  BOOST_CHECK(a.has_molecule("A"));
  BOOST_CHECK(a.has_ss("A"));
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_identifier() == "MET");

  a = che::io::file_fasta::read("../test/data/1UBI_A_DSSP3.fasta");
  BOOST_CHECK(a.get_chain_id_list().size() == 2);
  BOOST_CHECK(a.has_molecule("A"));
  BOOST_CHECK(a.has_molecule("B"));
  BOOST_CHECK(a.has_ss("A"));
  BOOST_CHECK(a.has_ss("B"));
  BOOST_CHECK(a.get_molecule("A").get_ps()[0].get_identifier() == "UNK");
  BOOST_CHECK(a.get_molecule("B").get_ps()[0].get_identifier() == "UNK");
}

BOOST_AUTO_TEST_SUITE_END()
