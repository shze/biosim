#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "che/io/file_assembly.h"
#include "che/algo/aligner_dp.h"

using namespace biosim;

int main(int argc, char *argv[]) { // main must not be in a namespace
  // declare all possible program options
  boost::program_options::options_description po_desc_main("Main options"); // the non-hidden options
  po_desc_main.add_options()("seq,s", boost::program_options::value<std::vector<std::string>>(),
                             "sequence file")("help,h", "print help")("debug,d", "print debug output");
  // declare positional options
  boost::program_options::positional_options_description po_desc_pos;
  po_desc_pos.add("seq", -1);

  // declare map to store the values, parse command line and store values
  boost::program_options::variables_map po_var_map;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).options(po_desc_main).positional(po_desc_pos).run(),
      po_var_map);
  boost::program_options::notify(po_var_map);

  // print help if requested by the user, or if no reference or no other secondary structure definitions given, and exit
  if(po_var_map.count("help") || po_var_map.count("seq") == 0) {
    LOG << "This program compares sequences with the Blosum62 matrix, the gap score set to the minimum matrix value.";
    LOG << "Supported file types are .fasta, .psipred_ss, .psipred_ss2, .rdbProf, .jufo9d_ss, and .pdb;";
    LOG << po_desc_main << "\n";
    return 0;
  } // if

  // set requested message level
  tools::log_message::set_output_level(po_var_map.count("debug") ? tools::log_message::level_type::debug
                                                                 : tools::log_message::level_type::info);

  try {
    // read all files into molecule
    che::io::file_assembly reader;
    std::vector<che::molecule> molecules;
    for(auto filename : po_var_map["seq"].as<std::vector<std::string>>()) {
      LOG << "Reading sequence from " << filename;
      std::vector<che::molecule> new_molecules(reader.read(filename).get_molecules());
      molecules.reserve(molecules.size() + new_molecules.size());
      molecules.insert(molecules.end(), new_molecules.begin(), new_molecules.end());
      LOG << "Found " << new_molecules.size() << " molecules";
      for(auto m : new_molecules) {
        LOG << m;
      } // for
    } // for

    // convert molecule into alignment
    std::vector<che::alignment> input_alignments;
    for(auto m : molecules) {
      input_alignments.emplace_back(m);
    } // for

    che::algo::aligner_dp aligner;
    std::list<che::scored_alignment> best_alignments(aligner.align_multiple(input_alignments));
    LOG << "Found " << best_alignments.size() << " best scoring alignments:";
    for(che::scored_alignment a : best_alignments) {
      LOG << a;
    } // for
  } catch(std::exception &e) {
    ERROR << e.what();
  } // catch

  return 0;
} // main()
