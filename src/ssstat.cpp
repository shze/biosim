#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "che/io/file_assembly.h"
#include "score/cm.h"

using namespace biosim;

int main(int argc, char *argv[]) { // main must not be in a namespace
  // declare all possible program options
  boost::program_options::options_description po_desc_main("Main options"); // the non-hidden options
  po_desc_main.add_options()("ref,r", boost::program_options::value<std::string>(),
                             "reference secondary structure file")(
      "ss,s", boost::program_options::value<std::vector<std::string>>(),
      "secondary structure file")("help,h", "print help")("debug,d", "print debug output");
  // declare positional options
  boost::program_options::positional_options_description po_desc_pos;
  po_desc_pos.add("ref", -1);

  // declare map to store the values, parse command line and store values
  boost::program_options::variables_map po_var_map;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).options(po_desc_main).positional(po_desc_pos).run(),
      po_var_map);
  boost::program_options::notify(po_var_map);

  // print help if requested by the user, or if no reference or no other secondary structure definitions given, and exit
  if(po_var_map.count("help") || (po_var_map.count("ref") == 0 || po_var_map.count("ss") == 0)) {
    LOG << "This program compares 1..n ss definitions to one reference/native ss by calculating different measures.";
    LOG << "Supported file types are .fasta, .psipred_ss, .psipred_ss2, .rdbProf, .jufo9d_ss, and .pdb;";
    LOG << "Only for non-reference ss definitions, .pool files are supported additionally.";
    LOG << po_desc_main << "\n";
    return 0;
  } // if

  // set requested message level
  tools::log_message::set_output_level(po_var_map.count("debug") ? tools::log_message::level_debug
                                                                 : tools::log_message::level_info);

  try {
    che::io::file_assembly reader;
    // read reference
    LOG << "Reading reference ss from " << po_var_map["ref"].as<std::string>();
    che::assembly ref = reader.read(po_var_map["ref"].as<std::string>());
    DEBUG << ref;
    reader.add_sse_pool_reader(ref); // add reader for .pool files with reference
    // read other ss definitions
    std::vector<che::assembly> ss;
    for(auto filename : po_var_map["ss"].as<std::vector<std::string>>()) {
      LOG << "Reading ss from " << filename;
      ss.emplace_back(reader.read(filename));
      DEBUG << ss.back();
    } // for

    for(auto a : ss) {
      for(auto p : score::assembly_compares::get_instances()) {
        LOG << p.first << "=" << p.second.get_object()->compare(ref, a);
      } // for
    } // for
  } catch(std::exception &e) {
    ERROR << e.what();
  } // catch

  return 0;
} // main()
