#include "che/io/file_sse_pool.h"
#include "tools/file.h"
#include "tools/log.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // stores all chain-dependent information that a pool file contains
      struct chain_data {
        // default ctor
        chain_data()
            : last_helix_serial_number(0),
              last_seqres_serial_number(0),
              last_seqres_residue_count(0),
              cc_sequence(),
              pool() {}

        size_t last_helix_serial_number; // starts with 1, and will be increased by 1 for each helix
        size_t last_seqres_serial_number; // starts with 1, and will be increased by 1 for each SEQRES lines
        size_t last_seqres_residue_count; // 0=not set
        ps cc_sequence; // amino acid sequence
        std::set<cchb_dssp_interval> pool; // secondary structure definition lines
      }; // struct chain_data

      // read ps and ss from a given file
      assembly file_sse_pool::read(const std::string &__filename) { return read(__filename, assembly()); }
      // read ps and ss from a given file and extends the ps to the given reference length
      assembly file_sse_pool::read(const std::string &__filename, assembly const &__reference) {
        // regex for a single data line of the sse pool format; groups will be used to construct sequences;
        // this regex will not match the begin ('bcl::assemble::SSEPool') and end lines ('END') of a sse pool file to
        // allow reading pdb files as well;
        // for the line specification, see: http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        // and http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html;
        // all pdb lines have to be 80 char long, see http://deposit.rcsb.org/adit/docs/pdb_atom_format.html;
        // read_to_string_list() trims the whitespace at begin and end, so if the line length is less than 80 chars,
        // a line will be extended to that length; if the line has more whitespace than 80 chars, it will be ignored.
        // HELIX line groups explanation:
        // 1=HELIX; 2=Helix serial number; 3=Helix id; 4=Initial residue name; 5=Chain id; 6=Residue sequence number;
        // 7=Code for insertions of residues; 8=Terminal residue name; 9=Chain id; 10=Residue sequence number;
        // 11=Code for insertions of residues; 12=Helix type; 13=Comment; 14=Helix length.
        // SHEET line groups explanation:
        // 15=SHEET; 16=Strand number in current sheet; 17=Sheet id; 18=Number of strands in current sheet;
        // 19=Initial residue name; 20=Chain id; 21=Residue sequence number; 22=Code for insertions of residues;
        // 23=Terminal residue name; 24=Chain id; 25=Residue sequence number; 26=Code for insertions of residues;
        // 27=Strand sense with respect to previous;
        // there is additional information on a SHEET line that is not in explicit groups, but in a generic
        // placeholder at the end of the SHEET part of the regex.
        static boost::regex const regex_sse_pool_line(
            // 1      2           3     4     5   6          7   8     9   10         11 12         13      14
            "(HELIX)  ([0-9 ]{3}) (...) (...) (.) ([0-9 ]{4})(.) (...) (.) ([0-9 ]{4})(.)([0-9 ]{2})(.{30}) ([0-9 ]{5})"
            ".*|" // match any space that go beyond the 76 chars specified in the regex, or match the SHEET line
            // 15     16          17   18          19    20 21         22  23    24 25         26 27
            "(SHEET)  ([0-9 ]{3}) (...)([0-9 ]{2}) (...) (.)([0-9 ]{4})(.) (...) (.)([0-9 ]{4})(.)([0-9 ]{2}).*|");

        std::list<std::string> file_lines(tools::file::read_to_string_list(__filename));
        DEBUG << "Read sse pool file content:\n" << boost::algorithm::join(file_lines, "\n");

        std::map<char, chain_data> chains; // local variable to store the data read from the file

        DEBUG << "Parsing sse pool file";
        for(auto &line : file_lines) { // process each line from the file
          // make all lines have length >=80 by filling up with spaces if length <80, to ensure the regex will match
          line.append(std::max(80 - line.length(), (size_t)0), ' ');

          boost::smatch match;
          // regex_match returns true if all chars of line were matched
          if(!boost::regex_match(line, match, regex_sse_pool_line)) {
            DEBUG << "No match: " << line;
            continue;
          } // if

          // match[0], same as match.str(), contains the whole string, see
          // http://www.boost.org/doc/libs/1_55_0/libs/regex/example/snippets/regex_iterator_example.cpp
          std::string complete_match(match[0]);
          boost::trim(complete_match); // remove newline at the end (the regex above matches an optional newline)
          DEBUG << "Match: " << complete_match;

          if(match[1].matched) { // if the 'HELIX' part matches, a helix line was found
            DEBUG << "Found HELIX line";

            // get strings from all submatches, then check if contents are valid for conversion;
            // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
            // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
            // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
            std::string const helix_serial_number_string(boost::trim_copy(std::string(match[2])));
            std::string const helix_id(boost::trim_copy(std::string(match[3])));
            std::string const initial_residue_name(match[4]);
            char const initial_residue_chain_id(boost::lexical_cast<char>(match[5]));
            size_t const initial_residue_sequence_number(
                boost::lexical_cast<size_t>(boost::trim_copy(std::string(match[6]))));
            char const initial_residue_insertion_code(boost::lexical_cast<char>(match[7]));
            std::string const terminal_residue_name(match[8]);
            char const terminal_residue_chain_id(boost::lexical_cast<char>(match[9]));
            size_t const terminal_residue_sequence_number(
                boost::lexical_cast<size_t>(boost::trim_copy(std::string(match[10]))));
            char const terminal_residue_insertion_code(boost::lexical_cast<char>(match[11]));
            std::string const helix_type_string(boost::trim_copy(std::string(match[12])));
            std::string const comment(match[13]);
            std::string const helix_length_string(boost::trim_copy(std::string(match[14])));

            // print the submatches as strings before consistency checks and converting to size_t
            DEBUG << "Submatches: helix_serial_number='" << helix_serial_number_string << "' helix_id='" << helix_id
                  << "' initial_residue='" << initial_residue_name << "|" << initial_residue_chain_id << "|"
                  << initial_residue_sequence_number << "|" << initial_residue_insertion_code << "' terminal_residue='"
                  << terminal_residue_name << "|" << terminal_residue_chain_id << "|"
                  << terminal_residue_sequence_number << "|" << terminal_residue_insertion_code << "' helix_type='"
                  << helix_type_string << "' comment='" << boost::trim_copy(comment) << "' helix_length='"
                  << helix_length_string << "'";

            // both chain ids should be identical
            if(initial_residue_chain_id != terminal_residue_chain_id) {
              DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                    << "' differing from terminal_residue_chain_id='" << terminal_residue_chain_id << "'";
              continue;
            } // if
            // add chain id if it isn't found in the map
            if(chains.find(initial_residue_chain_id) == chains.end()) {
              chains[initial_residue_chain_id] = chain_data();
            } // if

            // convert strings to size_t and catch exceptions for ones which have a default or a way of calculating
            size_t helix_serial_number;
            try {
              helix_serial_number = boost::lexical_cast<size_t>(helix_serial_number_string);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              helix_serial_number = chains[initial_residue_chain_id].last_helix_serial_number + 1;
              DEBUG << "Could not convert helix_serial_number='" << helix_serial_number_string
                    << "' to number, using default helix_serial_number='" << helix_serial_number << "'";
            } // catch
            // helix serial number should be incremented by 1 from last helix serial number
            if(helix_serial_number != chains[initial_residue_chain_id].last_helix_serial_number + 1) {
              DEBUG << "Found helix_serial_number='" << helix_serial_number
                    << "' differing from expected helix_serial_number='"
                    << (chains[initial_residue_chain_id].last_helix_serial_number + 1)
                    << "', using helix_serial_number from file.";
            } // if
            // update the helix counting variable
            chains[initial_residue_chain_id].last_helix_serial_number = helix_serial_number;

            // we could convert helix_type_string into size_t, but realistically no checks can be done b/c a pool
            // often comes from prediction and the helix_type is not known

            size_t helix_length;
            try {
              helix_length = boost::lexical_cast<size_t>(helix_length_string);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              // + 1, b/c the initial_residue_sequence_number is part of the helix
              helix_length = terminal_residue_sequence_number - initial_residue_sequence_number + 1;
              DEBUG << "Could not convert helix_length='" << helix_length_string
                    << "' to number, using calculated helix_length='" << helix_length << "'";
            } // catch
            // check helix_length is consist with residue sequence numbers
            if(helix_length != terminal_residue_sequence_number - initial_residue_sequence_number + 1) {
              size_t const expected_helix_length(terminal_residue_sequence_number - initial_residue_sequence_number +
                                                 1);
              DEBUG << "Found helix_length='" << helix_length << "' differing expected helix_length='"
                    << expected_helix_length << "', ignoring helix_length.";
              helix_length = expected_helix_length;
            } // if

            // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
            chains[initial_residue_chain_id].cc_sequence.resize(
                std::max(chains[initial_residue_chain_id].cc_sequence.size(), terminal_residue_sequence_number),
                cc(cc::specificity_type::unknown, cc::monomer_type::l_peptide_linking));
            chains[initial_residue_chain_id].cc_sequence[initial_residue_sequence_number - 1] =
                cc(initial_residue_name);
            chains[initial_residue_chain_id].cc_sequence[terminal_residue_sequence_number - 1] =
                cc(terminal_residue_name);

            // create and insert sequence_interval into the pool
            chains[initial_residue_chain_id].pool.insert(cchb_dssp_interval(
                initial_residue_sequence_number - 1, terminal_residue_sequence_number - 1, cchb_dssp('H')));
          } // if
          else if(match[15].matched) { // if the 'SHEET' part matches, a sheet line was found
            DEBUG << "Found SHEET line";

            // get strings from all submatches, then check if contents are valid for conversion;
            // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
            // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
            // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
            std::string const strand_number_in_sheet_string(boost::trim_copy(std::string(match[16])));
            std::string const sheet_id(boost::trim_copy(std::string(match[17])));
            std::string const number_strands_in_sheet_string(boost::trim_copy(std::string(match[18])));
            std::string const initial_residue_name(match[19]);
            char const initial_residue_chain_id(boost::lexical_cast<char>(match[20]));
            size_t const initial_residue_sequence_number(
                boost::lexical_cast<size_t>(boost::trim_copy(std::string(match[21]))));
            char const initial_residue_insertion_code(boost::lexical_cast<char>(match[22]));
            std::string const terminal_residue_name(match[23]);
            char const terminal_residue_chain_id(boost::lexical_cast<char>(match[24]));
            size_t const terminal_residue_sequence_number(
                boost::lexical_cast<size_t>(boost::trim_copy(std::string(match[25]))));
            char const terminal_residue_insertion_code(boost::lexical_cast<char>(match[26]));
            std::string const strand_sense_string(boost::trim_copy(std::string(match[27])));

            // print the submatches as strings before consistency checks and converting to size_t
            DEBUG << "Submatches: strand_number_in_sheet='" << strand_number_in_sheet_string << "' sheet_id='"
                  << sheet_id << "' number_strands_in_sheet='" << number_strands_in_sheet_string
                  << "' initial_residue='" << initial_residue_name << "|" << initial_residue_chain_id << "|"
                  << initial_residue_sequence_number << "|" << initial_residue_insertion_code << "' terminal_residue='"
                  << terminal_residue_name << "|" << terminal_residue_chain_id << "|"
                  << terminal_residue_sequence_number << "|" << terminal_residue_insertion_code << "' strand_sense='"
                  << strand_sense_string << "'";

            // we could try to convert strand_number_in_sheet_string, number_strands_in_sheet_string, and
            // strand_sense_string into size_t, but no checks could be done realistically b/c the pool likely comes
            // from secondary structure prediction from sequence, and won't have this information.

            // both chain ids should be identical; no fix needed, chain id is always ignored.
            if(initial_residue_chain_id != terminal_residue_chain_id) {
              DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                    << "' differing from terminal_residue_chain_id='" << terminal_residue_chain_id << "'";
              continue;
            } // if
            // add chain id if it isn't found in the map
            if(chains.find(initial_residue_chain_id) == chains.end()) {
              chains[initial_residue_chain_id] = chain_data();
            } // if

            // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
            chains[initial_residue_chain_id].cc_sequence.resize(
                std::max(chains[initial_residue_chain_id].cc_sequence.size(), terminal_residue_sequence_number),
                cc(cc::specificity_type::unknown, cc::monomer_type::l_peptide_linking));
            chains[initial_residue_chain_id].cc_sequence[initial_residue_sequence_number - 1] =
                cc(initial_residue_name);
            chains[initial_residue_chain_id].cc_sequence[terminal_residue_sequence_number - 1] =
                cc(terminal_residue_name);

            // create and insert sequence_interval into the pool
            chains[initial_residue_chain_id].pool.insert(cchb_dssp_interval(
                initial_residue_sequence_number - 1, terminal_residue_sequence_number - 1, cchb_dssp('E')));
          } // else if
        } // for

        assembly a; // final result

        size_t sequence_no(1); // start with 1
        for(auto chain_pair : chains) {
          // extend the sequences to the reference length
          std::string chain_id_string(1, chain_pair.first);
          if(__reference.has_ss(chain_id_string)) {
            chain_pair.second.cc_sequence.resize(
                std::max(chain_pair.second.cc_sequence.size(),
                         __reference.get_ss(chain_id_string).get_sequence().size()),
                cc(cc::specificity_type::unknown, cc::monomer_type::l_peptide_linking));
          } // if

          // add to final result
          a.set(chain_id_string, molecule(__filename + "/" + std::to_string(sequence_no), ">lcl|sequence",
                                          chain_pair.second.cc_sequence),
                ss(chain_pair.second.pool, chain_pair.second.cc_sequence.size()));

          ++sequence_no; // increase sequence_no, b/c it's not done in the for loop header
        } // for

        return a;
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
