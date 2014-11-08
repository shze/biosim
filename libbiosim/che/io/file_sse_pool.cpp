#include "che/io/file_sse_pool.h"
#include "tools/file.h"
#include "tools/log.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // creates a quarternary structure from a given file
      qs file_sse_pool::read(const std::string &__filename) { return read(__filename, qs()); }
      // creates a quarternary structure from a given file and extends the sequence to the given reference length
      qs file_sse_pool::read(const std::string &__filename, qs const &__reference) {
        // regex for a single data line of the sse pool format; groups will be used to construct sequences;
        // this regex will not match the begin ('bcl::assemble::SSEPool') and end lines ('END') of a sse pool file to
        // allow reading pdb files as well;
        // for the line specification, see: http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        // and http://www.wwpdb.org/documentation/format23/sect3.html;
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
        // SEQRES line groups explanation:
        // 28=SEQRES; 29=Serial number of SEQRES record for current chain; 30=Chain id; 31=Number of residues in chain;
        // 31+N=Residue name N for N=1..13;
        static boost::regex const regex_sse_pool_line(
            // 1      2           3     4     5   6          7   8     9   10         11 12         13      14
            "(HELIX)  ([0-9 ]{3}) (...) (...) (.) ([0-9 ]{4})(.) (...) (.) ([0-9 ]{4})(.)([0-9 ]{2})(.{30}) ([0-9 ]{5})"
            ".*|" // match any space that go beyond the 76 chars specified in the regex, or match the SHEET line
            // 15     16          17   18          19    20 21         22  23    24 25         26 27
            "(SHEET)  ([0-9 ]{3}) (...)([0-9 ]{2}) (...) (.)([0-9 ]{4})(.) (...) (.)([0-9 ]{4})(.)([0-9 ]{2}).*|"
            // 28      29          30  31           32(1) 33(2) 34(3) 35(4) 36(5) 37(6) 38(7) 39(8) 40(9) 41(10)42(11)
            "(SEQRES)  ([0-9 ]{2}) (.) ([0-9 ]{4})  (...) (...) (...) (...) (...) (...) (...) (...) (...) (...) (...) "
            //   43(12)44(13)
            "(...) (...).*");

        std::list<std::string> file_lines(tools::file::read_to_string_list(__filename));
        DEBUG << "Read sse pool file content:\n" << boost::algorithm::join(file_lines, "\n");

        // create variables to store the data read from the file
        ps cc_sequence;
        std::set<cchb_dssp_interval> pool;

        // variables to keep track of the numbers of helices; don't track strand# b/c it may be unknown when predicted
        size_t last_helix_serial_number(0); // starts with 1, and will be increased by 1 for each helix
        size_t last_seqres_serial_number(0); // starts with 1, and will be increased by 1 for each SEQRES lines
        size_t last_seqres_residue_count(0); // 0=not set
        char chain_id;
        bool chain_id_set(false);

        DEBUG << "Parsing sse pool file";
        for(auto &line : file_lines) { // process each line from the file
          // make all lines have length >=80 by filling up with spaces if length <80, to ensure the regex will match
          line.append(std::max(80 - line.length(), (size_t)0), ' ');

          boost::smatch match;
          // regex_match returns true if all chars of line were matched;
          // if that is true, the match size should be 45 (all match + groups), i.e. the second part is likely redundant
          if(!boost::regex_match(line, match, regex_sse_pool_line) || match.size() != 45) {
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

            // convert strings to size_t and catch exceptions for ones which have a default or way the calculating
            size_t helix_serial_number;
            try {
              helix_serial_number = boost::lexical_cast<size_t>(helix_serial_number_string);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              helix_serial_number = last_helix_serial_number + 1;
              DEBUG << "Could not convert helix_serial_number='" << helix_serial_number_string
                    << "' to number, using default helix_serial_number='" << helix_serial_number << "'";
            } // catch

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

            // helix serial number should be incremented by 1 from last helix serial number
            if(helix_serial_number != last_helix_serial_number + 1) {
              DEBUG << "Found helix_serial_number='" << helix_serial_number
                    << "' differing from expected helix_serial_number='" << (last_helix_serial_number + 1)
                    << "', using helix_serial_number from file.";
            } // if
            // update the helix counting variable
            last_helix_serial_number = helix_serial_number;

            // both chain ids should be identical
            if(initial_residue_chain_id != terminal_residue_chain_id) {
              DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                    << "' differing from terminal_residue_chain_id='" << terminal_residue_chain_id << "'";
              continue;
            } // if

            // check helix_length is consist with residue sequence numbers
            if(helix_length != terminal_residue_sequence_number - initial_residue_sequence_number + 1) {
              size_t const expected_helix_length(terminal_residue_sequence_number - initial_residue_sequence_number +
                                                 1);
              DEBUG << "Found helix_length='" << helix_length << "' differing expected helix_length='"
                    << expected_helix_length << "', ignoring helix_length.";
              helix_length = expected_helix_length;
            } // if

            if(chain_id_set) {
              if(initial_residue_chain_id != chain_id) {
                DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                      << "' differing from previous chain_id='" << chain_id << "'";
                continue;
              } // if
            } // if
            else {
              chain_id = initial_residue_chain_id;
              chain_id_set = true;
            } // else

            // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
            cc_sequence.resize(std::max(cc_sequence.size(), terminal_residue_sequence_number),
                               cc(cc::specificity_unknown, cc::l_peptide_linking));
            cc_sequence[initial_residue_sequence_number - 1] = cc(initial_residue_name);
            cc_sequence[terminal_residue_sequence_number - 1] = cc(terminal_residue_name);

            // create and insert sequence_interval into the pool
            pool.insert(cchb_dssp_interval(initial_residue_sequence_number - 1, terminal_residue_sequence_number - 1,
                                           cchb_dssp('H')));
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

            if(chain_id_set) {
              if(initial_residue_chain_id != chain_id) {
                DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                      << "' differing from previous chain_id='" << chain_id << "'";
                continue;
              } // if
            } // if
            else {
              chain_id = initial_residue_chain_id;
              chain_id_set = true;
            } // else

            // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
            cc_sequence.resize(std::max(cc_sequence.size(), terminal_residue_sequence_number),
                               cc(cc::specificity_unknown, cc::l_peptide_linking));
            cc_sequence[initial_residue_sequence_number - 1] = cc(initial_residue_name);
            cc_sequence[terminal_residue_sequence_number - 1] = cc(terminal_residue_name);

            // create and insert sequence_interval into the pool
            pool.insert(cchb_dssp_interval(initial_residue_sequence_number - 1, terminal_residue_sequence_number - 1,
                                           cchb_dssp('E')));
          } // else if
          else if(match[28].matched) { // if the 'SEQRES' part matches, a sequence residue line was found
            DEBUG << "Found SEQRES line";

            // get strings from all submatches except the residues
            std::string const seqres_serial_number_string(boost::trim_copy(std::string(match[29])));
            char const seqres_chain_id(boost::lexical_cast<char>(match[30]));
            std::string const seqres_residue_count_string(boost::trim_copy(std::string(match[31])));

            // print the submatches as strings before consistency checks and converting to size_t
            DEBUG << "Submatches: seqres_serial_number='" << seqres_serial_number_string << "' seqres_chain_id='"
                  << seqres_chain_id << "' seqres_residue_count='" << seqres_residue_count_string << "' residue_names='"
                  << match[32] << "|" << match[33] << "|" << match[34] << "|" << match[35] << "|" << match[36] << "|"
                  << match[37] << "|" << match[38] << "|" << match[39] << "|" << match[40] << "|" << match[41] << "|"
                  << match[42] << "|" << match[43] << "|" << match[44] << "'";

            size_t seqres_serial_number;
            try {
              seqres_serial_number = boost::lexical_cast<size_t>(seqres_serial_number_string);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              seqres_serial_number = last_seqres_serial_number + 1;
              DEBUG << "Could not convert seqres_serial_number='" << seqres_serial_number_string
                    << "' to number, using calculated seqres_serial_number='" << seqres_serial_number << "'";
            } // catch

            // check helix_length is consist with residue sequence numbers
            if(seqres_serial_number != last_seqres_serial_number + 1) {
              size_t const expected_seqres_serial_number(last_seqres_serial_number + 1);
              DEBUG << "Found seqres_serial_number='" << seqres_serial_number << "' differing from expected "
                    << "seqres_serial_number='" << expected_seqres_serial_number << "', ignoring seqres_serial_number.";
              seqres_serial_number = expected_seqres_serial_number;
            } // if

            last_seqres_serial_number = seqres_serial_number; // update the SEQRES serial number variable

            size_t seqres_residue_count(0);
            try {
              seqres_residue_count = boost::lexical_cast<size_t>(seqres_residue_count_string);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              // keep seqres_residue_count set to 0
              DEBUG << "Could not convert seqres_residue_count='" << seqres_residue_count_string << "' to number";
            } // catch

            if(last_seqres_residue_count != 0) {
              if(last_seqres_residue_count != seqres_residue_count) {
                // just ignore the seqres_residue_count if it's different, do not skip this line
                DEBUG << "Found seqres_residue_count='" << seqres_residue_count << "' differing from expected "
                      << "seqres_residue_count='" << last_seqres_residue_count << "', ignoring seqres_residue_count.";
              } // if
            } // if
            else {
              last_seqres_residue_count = seqres_residue_count; // if last_seqres_residue_count is still 0, save it
            }

            if(chain_id_set) {
              if(seqres_chain_id != chain_id) {
                // if a previous chain_id was set, and this one is different, ignore line, it could be a different chain
                DEBUG << "Ignoring line, seqres_chain_id='" << seqres_chain_id << "' differing from previous chain_id='"
                      << chain_id << "'";
                continue;
              } // if
            } // if
            else {
              chain_id = seqres_chain_id;
              chain_id_set = true;
            } // else

            // loop through all 13 residue name groups and check if there's a residue or spaces
            for(size_t pos(32); pos <= 44; ++pos) {
              std::string const residue_name(boost::trim_copy(std::string(match[pos])));
              if(residue_name.empty()) {
                continue; // ignore empty names, likely at the end of the SEQRES definition
              } // if

              cc_sequence.emplace_back(cc(residue_name));
            } // for
          } // else if
        } // for

        // extend the sequences to the reference length
        if(!__reference.get_chain_id_list().empty()) {
          cc_sequence.resize(std::max(cc_sequence.size(), __reference.get_ss("A").get_sequence().size()),
                             cc(cc::specificity_unknown, cc::l_peptide_linking));
        }

        return cc_sequence.empty() ? qs()
                                   : qs(ts(__filename, ">lcl|sequence", cc_sequence), ss(pool, cc_sequence.size()));
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
