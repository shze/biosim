#include "che/io/file_pdb.h"
#include "tools/file.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // stores the caveat data
      class caveat_data {
      public:
        // processes the input strings
        void process_caveat() {}
      }; // class caveat_data

      // converts string input from pdb files into internal data structures and stores it
      class seqres_data {
      private:
        // stores the seqres data for a single chain
        struct seqres_chain_data {
          size_t _last_serial_number; // last seen serial number; starts at 1 for each chain, + 1 for each SEQRES line
          size_t _last_residue_count; // last seen residue count
          ps _cc_sequence; // converted sequence data
        }; // struct seqres_chain_data

        std::map<char, seqres_chain_data> _chains; // map chain_id->chain_data

      public:
        // processes the input strings
        void process(std::string const &__serial_number_string, std::string const &__chain_id,
                     std::string const &__residue_count_string, std::vector<std::string> const &__residues) {
          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: seqres_serial_number='" << __serial_number_string << "' seqres_chain_id='" << __chain_id
                << "' seqres_residue_count='" << __residue_count_string << "' residue_names='"
                << boost::algorithm::join(__residues, "|") << "'";

          char const chain_id(boost::lexical_cast<char>(__chain_id));
          auto chain_itr(_chains.find(chain_id)); // if equal to _chains.end(), this chain_id is not in the map

          size_t serial_number;
          try {
            serial_number = boost::lexical_cast<size_t>(__serial_number_string);
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            serial_number = chain_itr == _chains.end() ? 1 : chain_itr->second._last_serial_number + 1;
            DEBUG << "Could not convert seqres_serial_number='" << __serial_number_string
                  << "' to number, using calculated seqres_serial_number='" << serial_number << "'";
          } // catch

          // check helix_length is consist with residue sequence numbers
          if(chain_itr != _chains.end() && serial_number != chain_itr->second._last_serial_number + 1) {
            size_t const expected_serial_number(chain_itr->second._last_serial_number + 1);
            DEBUG << "Found seqres_serial_number='" << serial_number << "' differing from expected "
                  << "seqres_serial_number='" << expected_serial_number << "', ignoring seqres_serial_number.";
            serial_number = expected_serial_number;
          } // if

          size_t residue_count(0); // initialize to 0, b/c it might not be set if an exception occurs
          try {
            residue_count = boost::lexical_cast<size_t>(__residue_count_string);
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            // keep seqres_residue_count set to 0
            DEBUG << "Could not convert seqres_residue_count='" << __residue_count_string << "' to number";
          } // catch

          if(chain_itr != _chains.end()) {
            size_t expected_residue_count(chain_itr->second._last_residue_count); // last residue count
            if(residue_count == 0) {
              DEBUG << "Found seqres_residue_count='" << residue_count << "' differing from expected "
                    << "seqres_residue_count='" << expected_residue_count << "', ignoring seqres_residue_count.";
              residue_count = expected_residue_count;
            } // if
            else if(residue_count != chain_itr->second._last_residue_count) {
              // just ignore the residue_count if it's different, do not skip this SEQRES line
              DEBUG << "Found seqres_residue_count='" << residue_count << "' differing from expected "
                    << "seqres_residue_count='" << expected_residue_count << "', ignoring seqres_residue_count.";
            } // else if
          } // if

          // save everything
          _chains[chain_id]._last_serial_number = serial_number;
          _chains[chain_id]._last_residue_count = residue_count;
          // loop through all 13 residue name groups and check if there's a residue or space, and save the residues
          for(auto r : __residues) {
            std::string const residue_name(boost::trim_copy(r));
            if(residue_name.empty()) {
              continue; // ignore empty names, likely at the end of the SEQRES definition
            } // if

            _chains[chain_id]._cc_sequence.emplace_back(cc(residue_name));
          } // for
        } // process()

        // get chain ids
        std::list<char> get_chain_ids() {
          std::list<char> chain_ids;
          for(auto p : _chains) {
            chain_ids.emplace_back(p.first);
          } // for
          return chain_ids;
        } // get_chain_ids()

        // returns the ps
        ps get_sequence(char const &__chain_id) { return _chains[__chain_id]._cc_sequence; } // get_sequence()
      }; // class seqres_data

      // stores the ss definition data
      class ssdef_data {
      public:
        // processes the input strings
        void process_helix() {}
        // processes the input strings
        void process_strand() {}
      }; // class ssdef_data

      // stores the model data
      class model_data {
      public:
        // processes the input strings
        void process_number_models() {}
        // processes the input strings
        void process_model_begin() {}
        // processes the input strings
        void process_model_end() {}
        // processes the input strings
        void process_atom() {}
        // processes the input strings
        void process_hetatm() {}
        // processes the input strings
        void process_terminate() {}
      }; // class model_data

      // read ps and ss from a given file
      assembly file_pdb::read(const std::string &__filename) {
        // regex for a single data line of the pdb format; groups will be used to construct sequences;
        // for the line specification, see: http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        // and http://www.rcsb.org/pdb/static.do?p=file_formats/pdb/index.html;
        // all pdb lines have to be 80 char long, see http://deposit.rcsb.org/adit/docs/pdb_atom_format.html;
        // read_to_string_list() trims the whitespace at begin and end, so if the line length is less than 80 chars,
        // a line will be extended to that length; if the line has more whitespace than 80 chars, it will be ignored.
        static boost::regex const regex_pdb_line(
            // title section
            // 1
            "(HEADER).*|"
            // 2
            "(OBSLTE).*|"
            // 3
            "(TITLE) .*|"
            // 4
            "(SPLIT) .*|"
            // 5       6           7         8
            "(CAVEAT)  ([0-9 ]{2}) (....)    (.{60}).*|"
            // 9
            "(COMPND).*|"
            // 10
            "(SOURCE).*|"
            // 11
            "(KEYWDS).*|"
            // 12
            "(EXPDTA).*|"
            // 13        14
            "(NUMMDL)    ([0-9 ]{4}).*|"
            // 15
            "(MDLTYP).*|"
            // 16
            "(AUTHOR).*|"
            // 17
            "(REVDAT).*|"
            // 18
            "(SPRSDE).*|"
            // 19
            "(JRNL)  .*|"
            // 20     21          22
            "(REMARK) ([0-9 ]{3}) (.{68}).*|"
            // primary structure section
            // 23
            "(DBREF)..*|"
            // 24
            "(SEQADV).*|"
            // 25      26          27  28           29(1) 30(2) 31(3) 32(4) 33(5) 34(6) 35(7) 36(8) 37(9) 38(10)39(11)
            "(SEQRES)  ([0-9 ]{2}) (.) ([0-9 ]{4})  (...) (...) (...) (...) (...) (...) (...) (...) (...) (...) (...) "
            // 40(12)41(13)
            "(...) (...).*|"
            // 42     43     44    45  46         47  48     49
            "(MODRES) (....) (...) (.) ([0-9 ]{4})(.) (...)  (.{41}).*|"
            // heterogen section
            // 50
            "(HET)   .*|"
            // 51
            "(HETNAM).*|"
            // 52
            "(HETSYN).*|"
            // 53
            "(FORMUL).*|"
            // secondary structure section
            // 54     55          56    57    58  59         60  61    62  63         64 65         66      67
            "(HELIX)  ([0-9 ]{3}) (...) (...) (.) ([0-9 ]{4})(.) (...) (.) ([0-9 ]{4})(.)([0-9 ]{2})(.{30}) ([0-9 ]{5})"
            ".*|"
            // 68     69          70   71          72    73 74         75  76    77 78         79 80
            "(SHEET)  ([0-9 ]{3}) (...)([0-9 ]{2}) (...) (.)([0-9 ]{4})(.) (...) (.)([0-9 ]{4})(.)([0-9 ]{2}).*|"
            // connectivity annotation section
            // 81
            "(SSBOND).*|"
            // 82
            "(LINK)  .*|"
            // 83
            "(CISPEP).*|"
            // miscellaneous features section
            // 84
            "(SITE)  .*|"
            // crystallographic and coordinate transformation section
            // 85
            "(CRYST1).*|"
            // 86
            "(ORIGX)..*|"
            // 87
            "(SCALE)..*|"
            // 88
            "(MTRIX)..*|"
            // coordinate section
            // 89        90
            "(MODEL)     ([0-9 ]{4}).*|"
            // 91    92          93    94 95    96 97    98    99                          100
            "(ATOM)  ([0-9 ]{5}) (....)(.)(...) (.)(....)(.)   ([0-9 -]{3}[0-9][.][0-9]{3})([0-9 -]{3}[0-9][.][0-9]{3})"
            // 101                       102                         103                                   104 105
            "([0-9 -]{3}[0-9][.][0-9]{3})([0-9 -]{2}[0-9][.][0-9]{2})([0-9 -]{2}[0-9][.][0-9]{2})          (..)(..).*|"
            // 106
            "(ANISOU).*|"
            // 107   108              109   110 111  112
            "(TER)   ([0-9 ]{5})      (...) (.)(....)(.).*|"
            // 113   114         115   116 117  118 119  120   121                         122
            "(HETATM)([0-9 ]{5}) (....)(.)(...) (.)(....)(.)   ([0-9 -]{3}[0-9][.][0-9]{3})([0-9 -]{3}[0-9][.][0-9]{3})"
            // 123                       124                         125                                   126 127
            "([0-9 -]{3}[0-9][.][0-9]{3})([0-9 -]{2}[0-9][.][0-9]{2})([0-9 -]{2}[0-9][.][0-9]{2})          (..)(..)..*|"
            // 128
            "(ENDMDL).*|"
            // connectivity section
            // 129
            "(CONECT).*|"
            // bookkeeping section
            // 130
            "(MASTER).*|"
            // 131
            "(END).*");

        std::list<std::string> file_lines(tools::file::read_to_string_list(__filename));
        DEBUG << "Read pdb file content:\n" << boost::algorithm::join(file_lines, "\n");

        // to temporarily store data
        caveat_data caveat;
        seqres_data seqres;
        ssdef_data ssdef;
        model_data model;

        DEBUG << "Parsing pdb file";
        for(auto &line : file_lines) { // process each line from the file
          // make all lines have length >=80 by filling up with spaces if length <80, to ensure the regex will match
          line.append(std::max(80 - line.length(), (size_t)0), ' ');

          boost::smatch match;
          // regex_match returns true if all chars of line were matched
          if(!boost::regex_match(line, match, regex_pdb_line)) {
            DEBUG << "No match: " << line;
            continue;
          } // if

          // match[0], same as match.str(), contains the whole string, see
          // http://www.boost.org/doc/libs/1_55_0/libs/regex/example/snippets/regex_iterator_example.cpp
          std::string complete_match(match[0]);
          boost::trim(complete_match); // remove newline at the end (the regex above matches an optional newline)
          DEBUG << "Match: " << complete_match;

          if(match[5].matched) {
            DEBUG << "Found CAVEAT line: " << match[5];
            caveat.process_caveat();
          } // if
          else if(match[13].matched) {
            DEBUG << "Found NUMMDL line: " << match[13];
            model.process_number_models();
          } // if
          else if(match[20].matched) {
            DEBUG << "Found REMARK line: " << match[20] << " " << match[21];
          } // else if
          else if(match[25].matched) {
            DEBUG << "Found SEQRES line: " << match[25];
            std::vector<std::string> residues;
            for(size_t pos(29); pos <= 41; ++pos) {
              residues.push_back(std::string(match[pos]));
            } // for
            seqres.process(match[26], match[27], match[28], residues);
          } // else if
          else if(match[42].matched) {
            DEBUG << "Found MODRES line: " << match[42];
          } // else if
          else if(match[54].matched) {
            DEBUG << "Found HELIX line: " << match[54];
            ssdef.process_helix();
          } // else if
          else if(match[68].matched) {
            DEBUG << "Found SHEET line: " << match[68];
            ssdef.process_strand();
          } // else if
          else if(match[89].matched) {
            DEBUG << "Found MODEL line: " << match[89];
            model.process_model_begin();
          } // else if
          else if(match[91].matched) {
            DEBUG << "Found ATOM line: " << match[91];
            model.process_atom();
          } // else if
          else if(match[107].matched) {
            DEBUG << "Found TER line: " << match[107];
            model.process_terminate();
          } // else if
          else if(match[113].matched) {
            DEBUG << "Found HETATM line: " << match[113];
            model.process_hetatm();
          } // else if
          else if(match[128].matched) {
            DEBUG << "Found ENDMDL line: " << match[128];
            model.process_model_end();
          } // else if
        } // for

        assembly a; // final result

        size_t sequence_no(1); // start with 1
        for(auto c : seqres.get_chain_ids()) {
          a.set(std::string(c, 1),
                molecule(__filename + "/" + std::to_string(sequence_no), ">lcl|sequence", seqres.get_sequence(c)),
                ss(sequence<cchb_dssp>()));
          ++sequence_no; // increase sequence_no, b/c it's not done in the for loop header
        }

        return a;
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
