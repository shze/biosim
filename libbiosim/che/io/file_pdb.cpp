#include "che/io/file_pdb.h"
#include "tools/file.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include "che/algo/aligner_dp.h"
#include "che/io/file_pdb_data.h"

namespace biosim {
  namespace che {
    namespace io {
      // adjust seqres and ssdef data based on the given alignment
      void adjust(pdb_seqres_data &__seqres, pdb_ssdef_data &__ssdef, char const &__chain_id,
                  che::scored_alignment const &__alignment) {
        // save original lengths before we change seqres or ssdef
        size_t original_cc_seq_length(__seqres.get_sequence(__chain_id).size());
        size_t original_ss_seq_length(__ssdef.get_sequence(__chain_id).size());

        // check unaligned sequences at the front (there should be exactly two structures in this alignment)
        size_t cc_unaligned_length_begin(__alignment.get_alignment().get_begins()[0]);
        size_t ss_unaligned_length_begin(__alignment.get_alignment().get_begins()[1]);
        int begin_diff(cc_unaligned_length_begin - ss_unaligned_length_begin);
        // if begin_diff == 0, do nothing
        if(begin_diff > 0) { // cc_seq has cc that are not in ss_seq
          LOG << "SEQRES and HELIX/SHEET positions differ for chain " << __chain_id
              << "; increase HELIX/SHEET positions by " << begin_diff << ".";
          __ssdef.shift(__chain_id, begin_diff);
        } // if
        else if(begin_diff < 0) { // ss_seq has cc that are not in cc_seq
          size_t first_sse_begin(__ssdef.get_pool(__chain_id).begin()->get_min());
          // if all 0..abs(begin_diff) cc of ss_seq are unknown (i.e. no sse starts there), remove these cc from ss_seq
          if(std::abs(begin_diff) <= first_sse_begin) {
            LOG << "SEQRES and HELIX/SHEET positions differ for chain " << __chain_id
                << "; decrease HELIX/SHEET positions by " << std::abs(begin_diff) << ".";
            __ssdef.shift(__chain_id, begin_diff);
          } // if
          else { // an sse starts there, copy 0..abs(begin_diff) cc from ss_seq to cc_seq
            // it would be possible to only copy from first_sse_begin..abs(begin_diff), but then additionally all sse
            // would have to shift back by begin_diff
            LOG << "HELIX/SHEET starts before SEQRES for chain " << __chain_id << "; increase SEQRES positions by "
                << std::abs(begin_diff) << " and insert unknown amino acids at the beginning.";
            __seqres.insert_front(__chain_id, __ssdef.get_sequence(__chain_id).begin(),
                                  __ssdef.get_sequence(__chain_id).begin() - begin_diff);
          } // else
        } // else if

        // check unaligned sequence in the back
        size_t cc_unaligned_length_end(original_cc_seq_length - __alignment.get_alignment().get_ends()[0]);
        size_t ss_unaligned_length_end(original_ss_seq_length - __alignment.get_alignment().get_ends()[1]);
        int end_diff(cc_unaligned_length_end - ss_unaligned_length_end);
        // if end_diff == 0, do nothing
        // if end_diff > 0, cc_seq has cc that are not in ss_seq: add unknown at back of ss_seq by setting total length
        if(end_diff < 0) { // ss_seq has cc that are not in cc_seq: add unknown at back of cc_seq
          // no check needed of end of ss_seq is just unknown b/c ss_seq ends with the last sse
          LOG << "HELIX/SHEET ends after SEQRES for chain " << __chain_id << "; insert " << std::abs(end_diff)
              << " unknown amino acids at the end of SEQRES.";
          __seqres.insert_back(__chain_id, __ssdef.get_sequence(__chain_id).end() - std::abs(end_diff),
                               __ssdef.get_sequence(__chain_id).end());
        } // if
      } // adjust()

      // count mismatching known cc and total known cc
      std::pair<size_t, size_t> count_mismatches(ps const &__cc_seq, ps const &__ss_seq) {
        size_t min_length(std::min(__cc_seq.size(), __ss_seq.size()));
        size_t known_cc(0), mismatched_known_cc(0);
        for(size_t pos(0); pos < min_length; ++pos) {
          if(!__ss_seq[pos].is_unknown()) {
            ++known_cc;
            if(__cc_seq[pos].get_identifier_char() != __ss_seq[pos].get_identifier_char()) {
              ++mismatched_known_cc;
            } // if
          } // if
        } // for
        return std::pair<size_t, size_t>(mismatched_known_cc, known_cc);
      } // count_mismatches()

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
            "(SHEET)  ([0-9 ]{3}) (...)([0-9 ]{2}) (...) (.)([0-9 ]{4})(.) (...) (.)([0-9 ]{4})(.)([0-9 -]{2}).*|"
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
            "([0-9 -]{3}[0-9][.][0-9]{3})([0-9 -]{2}[0-9][.][0-9]{2})([0-9 -]{2}[0-9][.][0-9]{2})          (..)(..).*|"
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

        // to temporarily store parsed data
        using caveat_value_array = std::array<std::string, 3>;
        using modres_value_array = std::array<std::string, 7>;
        using seqres_value_array = std::array<std::string, 16>;
        using helix_value_array = std::array<std::string, 13>;
        using sheet_value_array = std::array<std::string, 12>;
        using atom_value_array = std::array<std::string, 15>;
        using ter_value_array = std::array<std::string, 5>;
        // for non-model lines, save every value except the keyword
        std::list<caveat_value_array> caveat_lines;
        std::list<modres_value_array> modres_lines;
        std::list<seqres_value_array> seqres_lines;
        std::list<helix_value_array> helix_lines;
        std::list<sheet_value_array> sheet_lines;
        // for model lines (serial number dependent lines), the vector of values includes the keyword
        std::list<std::vector<std::string>> model_lines;

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
            caveat_lines.push_back({match[6], match[7], match[8]});
          } // if
          else if(match[13].matched) {
            DEBUG << "Found NUMMDL line: " << match[13];
            model_lines.push_back({match[13], match[14]});
          } // if
          else if(match[20].matched) {
            DEBUG << "Found REMARK line: " << match[20] << " " << match[21];
          } // else if
          else if(match[25].matched) {
            DEBUG << "Found SEQRES line: " << match[25];
            seqres_lines.push_back({match[26], match[27], match[28], match[29], match[30], match[31], match[32],
                                    match[33], match[34], match[35], match[36], match[37], match[38], match[39],
                                    match[40], match[41]});
          } // else if
          else if(match[42].matched) {
            DEBUG << "Found MODRES line: " << match[42];
            modres_lines.push_back({match[43], match[44], match[45], match[46], match[47], match[48], match[49]});
          } // else if
          else if(match[54].matched) {
            DEBUG << "Found HELIX line: " << match[54];
            helix_lines.push_back({match[55], match[56], match[57], match[58], match[59], match[60], match[61],
                                   match[62], match[63], match[64], match[65], match[66], match[67]});
          } // else if
          else if(match[68].matched) {
            DEBUG << "Found SHEET line: " << match[68];
            sheet_lines.push_back({match[69], match[70], match[71], match[72], match[73], match[74], match[75],
                                   match[76], match[77], match[78], match[79], match[80]});
          } // else if
          else if(match[89].matched) {
            DEBUG << "Found MODEL line: " << match[89];
            model_lines.push_back({match[89], match[90]});
          } // else if
          else if(match[91].matched) {
            DEBUG << "Found ATOM line: " << match[91];
            model_lines.push_back({match[91], match[92], match[93], match[94], match[95], match[96], match[97],
                                   match[98], match[99], match[100], match[101], match[102], match[103], match[104],
                                   match[105]});
          } // else if
          else if(match[107].matched) {
            DEBUG << "Found TER line: " << match[107];
            model_lines.push_back({match[107], match[108], match[109], match[110], match[111], match[112]});
          } // else if
          else if(match[113].matched) {
            DEBUG << "Found HETATM line: " << match[113];
            model_lines.push_back({match[113], match[114], match[115], match[116], match[117], match[118], match[119],
                                   match[120], match[121], match[122], match[123], match[124], match[125], match[126],
                                   match[127]});
          } // else if
          else if(match[128].matched) {
            DEBUG << "Found ENDMDL line: " << match[128];
            model_lines.push_back({match[128]});
          } // else if
        } // for
        DEBUG << "Parsing found " << caveat_lines.size() << " CAVEAT lines, " << seqres_lines.size()
              << " SEQRES lines, " << modres_lines.size() << " MODRES lines, " << helix_lines.size() << " HELIX lines, "
              << sheet_lines.size() << " SHEET lines, " << model_lines.size()
              << " model lines (NUMMDL, MODEL, ATOM, TER, HETATM, ENDMDL)";

        // to temporarily store processed data
        pdb_caveat_data caveat;
        pdb_modres_data modres;
        pdb_seqres_data seqres;
        pdb_ssdef_data ssdef;
        pdb_model_data model;

        DEBUG << "Processing pdb file";
        // process non-model lines
        std::for_each(caveat_lines.begin(), caveat_lines.end(),
                      [&](caveat_value_array const &__a) { caveat.process(); });
        std::for_each(modres_lines.begin(), modres_lines.end(),
                      [&](modres_value_array const &__a) { modres.process(__a); });
        std::for_each(seqres_lines.begin(), seqres_lines.end(),
                      [&](seqres_value_array const &__a) { seqres.process(__a, modres); });
        std::for_each(helix_lines.begin(), helix_lines.end(),
                      [&](helix_value_array const &__a) { ssdef.process_helix(__a, modres); });
        std::for_each(sheet_lines.begin(), sheet_lines.end(),
                      [&](sheet_value_array const &__a) { ssdef.process_strand(__a, modres); });
        // process model lines
        std::for_each(model_lines.begin(), model_lines.end(), [&](std::vector<std::string> const &__v) {
          if(__v[0] == "NUMMDL") {
            model.process_number_models(__v[1]);
          } else if(__v[0] == "MODEL") {
            model.process_model_begin(__v[1]);
          } else if(__v[0] == "ATOM" || __v[0] == "HETATM") {
            model.process_atom(atom_value_array{__v[0], __v[1], __v[2], __v[3], __v[4], __v[5], __v[6], __v[7], __v[8],
                                                __v[9], __v[10], __v[11], __v[12], __v[13], __v[14]},
                               modres);
          } else if(__v[0] == "TER") {
            model.process_terminate(ter_value_array{__v[1], __v[2], __v[3], __v[4], __v[5]}, modres);
          } else if(__v[0] == "ENDMDL") {
            model.process_model_end();
          }
        });

        assembly a; // final result

        che::score::ev_alignment::cc_cm_function_ptr b62_ptr(new che::score::cm_cc_blosum(true));
        che::algo::aligner_dp aligner(
            che::score::ev_alignment(b62_ptr, -b62_ptr->get_unknown_score(), b62_ptr->get_min_score()));
        size_t sequence_no(1); // start with 1
        for(auto c : seqres.get_chain_ids()) {
          ss this_ss;

          bool is_peptide_chain(seqres.get_sequence(c).front().get_monomer_type() ==
                                cc::monomer_type::l_peptide_linking);
          if(is_peptide_chain) {
            this_ss = ss(std::set<cchb_dssp_interval>(), seqres.get_sequence(c).size()); // coil ss with correct length

            std::list<char> ssdef_chain_ids(ssdef.get_chain_ids());
            bool chain_has_ss(std::find(ssdef_chain_ids.begin(), ssdef_chain_ids.end(), c) != ssdef_chain_ids.end());
            if(chain_has_ss) {
              structure cc_seq("", "", seqres.get_sequence(c));
              // pass the lenth into ss in case the pool has overlapping intervals that cause the last one to be removed
              structure ss_seq("", "", ssdef.get_sequence(c), ss(ssdef.get_pool(c), ssdef.get_sequence(c).size()));
              che::scored_alignment best_alignment(
                  aligner.align_multiple({che::alignment(cc_seq), che::alignment(ss_seq)}).front());
              adjust(seqres, ssdef, c, best_alignment);

              std::pair<size_t, size_t> mismatches(count_mismatches(seqres.get_sequence(c), ssdef.get_sequence(c)));
              if(mismatches.first > 0) {
                LOG << "SEQRES and HELIX/SHEET mismatch in " << mismatches.first << "/" << mismatches.second
                    << " known positions for chain " << c << ", ignoring structure: " << __filename << "/"
                    << std::to_string(sequence_no);
                continue;
              } // if
              this_ss = ss(ssdef.get_pool(c), seqres.get_sequence(c).size()); // set ss from pdb file
            } // if
          } // if

          std::string storage(__filename + "/" + std::to_string(sequence_no));
          structure this_structure = is_peptide_chain
                                         ? structure(storage, ">lcl|sequence", seqres.get_sequence(c), this_ss)
                                         : structure(storage, ">lcl|sequence", seqres.get_sequence(c));
          a.set(std::string(1, c), this_structure);
          ++sequence_no; // increase sequence_no, b/c it's not done in the for loop header
        }

        return a;
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
