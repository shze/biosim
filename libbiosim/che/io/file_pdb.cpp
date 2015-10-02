#include "che/io/file_pdb.h"
#include "tools/file.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include "che/algo/aligner_dp.h"
#include "che/io/file_pdb_data.h"

namespace biosim {
  namespace che {
    namespace io {
      // class to adjust positions of two structures based on an alignment created with the passed aligner
      class adjuster {
      public:
        // update structures to match positions using an alignment created with the aligner; labels are used for output
        static void update(che::structure &__structure1, std::string const &__label1, che::structure &__structure2,
                           std::string const &__label2, che::algo::aligner_dp const &__aligner) {
          // create intro output string
          std::string const intro(__label1 + " and " + __label2 + " differ; ");
          // save original lengths before we change them
          size_t original_length1(__structure1.get_length()), original_length2(__structure2.get_length());

          // create alignment
          che::scored_alignment alignment(
              __aligner.align_multiple({che::alignment(__structure1), che::alignment(__structure2)}).front());

          // begin cases
          size_t const &unaligned_length_begin1(alignment.get_alignment().get_begins()[0]);
          size_t const &unaligned_length_begin2(alignment.get_alignment().get_begins()[1]);
          int begin_diff(unaligned_length_begin1 - unaligned_length_begin2);
          if(begin_diff > 0) { // structure1 starts before structure2
            if(is_undetermined(__structure1, 0, begin_diff)) { // if structure1 begin contains undetermined ps/ss
              LOG << intro << "remove 0.." << (begin_diff - 1) << " from " << __label1;
              remove(__structure1, 0, begin_diff);
            } // if
            else { // if structure1 begin contains at least some determined ps/ss
              LOG << intro << "insert 0.." << (begin_diff - 1) << " from " << __label1 << " at begin of " << __label2;
              insert(__structure2, 0, __structure1, 0, begin_diff);
            } // else
          } // if
          else if(begin_diff < 0) { // structure2 starts before structure1
            if(is_undetermined(__structure2, 0, -begin_diff)) { // if structure2 begin contains undetermined ps/ss
              LOG << intro << "remove 0.." << (-begin_diff - 1) << " from " << __label2;
              remove(__structure2, 0, -begin_diff);
            } // if
            else { // if structure2 begin contains at least some determined ps/ss
              LOG << intro << "insert 0.." << (-begin_diff - 1) << " from " << __label2 << " at begin of " << __label1;
              insert(__structure1, 0, __structure2, 0, -begin_diff);
            } // else
          } // else if

          // end cases
          size_t unaligned_length_end1(original_length1 - alignment.get_alignment().get_ends()[0]);
          size_t unaligned_length_end2(original_length2 - alignment.get_alignment().get_ends()[1]);
          int end_diff(unaligned_length_end1 - unaligned_length_end2);
          if(end_diff > 0) { // structure1 ends after structure2
            size_t begin(__structure1.get_length() - end_diff), end(__structure1.get_length());
            if(is_undetermined(__structure1, begin, end)) { // if structure1 end contains undetermined ps/ss
              LOG << intro << "remove " << begin << ".." << (end - 1) << " from " << __label1;
              remove(__structure1, begin, end);
            } // if
            else { // if structure1 end contains at least some determined ps/ss
              LOG << intro << "insert " << begin << ".." << (end - 1) << " from " << __label1 << " at end of "
                  << __label2;
              insert(__structure2, __structure2.get_length(), __structure1, begin, end);
            } // else
          } // if
          else if(end_diff < 0) { // structure2 ends after structure1
            size_t begin(__structure2.get_length() + end_diff), end(__structure2.get_length());
            if(is_undetermined(__structure2, begin, end)) { // if structure2 end contains undetermined ps/ss
              LOG << intro << "remove " << begin << ".." << (end - 1) << " from " << __label2;
              remove(__structure2, begin, end);
            } // if
            else { // if structure2 end contains determined ps/ss
              LOG << intro << "insert " << begin << ".." << (end - 1) << " from " << __label2 << " at end of "
                  << __label1;
              insert(__structure1, __structure1.get_length(), __structure2, begin, end);
            } // else
          } // else if
        } // update()
        // return if positions __begin..__end are undetermined in __s
        static bool is_undetermined(che::structure &__s, size_t __begin, size_t __end) {
          bool undetermined(true);
          for(; __begin < __end; ++__begin) {
            if(!__s.get_ps().at(__begin).is_unknown() ||
               (__s.get_ss().defined() && !__s.get_ss().get_sequence().at(__begin).is_unknown())) {
              undetermined = false;
              break;
            } // if
          } // for
          return undetermined;
        } // is_undetermined()
        // remove the positions __begin..__end - 1 from __s
        static void remove(che::structure &__s, size_t __begin, size_t __end) {
          // sequence
          ps new_ps;
          new_ps.insert(new_ps.end(), __s.get_ps().begin(), __s.get_ps().begin() + __begin);
          new_ps.insert(new_ps.end(), __s.get_ps().begin() + __end, __s.get_ps().end());
          // pool
          std::set<cchb_dssp_interval> new_sses;
          for(cchb_dssp_interval const &sse : __s.get_ss().get_sses()) {
            size_t new_min(sse.get_min() >= __end ? sse.get_min() - __end + __begin : sse.get_min());
            size_t new_max(sse.get_max() >= __end ? sse.get_max() - __end + __begin : sse.get_max());
            new_sses.emplace(new_min, new_max, sse.get_type());
          } // for
          // save
          ss new_ss(new_sses, new_ps.size());
          __s = __s.get_ss().defined() ? che::structure("", "", new_ps, new_ss) : che::structure(new_ps);
        } // remove
        // insert the positions __from_begin..__from_end - 1 at __to_pos
        static void insert(che::structure &__to, size_t __to_pos, che::structure &__from, size_t __from_begin,
                           size_t __from_end) {
          // sequence
          ps new_ps;
          new_ps.insert(new_ps.end(), __to.get_ps().begin(), __to.get_ps().begin() + __to_pos);
          new_ps.insert(new_ps.end(), __from.get_ps().begin() + __from_begin, __from.get_ps().begin() + __from_end);
          new_ps.insert(new_ps.end(), __to.get_ps().begin() + __to_pos, __to.get_ps().end());
          // pool
          std::set<cchb_dssp_interval> new_sses;
          for(cchb_dssp_interval const &sse : __to.get_ss().get_sses()) {
            size_t new_min(sse.get_min() >= __to_pos ? sse.get_min() + __from_end - __from_begin : sse.get_min());
            size_t new_max(sse.get_max() >= __to_pos ? sse.get_max() + __from_end - __from_begin : sse.get_max());
            new_sses.emplace(new_min, new_max, sse.get_type());
          } // for
          // save
          ss new_ss(new_sses, new_ps.size());
          __to = __to.get_ss().defined() ? che::structure("", "", new_ps, new_ss) : che::structure(new_ps);
        } // insert()
      }; // class adjuster

      // count mismatching known cc and total known cc
      std::pair<size_t, size_t> count_mismatches(ps const &__cc_seq, ps const &__ss_seq) {
        size_t min_length(std::min(__cc_seq.size(), __ss_seq.size()));
        size_t known_cc(0), mismatched_known_cc(0);
        for(size_t pos(0); pos < min_length; ++pos) {
          if(!__cc_seq[pos].is_unknown() && !__ss_seq[pos].is_unknown()) {
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

        // collect all chain ids
        std::list<char> seqres_chain_ids(seqres.get_chain_ids());
        std::list<char> model_chain_ids(model.get_chain_ids());
        std::list<char> ssdef_chain_ids(ssdef.get_chain_ids());
        // combine all chain ids into all_chain_ids; every chain id in this set has at least one kind of data
        std::set<char> all_chain_ids(seqres_chain_ids.begin(), seqres_chain_ids.end());
        all_chain_ids.insert(model_chain_ids.begin(), model_chain_ids.end());
        all_chain_ids.insert(ssdef_chain_ids.begin(), ssdef_chain_ids.end());

        // create before the for loop
        che::score::ev_alignment::cc_cm_function_ptr b62_ptr(new che::score::cm_cc_blosum(true));
        che::algo::aligner_dp aligner(
            che::score::ev_alignment(b62_ptr, -b62_ptr->get_unknown_score(), b62_ptr->get_min_score()));

        size_t sequence_no(0); // start with 0, increase at the beginning of each loop
        for(auto const &c : all_chain_ids) {
          ++sequence_no; // increase sequence_no, b/c it's not done in the for loop header

          // prepare shared data
          std::string storage(__filename + "/" + std::to_string(sequence_no));
          std::string id(">lcl|sequence");
          std::string chain_id(1, c);

          // check which data exist for this chain
          bool has_seqres(std::find(seqres_chain_ids.begin(), seqres_chain_ids.end(), c) != seqres_chain_ids.end());
          bool has_model(std::find(model_chain_ids.begin(), model_chain_ids.end(), c) != model_chain_ids.end());
          bool has_ssdef(std::find(ssdef_chain_ids.begin(), ssdef_chain_ids.end(), c) != ssdef_chain_ids.end());
          LOG << "Processing chain " << c << " (" << (has_seqres ? "" : "no ") << "seqres data, "
              << (has_model ? "" : "no ") << "model data, " << (has_ssdef ? "" : "no ") << "ss data)";

          if(has_seqres + has_model + has_ssdef == 1) { // if only one kind of data exists, save and done
            std::list<che::structure> this_chain_ensemble;
            if(has_seqres) {
              this_chain_ensemble.emplace_back(storage, id, seqres.get_sequence(c));
            } // if
            else if(has_model) {
              for(auto &this_ps : model.get_sequences(c)) {
                this_chain_ensemble.emplace_back(che::structure(storage, id, this_ps));
              } // for
            } // else if
            else /*if(has_ssdef)*/ {
              this_chain_ensemble.emplace_back(storage, id, ssdef.get_sequence(c),
                                               ss(ssdef.get_pool(c), ssdef.get_sequence(c).size()));
            } // else
            a.set(chain_id, che::structure_ensemble(this_chain_ensemble.begin(), this_chain_ensemble.end()));
            continue;
          } // if

          // at this point, at least 2 kinds of data exist; combine all three while checking for missing data
          // combine seqres and ssdef first, b/c each model needs to be aligned b/c it could have different atom data
          // check if the current chain is a peptide chain, b/c ss can only handle peptide secondary structure
          bool is_peptide(seqres.get_sequence(c).front().get_monomer_type() == cc::monomer_type::l_peptide_linking);
          structure seqres_ssdef_structure;
          if(has_seqres && has_ssdef && is_peptide) {
            structure seqres_structure(seqres.get_sequence(c));
            // pass the lenth into ss in case the pool has overlapping intervals that cause the last one to be removed
            structure ssdef_structure("", "", ssdef.get_sequence(c),
                                      ss(ssdef.get_pool(c), ssdef.get_sequence(c).size()));

            // align and adjust structures
            adjuster::update(seqres_structure, "SEQRES", ssdef_structure, "HELIX/SHEET", aligner);
            // check mismatches
            std::pair<size_t, size_t> mismatches(count_mismatches(seqres_structure.get_ps(), ssdef_structure.get_ps()));
            if(mismatches.first > 0) {
              LOG << "SEQRES and HELIX/SHEET mismatch in " << mismatches.first << "/" << mismatches.second
                  << " known positions for chain " << c << ", ignoring structure: " << __filename << "/"
                  << std::to_string(sequence_no);
              continue;
            } // if

            ss this_ss(ssdef_structure.get_ss().get_sses(), seqres_structure.get_length());
            seqres_ssdef_structure = structure(storage, id, seqres_structure.get_ps(), this_ss);
          } // if
          else if(!has_seqres && has_ssdef) {
            seqres_ssdef_structure = structure(storage, id, ssdef.get_sequence(c), che::ss(ssdef.get_pool(c)));
          } // else if
          else {
            seqres_ssdef_structure = structure(storage, id, seqres.get_sequence(c));
          } // else

          // combine model sequence and seqres/ssdef sequence for each model
          if(model.get_sequences(c).empty()) {
            a.set(chain_id, che::structure_ensemble(seqres_ssdef_structure));
            continue;
          } // if
          std::list<che::structure> this_chain_ensemble;
          for(auto const &this_ps : model.get_sequences(c)) {
            structure model_structure(this_ps);

            // align and adjust structures
            adjuster::update(model_structure, "ATOM", seqres_ssdef_structure, "SEQRES/HELIX/SHEET", aligner);
            // check mismatches
            std::pair<size_t, size_t> mismatches(
                count_mismatches(model_structure.get_ps(), seqres_ssdef_structure.get_ps()));
            if(mismatches.first > 0) {
              LOG << "ATOM and SEQRES/HELIX/SHEET mismatch in " << mismatches.first << "/" << mismatches.second
                  << " known positions for chain " << c << ", ignoring structure: " << __filename << "/"
                  << std::to_string(sequence_no);
              continue;
            } // if

            // copy missing cc datan from seqres_ssdef_structure into model_structure (atoms that were not determined)
            ps model_ps(model_structure.get_ps());
            ps seqres_ssdef_ps(seqres_ssdef_structure.get_ps());
            for(size_t pos(0); pos < model_ps.size() && pos < seqres_ssdef_ps.size(); ++pos) {
              if(model_ps[pos].is_unknown()) {
                model_ps[pos] = seqres_ssdef_ps[pos];
              } // if
            } // for

            // save structure in ensemble, either with or without ss
            if(has_ssdef) {
              this_chain_ensemble.emplace_back(storage, id, model_ps, seqres_ssdef_structure.get_ss());
            } // if
            else {
              this_chain_ensemble.emplace_back(storage, id, model_ps);
            } // else
          } // for

          if(!this_chain_ensemble.empty()) { // if ensemble is not empty, save ensemble in assembly
            a.set(chain_id, che::structure_ensemble(this_chain_ensemble.begin(), this_chain_ensemble.end()));
          } // if
        } // for

        return a;
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
