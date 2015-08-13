#include "che/io/file_pssm.h"
#include "tools/file.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // (static) read ps from profile data of a given file
      assembly file_pssm::read(std::string const &__filename) {
        pssm_data data(read_pssm_ascii(__filename));

        std::string const id_chars = "ARNDCQEGHILKMFPSTWYV"; // defines the order, depends on the regex
        ps seq;
        for(pssm_amino_acid_data aa_data : data.aa_data) {
          cc::weight_map weights; // create a weight_map
          for(size_t pos(0); pos < 20; ++pos) { // and fill it
            cc c(id_chars[pos], cc::monomer_type::l_peptide_linking); // get the id from the id_char
            weights[c.get_identifier()] = aa_data.percent[pos];
          }
          seq.emplace_back(cc(weights));
        }
        return assembly(structure(__filename, ">lcl|sequence", seq));
      } // read()

      // (static) reads data from a pssm ascii file into an intermediate data structure
      file_pssm::pssm_data file_pssm::read_pssm_ascii(std::string const &__filename) {
        // regex for a single data line of the pssm ascii format; it will match header and amino acid lines only, no
        // comments, or K or Lambda value lines;
        // the regex was created from pssm ascii files from PSIBLAST 2.2.25+ and 2.2.30+;
        // header line groups explanation:
        // 1=Order of the amino acids;
        // amino acid line groups explanation:
        // 2=Amino acid id; 3=Amino acid single letter code; 4-23=Score for all 20 amino acids; 24-43=Percentage of all
        // amino acids; 44=information content; 45=relative weight of gapless real matches to pseudocounts.
        static boost::regex const regex_pssm_ascii_line(
            // 1
            " *(A {2,3}R {2,3}N {2,3}D {2,3}C {2,3}Q {2,3}E {2,3}G {2,3}H {2,3}I {2,3}L {2,3}K {2,3}M {2,3}F {2,3}P"
            " {2,3}S {2,3}T {2,3}W {2,3}Y {2,3}V   A {2,3}R {2,3}N {2,3}D {2,3}C {2,3}Q {2,3}E {2,3}G {2,3}H {2,3}I"
            " {2,3}L {2,3}K {2,3}M {2,3}F {2,3}P {2,3}S {2,3}T {2,3}W {2,3}Y {2,3}V).*|"
            // 2            3    4                   5                   6                   7
            " *([0-9]{1,5}) (.) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9])"
            // 8                   9                   10                  11                  12
            " +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9])"
            // 13                  14                  15                  16                  17
            " +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9])"
            // 18                  19                  20                  21                  22
            " +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9]) +([0-9-]{0,2}[0-9])"
            // 23                  24            25            26            27            28            29
            " +([0-9-]{0,2}[0-9]) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3})"
            // 30            31            32            33            34            35            36
            " +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3})"
            // 37            38            39            40            41            42            43
            " +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3}) +([0-9]{1,3})"
            // 44                45
            " +([0-9].[0-9]{2}) +((?:[0-9].[0-9]{2})|inf).*");

        std::list<std::string> file_lines(tools::file::read_to_string_list(__filename));
        DEBUG << "Read pssm_ascii file content:\n" << boost::algorithm::join(file_lines, "\n");

        DEBUG << "Parsing pssm_ascii file";
        bool found_header(false);
        pssm_data data;
        size_t last_id(0);
        for(auto &line : file_lines) { // process each line from the file
          boost::smatch match;
          // regex_match returns true if all chars of line were matched
          if(!boost::regex_match(line, match, regex_pssm_ascii_line)) {
            DEBUG << "No match: " << line;
            continue;
          } // if

          // match[0], same as match.str(), contains the whole string, see
          // http://www.boost.org/doc/libs/1_55_0/libs/regex/example/snippets/regex_iterator_example.cpp
          std::string complete_match(match[0]);
          boost::trim(complete_match); // remove newline at the end (the regex above matches an optional newline)
          DEBUG << "Match: " << complete_match;

          if(match[1].matched) {
            DEBUG << "Found table header";
            // dont print submatches, there is only one
            found_header = true;
          } // if
          else if(match[2].matched) {
            DEBUG << "Found table row";
            // trimming strings is not necessary, the regex does not match white space
            DEBUG << "Submatches: id='" << match[2] << "' aa='" << match[3] << "' score='" << match[4] << "|"
                  << match[5] << "|" << match[6] << "|" << match[7] << "|" << match[8] << "|" << match[9] << "|"
                  << match[10] << "|" << match[11] << "|" << match[12] << "|" << match[13] << "|" << match[14] << "|"
                  << match[15] << "|" << match[16] << "|" << match[17] << "|" << match[18] << "|" << match[19] << "|"
                  << match[20] << "|" << match[21] << "|" << match[22] << "|" << match[23] << "' percent='" << match[24]
                  << "|" << match[25] << "|" << match[26] << "|" << match[27] << "|" << match[28] << "|" << match[29]
                  << "|" << match[30] << "|" << match[31] << "|" << match[32] << "|" << match[33] << "|" << match[34]
                  << "|" << match[35] << "|" << match[36] << "|" << match[37] << "|" << match[38] << "|" << match[39]
                  << "|" << match[40] << "|" << match[41] << "|" << match[42] << "|" << match[43] << "' information='"
                  << match[44] << "' pseudocount_weight='" << match[45] << "'";

            if(!found_header) { // need to find the header first to ensure the assumed order of amino acids is correct
              DEBUG << "Skipping table row, table header was not found.";
              continue;
            } // if

            size_t id;
            try {
              id = boost::lexical_cast<size_t>(match[2]);
              // check previously stored value
              if(id != last_id + 1) {
                size_t const expected_id(last_id + 1);
                DEBUG << "Found id='" << id << "' differing from expected id='" << expected_id << "', skipping row.";
                continue;
              } // if
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              id = last_id + 1;
              DEBUG << "Could not convert id='" << match[2] << "' to number, using expected id='" << id << "'";
            } // catch
            last_id = id;

            char const type(boost::lexical_cast<char>(match[3]));

            double information(0.0);
            try {
              information = boost::lexical_cast<double>(match[44]);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              DEBUG << "Could not convert information='" << match[44]
                    << "' to number, using information=" << information;
            } // catch

            double pseudocount_weight(1.0);
            try {
              pseudocount_weight = boost::lexical_cast<double>(match[45]);
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              DEBUG << "Could not convert pseudocount_weight='" << match[45]
                    << "' to number, using pseudocount_weight=" << pseudocount_weight;
            } // catch

            pssm_amino_acid_data aa_data(id, type, information, pseudocount_weight);
            for(size_t pos(0); pos < 20; ++pos) {
              int score(0);
              try {
                score = boost::lexical_cast<int>(match[4 + pos]);
              } // try
              catch(boost::bad_lexical_cast const &__e) {
                DEBUG << "Could not convert score='" << match[4 + pos] << "' to number, using score=" << score;
              } // catch
              aa_data.score.emplace_back(score);

              size_t percent(0);
              try {
                percent = boost::lexical_cast<size_t>(match[24 + pos]);
              } // try
              catch(boost::bad_lexical_cast const &__e) {
                DEBUG << "Could not convert percent='" << match[24 + pos] << "' to number, using percent=" << percent;
              } // catch
              aa_data.percent.emplace_back(percent);
            } // for

            data.aa_data.emplace_back(aa_data);
          } // else if
        } // for

        return data;
      } // read_pssm_ascii()
    } // namespace io
  } // namespace che
} // namespace biosim
