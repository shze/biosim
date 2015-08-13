#include "che/io/file_psipredv.h"
#include "tools/file.h"
#include "tools/log.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // (static) reads ps and ss from a given file
      assembly file_psipredv::read(std::string const &__filename) {
        // regex for a single data line of the psipred vertical format; comment or empty lines do not match this regex;
        // groups will be used to construct sequences.
        static boost::regex const regex_psipredv_line(
            " *([[:digit:]]+) +([[:upper:]]) +([[:upper:]])" // line#, aa one letter code, ss one letter code
            " +([[:digit:]](?:\\.[[:digit:]]+)?)" // this and the two following lines are the probabilities
            " +([[:digit:]](?:\\.[[:digit:]]+)?)" // for the 3 secondary structure types C, H, E
            " +([[:digit:]](?:\\.[[:digit:]]+)?)");

        std::list<std::string> file_lines(tools::file::read_to_string_list(__filename));
        DEBUG << "Read psipred file content:\n" << boost::algorithm::join(file_lines, "\n");

        // create two sequences, one for each type
        ps cc_sequence;
        sequence<cchb_dssp> cchb_sequence;

        for(auto const &line : file_lines) {
          boost::smatch match;
          if(boost::regex_match(line, match, regex_psipredv_line)) {
            // match[0], same as match.str(), contains the whole string, see
            // http://www.boost.org/doc/libs/1_55_0/libs/regex/example/snippets/regex_iterator_example.cpp
            std::string line(match[0]);
            boost::trim(line); // remove newline at the end (the regex above matches an optional newline)
            DEBUG << "Match: " << line;

            // convert everything into the correct type
            size_t const pos(boost::lexical_cast<size_t>(match[1]));
            char const aa_char(boost::lexical_cast<char>(match[2]));
            char const ss_char(boost::lexical_cast<char>(match[3]));
            double const prob_c(boost::lexical_cast<double>(match[4]));
            double const prob_h(boost::lexical_cast<double>(match[5]));
            double const prob_e(boost::lexical_cast<double>(match[6]));

            DEBUG << "Submatches: pos=" << pos << "; aa=" << aa_char << "; ss=" << ss_char << "; prob_C=" << prob_c
                  << "; prob_H=" << prob_h << "; prob_E=" << prob_e;

            cc_sequence.emplace_back(aa_char);
            cchb_sequence.emplace_back(ss_char, cchb_dssp::weight_map({{'C', prob_c}, {'H', prob_h}, {'E', prob_e}}));
          } // if
        } // for

        return cc_sequence.empty() ? assembly()
                                   : assembly(structure(__filename, ">lcl|sequence", cc_sequence, ss(cchb_sequence)));
      } // read()
    } // namespace io
  } // namespace che
} // namespace biosim
