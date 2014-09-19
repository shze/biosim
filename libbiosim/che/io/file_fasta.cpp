#include "che/io/file_fasta.h"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include "tools/file.h"
#include "tools/log.h"

namespace biosim {
  namespace che {
    namespace io {
      // (static) reads sequences from a given file
      complex file_fasta::read(std::string const &__filename) {
        // regex string for matching a single fasta sequence (identifier/sequence pair) out of a multifasta
        // each iteration matches a single fasta with [0] the complete match, [1] the identifier, [2] the sequence
        static boost::regex const regex_multifasta(
            "(" // group all possible types of identifier within '(..)'
            "(?:;[^\n]*[[:space:]]*)+|" // multiline identifier starting with ';', use [^\n] to match until eol; or
            "(?:>[^\n]*[[:space:]]*)" // single line identifier starting with '>'; add '?:' to not add matches from
            ")?" // these capture groups '(..)' to the result; last ? to make identifier optional
            "([^>;*]+)\\*?[[:space:]]*"); // sequence, allow all characters except identifier begins and sequence ends

        std::string file_content(boost::algorithm::join(tools::file::read_to_string_list(__filename), "\n"));
        DEBUG << "Read fasta file content:\n" << file_content;

        complex all_molecules;

        // loop over all identifier/sequence pairs
        boost::sregex_iterator itr(file_content.begin(), file_content.end(), regex_multifasta), itr_end;
        for(; itr != itr_end; ++itr) {
          boost::smatch match(*itr);

          // remove newlines at the end
          std::string const match_identifier(boost::trim_copy(match.str(1))),
              match_sequence(boost::trim_copy(match.str(2)));

          DEBUG << "Match:\n" << boost::trim_copy(match.str());
          DEBUG << "Identifier:\n" << match_identifier;
          DEBUG << "Sequence:\n" << match_sequence;

          // try to detect which symbol type this sequence has, and based on that fill the sequence into the molecule
          std::set<char> charset_cc(get_unique_chars_ignore_whitespace(cc::get_identifier_char_string()));
          std::set<char> charset_match(get_unique_chars_ignore_whitespace(match_sequence));
          bool is_cc_sequence(std::includes(charset_cc.begin(), charset_cc.end(), charset_match.begin(),
                                            charset_match.end())); // true if charset_match is a subset of charset_ccd

          if(is_cc_sequence) {
            LOG << "Detected file type for " << __filename << ": fasta file, cc sequence.";
            // if file had no/empty identifier, use default
            molecule new_molecule(__filename, match_identifier.empty() ? ">lcl|sequence" : match_identifier);
            new_molecule.set_ps(convert_to_ps(match_sequence));
            all_molecules.add(new_molecule);
          } // if
          else { // do not throw exception here, the user of this class has to check how many molecules were found
            std::stringstream ss(match_identifier);
            std::string first_line;
            std::getline(ss, first_line); // read until linebreak from stream into string
            LOG << "Sequence type not identified, ignoring sequence." << first_line;
          } // else
        } // for

        return all_molecules;
      } // read()

      // (static) return the set of unique chars in the given string without whitespace
      std::set<char> file_fasta::get_unique_chars_ignore_whitespace(std::string const &__s) {
        std::set<char> unique;
        for(auto const &c : __s) {
          if(!std::isspace(c)) { // ignore whitespace
            unique.insert(c);
          } // if
        } // for

        return unique;
      } // get_unique_chars_ignore_whitespace()

      // (static) convert a string of id_chars to a ps of cc's
      ps file_fasta::convert_to_ps(std::string __s) {
        ps new_ps;
        for(auto const &id_char : __s) {
          if(!std::isspace(id_char)) { // ignore whitespace
            new_ps.emplace_back(id_char);
          } // if
        } // for
        return new_ps;
      } // convert_to_ps()
    } // namespace io
  } // namespace che
} // namespace biosim
