#include "che/io/file_fasta.h"
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include "tools/file.h"
#include "tools/log.h"
#include "tools/string.h"

namespace biosim {
  namespace che {
    namespace io {
      // (static) reads sequences from a given file
      assembly file_fasta::read(std::string const &__filename) {
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

        size_t sequence_no(1);
        // collect all sequences to try to match them after all of them are read
        std::list<std::pair<std::string, ss>> sss; // pairs of (id, ss); no map, so order is preserved
        std::list<molecule> molecules; // molecules containing id, ps

        // loop over all identifier/sequence pairs
        boost::sregex_iterator itr(file_content.begin(), file_content.end(), regex_multifasta), itr_end;
        for(; itr != itr_end; ++itr, ++sequence_no) { // also increment sequence_no
          boost::smatch match(*itr);

          // remove newlines at the end
          std::string match_identifier(boost::trim_copy(match.str(1))), match_sequence(boost::trim_copy(match.str(2)));

          DEBUG << "Match:\n" << boost::trim_copy(match.str());
          DEBUG << "Identifier:\n" << match_identifier;
          DEBUG << "Sequence:\n" << match_sequence;

          // try to detect which symbol type this sequence has, and based on that fill the sequence into the molecule
          tools::char_function char_isspace = (int (*)(int)) & std::isspace;
          std::set<char> charset_cc(tools::get_unique_char_set_ignore(cc::get_identifier_char_string(), char_isspace));
          std::set<char> charset_cchb_dssp(
              tools::get_unique_char_set_ignore(cchb_dssp::get_identifier_char_string(), char_isspace));
          std::set<char> charset_match(tools::get_unique_char_set_ignore(match_sequence, char_isspace));
          bool is_cc_sequence(std::includes(charset_cc.begin(), charset_cc.end(), charset_match.begin(),
                                            charset_match.end())); // true if charset_match is subset of charset_cc
          bool is_cchb_dssp_sequence(
              std::includes(charset_cchb_dssp.begin(), charset_cchb_dssp.end(), charset_match.begin(),
                            charset_match.end())); // true if charset_match is subset of charset_cchb_dssp

          std::string const file_id(__filename + "/" + std::to_string(sequence_no)); // used for all sequence types

          if(is_cchb_dssp_sequence) { // test for cchb first, b/c it's a subset of cc
            LOG << "Detected sequence type for " << __filename << "/" << sequence_no << ": dssp sequence.";

            // if gap chars are found, but no unknown chars, replace all gap chars with unknown chars
            cchb_dssp unknown(cchb_dssp::specificity_type::unknown), gap(cchb_dssp::specificity_type::gap);
            if(match_sequence.find(gap.get_identifier()) != std::string::npos &&
               match_sequence.find(unknown.get_identifier()) == std::string::npos) {
              std::replace(match_sequence.begin(), match_sequence.end(), gap.get_identifier(),
                           unknown.get_identifier());
              LOG << "Found '" << gap.get_identifier() << "' but no '" << unknown.get_identifier() << "'; replacing '"
                  << gap.get_identifier() << "' with '" << unknown.get_identifier() << "'.";
            } // if

            sss.emplace_back(file_id, convert_to_ss(match_sequence)); // store for now
          } else if(is_cc_sequence) {
            LOG << "Detected sequence type for " << __filename << "/" << sequence_no << ": cc sequence.";

            // store for now; if file had no/empty identifier, use default
            molecules.emplace_back(molecule(file_id, match_identifier.empty() ? ">lcl|sequence" : match_identifier,
                                            convert_to_ps(match_sequence)));
          } // if
          else { // do not throw exception here, the user of this class has to check how many assemblies were found
            std::stringstream ss(match_identifier);
            std::string first_line;
            std::getline(ss, first_line); // read until linebreak from stream into string
            LOG << "Sequence type not identified, ignoring sequence: " << __filename << "/" << sequence_no << ": "
                << first_line;
          } // else
        } // for

        assembly new_assembly;

        // iterate over sss and molecules to match them based on their length
        DEBUG << "Trying to match " << molecules.size() << " cc sequences with " << sss.size() << " dssp sequences.";
        for(auto &current_molecule : molecules) {
          // try to find a matching secondary structure sequence
          bool inserted(false);
          for(auto itr(sss.begin()), itr_end(sss.end()); itr != itr_end; ++itr) {
            if(itr->second.get_sequence().size() == current_molecule.get_length()) {
              LOG << "Assuming " << shorten_file_storage(itr->first) << " is the dssp assignment for "
                  << shorten_file_storage(current_molecule.get_storage()) << ".";
              new_assembly.add(molecule(current_molecule.get_storage(), current_molecule.get_identifier(),
                                        current_molecule.get_ps(), itr->second));
              inserted = true;
              sss.erase(itr);
              break;
            } // if
          } // for

          if(!inserted) { // if not found, just insert the cc sequence without secondary structure
            DEBUG << "No matching dssp assignment found for " << current_molecule.get_storage();
            new_assembly.add(current_molecule);
          } // if
        } // for

        DEBUG << "Inserting " << sss.size() << " unmatched dssp sequences.";
        for(auto p : sss) {
          ps new_ps;
          new_ps.resize(p.second.get_sequence().size(),
                        cc(cc::specificity_type::unknown, cc::monomer_type::l_peptide_linking));
          new_assembly.add(molecule(p.first, ">lcl|sequence", new_ps, p.second));
        } // for

        return new_assembly;
      } // read()

      // (static) convert a string of id_chars to a ps
      ps file_fasta::convert_to_ps(std::string __s) {
        ps new_ps;
        for(auto const &id_char : __s) {
          if(!std::isspace(id_char)) { // ignore whitespace
            new_ps.emplace_back(id_char);
          } // if
        } // for
        return new_ps;
      } // convert_to_ps()
      // (static) convert a string of id_chars to an ss
      ss file_fasta::convert_to_ss(std::string __s) {
        sequence<cchb_dssp> ss_sequence;
        for(auto const &id_char : __s) {
          if(!std::isspace(id_char)) { // ignore whitespace
            ss_sequence.emplace_back(id_char);
          } // if
        } // for
        return ss(ss_sequence);
      } // convert_to_ss()

      // (static) shorten the file storage location to the last two parts
      std::string file_fasta::shorten_file_storage(std::string __s) {
        std::vector<std::string> parts;
        boost::algorithm::split(parts, __s, boost::is_any_of("/"));
        return parts.size() < 2 ? __s : parts[parts.size() - 2] + "/" + parts[parts.size() - 1];
      } // shorten_file_storage()
    } // namespace io
  } // namespace che
} // namespace biosim
