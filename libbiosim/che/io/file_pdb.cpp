#include "che/io/file_pdb.h"
#include "tools/file.h"
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include "che/algo/aligner_dp.h"

namespace biosim {
  namespace che {
    namespace io {
      // stores the caveat data
      class caveat_data {
      public:
        // processes the input strings
        void process_caveat() {}
      }; // class caveat_data

      // stores the modres data
      class modres_data {
      public:
        // use only residue_name as key; chain_id, residue_sequence_number, residue_insertion_code are not always known
        using modres_map = std::map<std::string, std::string>; // map residue_name -> standard_residue_name

      private:
        modres_map _map; // residue name mapping
        mutable std::set<std::string> _used_set; // set of residue names for which the modres mapping was used

      public:
        // processes the input strings
        void process(std::vector<std::string> __line_values) {
          std::string pdb_id(__line_values[0]);
          std::string residue_name(boost::trim_copy(__line_values[1])), chain_id_string(__line_values[2]);
          std::string residue_sequence_number_string(__line_values[3]), residue_insertion_code_string(__line_values[4]);
          std::string standard_residue_name(__line_values[5]);
          std::string description(__line_values[6]);

          if(_map.find(residue_name) != _map.end() && _map[residue_name] != standard_residue_name) {
            LOG << "Found two mappings to standard residue names for " << residue_name << "; using " << residue_name
                << "->" << _map[residue_name] << ", ignoring " << residue_name << "->" << standard_residue_name;
            return;
          } // if

          _map[residue_name] = standard_residue_name;
        } // process()

        // get the map
        modres_map const &get_map() const { return _map; }
        // get the used set of mappings
        std::set<std::string> const &get_used_set() const { return _used_set; }
        // get the cc for the given id, taking modres mapping into account
        cc get_cc(std::string const &__id) const {
          cc this_cc('X'); // default unknown
          try {
            this_cc = cc(__id); // try to find the cc with __id
          } // try
          catch(cc_data_not_found const &__e) { // if no cc with __id is stored, then try modres mapping
            try {
              this_cc = cc(_map.at(__id));
              _used_set.insert(__id);
            } // try
            catch(std::out_of_range const &__e) { // if nothing is found in the map, use the unknown, and tell the user
              LOG << "Residue " << __id
                  << " is not known, and no MODRES information is provided by the pdb file; using unknown instead.";
            } // catch
            catch(cc_data_not_found const &__e) { // if standard_residue_name is not found
              LOG << "Residue " << __id << " is not known, and MODRES provided standard residue name " << _map.at(__id)
                  << " is not known; using unknown instead.";
            } // catch
          } // catch
          return this_cc;
        } // get_cc()
        // get the standard_residue_name for the given residue_name __id
        std::string get_standard_name(std::string const &__id) const { return get_cc(__id).get_identifier(); }
      }; // class modres_data

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
        void process(std::vector<std::string> __line_values, modres_data const &__modres) {
          // trim strings and convert chain id string to char
          std::string const trimmed_serial_number_string(boost::trim_copy(__line_values[0]));
          std::string const trimmed_residue_count_string(boost::trim_copy(std::string(__line_values[2])));
          std::vector<std::string> residues(__line_values.begin() + 3, __line_values.end());
          char const chain_id(boost::lexical_cast<char>(__line_values[1]));
          auto chain_itr(_chains.find(chain_id)); // if equal to _chains.end(), this chain_id is not in the map

          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: seqres_serial_number='" << trimmed_serial_number_string << "' seqres_chain_id='"
                << chain_id << "' seqres_residue_count='" << trimmed_residue_count_string << "' residue_names='"
                << boost::algorithm::join(residues, "|") << "'";

          size_t serial_number;
          try {
            serial_number = boost::lexical_cast<size_t>(trimmed_serial_number_string);
            // check previously stored value
            if(chain_itr != _chains.end() && serial_number != chain_itr->second._last_serial_number + 1) {
              size_t const expected_serial_number(chain_itr->second._last_serial_number + 1);
              DEBUG << "Found seqres_serial_number='" << serial_number << "' differing from expected "
                    << "seqres_serial_number='" << expected_serial_number << "', using expected seqres_serial_number.";
              serial_number = expected_serial_number;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            // default to 1, if we have not seen this chain_id previously
            serial_number = chain_itr == _chains.end() ? 1 : chain_itr->second._last_serial_number + 1;
            DEBUG << "Could not convert seqres_serial_number='" << trimmed_serial_number_string
                  << "' to number, using expected seqres_serial_number='" << serial_number << "'";
          } // catch

          size_t residue_count;
          try {
            residue_count = boost::lexical_cast<size_t>(trimmed_residue_count_string);
            // check previously stored value
            if(chain_itr != _chains.end() && residue_count != chain_itr->second._last_residue_count) {
              size_t const expected_residue_count(chain_itr->second._last_residue_count);
              DEBUG << "Found seqres_residue_count='" << residue_count << "' differing from expected "
                    << "seqres_residue_count='" << expected_residue_count << "', using expected seqres_residue_count.";
              residue_count = expected_residue_count;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            residue_count = chain_itr == _chains.end() ? 0 : chain_itr->second._last_residue_count;
            DEBUG << "Could not convert seqres_residue_count='" << trimmed_residue_count_string << "' to number, using"
                  << "expected residue_count='" << residue_count << "'";
          } // catch

          // save everything
          _chains[chain_id]._last_serial_number = serial_number;
          _chains[chain_id]._last_residue_count = residue_count;
          // loop through all 13 residue name groups and check if there's a residue or space, and save the residues
          for(auto r : residues) {
            std::string const residue_name(boost::trim_copy(r));
            if(residue_name.empty()) {
              continue; // ignore empty names, likely at the end of the SEQRES definition
            } // if

            _chains[chain_id]._cc_sequence.emplace_back(__modres.get_cc(residue_name));
          } // for
        } // process()

        // insert cc [begin..end) at the fron   t of the cc_sequence
        template <class InputIterator>
        void insert_front(char const &__chain_id, InputIterator __begin, InputIterator __end) {
          _chains[__chain_id]._cc_sequence.insert(_chains[__chain_id]._cc_sequence.begin(), __begin, __end);
        } // insert_front()

        // insert cc [begin..end) at the back of the cc_sequence
        template <class InputIterator>
        void insert_back(char const &__chain_id, InputIterator __begin, InputIterator __end) {
          _chains[__chain_id]._cc_sequence.insert(_chains[__chain_id]._cc_sequence.end(), __begin, __end);
        } // insert_back()

        // get chain ids
        std::list<char> get_chain_ids() const {
          std::list<char> chain_ids;
          for(auto p : _chains) {
            chain_ids.emplace_back(p.first);
          } // for
          return chain_ids;
        } // get_chain_ids()

        // returns the ps
        ps const &get_sequence(char const &__chain_id) const { return _chains.at(__chain_id)._cc_sequence; }
      }; // class seqres_data

      // stores the ss definition data
      class ssdef_data {
      private:
        // stores the ss definition data for a single chain
        struct ssdef_chain_data {
          size_t _last_helix_serial_number; // last seen helix serial number
          std::string _last_strand_id; // last seen strand id
          size_t _last_strand_number_in_sheet; // last seen strand number in this sheet (out of number strands in sheet)
          size_t _last_number_strands_in_sheet; // last seen number of strands in this sheet
          ps _cc_sequence; // converted sequence data
          std::set<cchb_dssp_interval> _pool; // converted ss definitions
        }; // struct ssdef_chain_data

        std::map<char, ssdef_chain_data> _chains; // map chain_id->chain_data

      public:
        // processes the input strings
        void process_helix(std::vector<std::string> __line_values, modres_data const &__modres) {
          // check if contents are valid for conversion;
          // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
          // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
          // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
          std::string const serial_number_string(boost::trim_copy(__line_values[0]));
          std::string const id(boost::trim_copy(__line_values[1]));
          std::string const initial_residue_name(boost::trim_copy(__line_values[2]));
          char const initial_residue_chain_id(boost::lexical_cast<char>(__line_values[3]));
          size_t const initial_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__line_values[4])));
          char const initial_residue_insertion_code(boost::lexical_cast<char>(__line_values[5]));
          std::string const terminal_residue_name(boost::trim_copy(__line_values[6]));
          char const terminal_residue_chain_id(boost::lexical_cast<char>(__line_values[7]));
          size_t const terminal_residue_sequence_number(
              boost::lexical_cast<size_t>(boost::trim_copy(__line_values[8])));
          char const terminal_residue_insertion_code(boost::lexical_cast<char>(__line_values[9]));
          std::string const type_string(boost::trim_copy(__line_values[10]));
          std::string const comment(__line_values[11]);
          std::string const length_string(boost::trim_copy(__line_values[12]));

          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: helix_serial_number='" << serial_number_string << "' helix_id='" << id
                << "' initial_residue='" << initial_residue_name << "|" << initial_residue_chain_id << "|"
                << initial_residue_sequence_number << "|" << initial_residue_insertion_code << "' terminal_residue='"
                << terminal_residue_name << "|" << terminal_residue_chain_id << "|" << terminal_residue_sequence_number
                << "|" << terminal_residue_insertion_code << "' helix_type='" << type_string << "' comment='"
                << boost::trim_copy(comment) << "' helix_length='" << length_string << "'";

          // both chain ids should be identical
          if(initial_residue_chain_id != terminal_residue_chain_id) {
            DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                  << "' differing from terminal_residue_chain_id='" << terminal_residue_chain_id << "'";
            return;
          } // if
          char chain_id(initial_residue_chain_id);
          auto chain_itr(_chains.find(chain_id)); // if equal to _chains.end(), this chain_id is not in the map

          // convert strings to size_t and catch exceptions for ones which have a default or way the calculating
          size_t serial_number;
          try {
            serial_number = boost::lexical_cast<size_t>(serial_number_string);
            // check previously stored value
            if(chain_itr != _chains.end() && serial_number != chain_itr->second._last_helix_serial_number + 1) {
              size_t const expected_serial_number(chain_itr->second._last_helix_serial_number + 1);
              DEBUG << "Found helix_serial_number='" << serial_number << "' differing from expected "
                    << "helix_serial_number='" << expected_serial_number << "', using expected helix_serial_number.";
              serial_number = expected_serial_number;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            serial_number = chain_itr == _chains.end() ? 1 : chain_itr->second._last_helix_serial_number + 1;
            DEBUG << "Could not convert helix_serial_number='" << serial_number_string
                  << "' to number, using calculated helix_serial_number='" << serial_number << "'";
          } // catch

          size_t length;
          try {
            length = boost::lexical_cast<size_t>(length_string);
            // check consistency with other values
            if(length != terminal_residue_sequence_number - initial_residue_sequence_number + 1) {
              size_t const expected_length(terminal_residue_sequence_number - initial_residue_sequence_number + 1);
              DEBUG << "Found helix_length='" << length << "' differing expected helix_length='" << length
                    << "', using expected helix_length.";
              length = expected_length;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            // + 1, b/c the initial_residue_sequence_number is part of the helix
            length = terminal_residue_sequence_number - initial_residue_sequence_number + 1;
            DEBUG << "Could not convert helix_length='" << length_string << "' to number, using calculated "
                  << "helix_length='" << length << "'";
          } // catch

          // we could convert type_string into size_t, but to check we'd need to calculate hydrogen bonds;
          // no way to check id and {initial,terminal}_residue_insertion_code

          // save everything
          _chains[chain_id]._last_helix_serial_number = serial_number;
          // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
          _chains[chain_id]._cc_sequence.resize(
              std::max(_chains[chain_id]._cc_sequence.size(), terminal_residue_sequence_number), cc('X'));
          _chains[chain_id]._cc_sequence[initial_residue_sequence_number - 1] = __modres.get_cc(initial_residue_name);
          _chains[chain_id]._cc_sequence[terminal_residue_sequence_number - 1] = __modres.get_cc(terminal_residue_name);
          // create and insert sequence_interval into the pool
          _chains[chain_id]._pool.insert(cchb_dssp_interval(initial_residue_sequence_number - 1,
                                                            terminal_residue_sequence_number - 1, cchb_dssp('H')));
        } // process_helix()
        void process_strand(std::vector<std::string> __line_values, modres_data const &__modres) {
          // check if contents are valid for conversion;
          // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
          // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
          // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
          std::string const strand_number_in_sheet_string(boost::trim_copy(__line_values[0]));
          std::string const id(boost::trim_copy(__line_values[1]));
          std::string const number_strands_in_sheet_string(boost::trim_copy(__line_values[2]));
          std::string const initial_residue_name(boost::trim_copy(__line_values[3]));
          char const initial_residue_chain_id(boost::lexical_cast<char>(__line_values[4]));
          size_t const initial_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__line_values[5])));
          char const initial_residue_insertion_code(boost::lexical_cast<char>(__line_values[6]));
          std::string const terminal_residue_name(boost::trim_copy(__line_values[7]));
          char const terminal_residue_chain_id(boost::lexical_cast<char>(__line_values[8]));
          size_t const terminal_residue_sequence_number(
              boost::lexical_cast<size_t>(boost::trim_copy(__line_values[9])));
          char const terminal_residue_insertion_code(boost::lexical_cast<char>(__line_values[10]));
          std::string const strand_sense_string(boost::trim_copy(__line_values[11]));

          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: strand_number_in_sheet='" << strand_number_in_sheet_string << "' sheet_id='" << id
                << "' number_strands_in_sheet='" << number_strands_in_sheet_string << "' initial_residue='"
                << initial_residue_name << "|" << initial_residue_chain_id << "|" << initial_residue_sequence_number
                << "|" << initial_residue_insertion_code << "' terminal_residue='" << terminal_residue_name << "|"
                << terminal_residue_chain_id << "|" << terminal_residue_sequence_number << "|"
                << terminal_residue_insertion_code << "' strand_sense='" << strand_sense_string << "'";

          // both chain ids should be identical
          if(initial_residue_chain_id != terminal_residue_chain_id) {
            DEBUG << "Ignoring line, initial_residue_chain_id='" << initial_residue_chain_id
                  << "' differing from terminal_residue_chain_id='" << terminal_residue_chain_id << "'";
            return;
          } // if
          char chain_id(initial_residue_chain_id);
          auto chain_itr(_chains.find(chain_id)); // if equal to _chains.end(), this chain_id is not in the map

          size_t strand_number_in_sheet;
          try {
            strand_number_in_sheet = boost::lexical_cast<size_t>(strand_number_in_sheet_string);
            // if we handled a ssdef with this chain_id before, compare data from this line to previously stored values
            if(chain_itr != _chains.end() && id == chain_itr->second._last_strand_id &&
               strand_number_in_sheet != chain_itr->second._last_strand_number_in_sheet + 1) {
              size_t const expected_strand_number_in_sheet(chain_itr->second._last_strand_number_in_sheet + 1);
              DEBUG << "Found strand_number_in_sheet='" << strand_number_in_sheet << "' differing from expected "
                    << "strand_number_in_sheet='" << expected_strand_number_in_sheet
                    << "', using expected strand_number_in_sheet.";
              strand_number_in_sheet = expected_strand_number_in_sheet;
            } // if (no else, because we can't do anything in this case, we don't have old data or we cannot use it)
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            strand_number_in_sheet = chain_itr == _chains.end() || id != chain_itr->second._last_strand_id
                                         ? 1
                                         : chain_itr->second._last_strand_number_in_sheet + 1;
            DEBUG << "Could not convert strand_number_in_sheet='" << strand_number_in_sheet_string
                  << "' to number, using calculated strand_number_in_sheet='" << strand_number_in_sheet << "'";
          } // catch

          size_t number_strands_in_sheet;
          try {
            number_strands_in_sheet = boost::lexical_cast<size_t>(number_strands_in_sheet_string);
            // if we handled a ssdef with this chain_id before, compare data from this line to previously stored values
            if(chain_itr != _chains.end() && id == chain_itr->second._last_strand_id &&
               number_strands_in_sheet != chain_itr->second._last_number_strands_in_sheet) {
              size_t const expected_number_strands_in_sheet(chain_itr->second._last_number_strands_in_sheet);
              DEBUG << "Found number_strands_in_sheet='" << number_strands_in_sheet << "' differing from expected "
                    << "number_strands_in_sheet='" << expected_number_strands_in_sheet
                    << "', using expected number_strands_in_sheet.";
              number_strands_in_sheet = expected_number_strands_in_sheet;
            } // if (no else, because we can't do anything in this case, we don't have old data or we cannot use it)
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            number_strands_in_sheet = chain_itr == _chains.end() || id != chain_itr->second._last_strand_id
                                          ? 1
                                          : chain_itr->second._last_number_strands_in_sheet;
            DEBUG << "Could not convert number_strands_in_sheet='" << number_strands_in_sheet_string
                  << "' to number, using expected number_strands_in_sheet='" << number_strands_in_sheet << "'";
          } // catch

          if(strand_number_in_sheet > number_strands_in_sheet) {
            DEBUG << "Found strand_number_in_sheet='" << strand_number_in_sheet << "' larger than "
                  << "number_strands_in_sheet='" << number_strands_in_sheet << "', increasing number_strands_in_sheet";
            number_strands_in_sheet = strand_number_in_sheet;
          } // if

          // we could try to convert strand_sense_string into size_t, but no checks could be done realistically b/c
          // we would need to evaluate atom positions; no checks are possible on initial_residue_insertion_code and
          // terminal_residue_insertion_code either

          // save everything
          _chains[chain_id]._last_strand_id = id;
          _chains[chain_id]._last_strand_number_in_sheet = strand_number_in_sheet;
          _chains[chain_id]._last_number_strands_in_sheet = number_strands_in_sheet;
          // resize cc_sequence, fill with unknown cc, and set the correct cc for begin and end of the sse
          _chains[chain_id]._cc_sequence.resize(
              std::max(_chains[chain_id]._cc_sequence.size(), terminal_residue_sequence_number), cc('X'));
          _chains[chain_id]._cc_sequence[initial_residue_sequence_number - 1] = __modres.get_cc(initial_residue_name);
          _chains[chain_id]._cc_sequence[terminal_residue_sequence_number - 1] = __modres.get_cc(terminal_residue_name);
          // create and insert sequence_interval into the pool
          _chains[chain_id]._pool.insert(cchb_dssp_interval(initial_residue_sequence_number - 1,
                                                            terminal_residue_sequence_number - 1, cchb_dssp('E')));
        } // process_strand()

        // shift the given chain by __insert_length
        void shift(char const &__chain_id, int __insert_length) {
          // shift sequence
          if(__insert_length >= 0) { // shift sequence to the right by inserting at the front
            _chains[__chain_id]._cc_sequence.insert(_chains[__chain_id]._cc_sequence.begin(), __insert_length, cc('X'));
          } // if
          else { // shift sequence to the left by removing from the front
            ps new_seq;
            new_seq.insert(new_seq.end(), _chains[__chain_id]._cc_sequence.begin() + std::abs(__insert_length),
                           _chains[__chain_id]._cc_sequence.end());
            _chains[__chain_id]._cc_sequence = new_seq;
          } // else
          // shift pool
          std::set<cchb_dssp_interval> new_pool;
          for(cchb_dssp_interval const &sse : _chains[__chain_id]._pool) {
            new_pool.insert(
                cchb_dssp_interval(sse.get_min() + __insert_length, sse.get_max() + __insert_length, sse.get_type()));
          } // for
          _chains[__chain_id]._pool = new_pool;
        } // shift()

        // get chain ids
        std::list<char> get_chain_ids() {
          std::list<char> chain_ids;
          for(auto p : _chains) {
            chain_ids.emplace_back(p.first);
          } // for
          return chain_ids;
        } // get_chain_ids()

        // returns the ps
        ps const &get_sequence(char const &__chain_id) const { return _chains.at(__chain_id)._cc_sequence; }
        // return the pool
        std::set<cchb_dssp_interval> const &get_pool(char const &__chain_id) const {
          return _chains.at(__chain_id)._pool;
        } // get_pool()
      }; // class ssdef_data

      // stores the model data
      class model_data {
      private:
        // stores the model data for a single chain
        struct model_chain_data {
          // default ctor
          model_chain_data()
              : _cc_sequence(), _last_residue_sequence_number(0), _terminated(false), _hetatm_after_ter(0) {}

          ps _cc_sequence; // atom/hetatm converted into a sequence
          size_t _last_residue_sequence_number; // last see residue sequence number
          bool _terminated; // if the chain was terminated with a TER line
          size_t _hetatm_after_ter; // how many atom/hetatm lines were processed after the TER line
        }; //
        std::vector<std::map<char, model_chain_data>> _model_ensemble; // ensemble of chains
        size_t _model_count; // number of models, either default or as specified by NUMMDL
        size_t _current_model_number; // number of the current model
        bool _current_model_complete; // if the current model was completed with a ENDMDL line
        size_t _last_serial_number; // last seen serial number

      public:
        // default ctor
        model_data()
            : _model_ensemble(),
              _model_count(0),
              _current_model_number(0),
              _current_model_complete(false),
              _last_serial_number(0) {}

        // processes the input strings
        void process_number_models(std::vector<std::string> __line_values) {
          std::string const model_count_string(boost::trim_copy(__line_values[0]));
          DEBUG << "Submatches: nummdl='" << model_count_string << "'"; // print the submatches

          try {
            _model_count = boost::lexical_cast<size_t>(model_count_string);
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            _model_count = 1;
            DEBUG << "Could not convert number_models='" << model_count_string
                  << "' to number, using default number_models=" << _model_count;
          } // catch
        } //
        // processes the input strings
        void process_model_begin(std::vector<std::string> __line_values) {
          std::string const model_number_string(boost::trim_copy(__line_values[0]));
          DEBUG << "Submatches: model='" << model_number_string << "'"; // print the submatches

          size_t model_number;
          try {
            model_number = boost::lexical_cast<size_t>(model_number_string);
            if(model_number != _current_model_number + 1) {
              size_t const expected_serial_number(_current_model_number + 1);
              DEBUG << "Found model_number='" << model_number << "' differing from expected model_number='"
                    << expected_serial_number << "', using expected model_number.";
              model_number = expected_serial_number;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            model_number = _current_model_number + 1;
            LOG << "Could not convert model_number='" << model_number_string
                << "' to number, using calculated model_number=" << model_number;
          } // catch

          _model_ensemble.emplace_back(std::map<char, model_chain_data>()); // add new model, initialize everything else
          _current_model_number = model_number;
          _current_model_complete = false;
          _last_serial_number = 0;
        } // process_model_begin()
        // processes the input strings
        void process_model_end() {
          DEBUG << "Processing ENDMDL.";
          _current_model_complete = true;
        } // process_model_end()
        // processes the input strings
        void process_atom(std::vector<std::string> __line_values, modres_data const &__modres) {
          // check if contents are valid for conversion;
          // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
          // exception 2: always try to convert residue sequence number, there's no other way of knowing them
          // (if converting the residue sequence number does not work, lexical_cast will throw an exception)
          std::string const line_type(__line_values[0]);
          std::string const serial_number_string(boost::trim_copy(__line_values[1]));
          std::string const atom_name(boost::trim_copy(__line_values[2]));
          std::string const alternate_location(boost::trim_copy(__line_values[3]));
          std::string const residue_name(boost::trim_copy(__line_values[4]));
          std::string chain_id_string(__line_values[5]);
          char const chain_id(boost::lexical_cast<char>(chain_id_string));
          std::string residue_sequence_number_string(__line_values[6]);
          size_t const residue_sequence_number(
              boost::lexical_cast<size_t>(boost::trim_copy(residue_sequence_number_string)));
          std::string residue_insertion_code_string(__line_values[7]);
          char const residue_insertion_code(boost::lexical_cast<char>(residue_insertion_code_string));
          std::string const x_string(boost::trim_copy(__line_values[8]));
          std::string const y_string(boost::trim_copy(__line_values[9]));
          std::string const z_string(boost::trim_copy(__line_values[10]));
          std::string const occupancy_string(boost::trim_copy(__line_values[11]));
          std::string const temperature_factor_string(boost::trim_copy(__line_values[12]));
          std::string const element(boost::trim_copy(__line_values[13]));
          std::string const charge_string(boost::trim_copy(__line_values[14]));

          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: " << line_type << "='" << serial_number_string << "|" << atom_name << "|"
                << alternate_location << "' residue='" << residue_name << "|" << chain_id << "|"
                << residue_sequence_number << "|" << residue_insertion_code << "' pos='" << x_string << "|" << y_string
                << "|" << z_string << "' occupancy='" << occupancy_string << "' temp_factor='"
                << temperature_factor_string << "' element='" << element << "' charge='" << charge_string << "'";

          // if this pdb contains a single model, there won't be any MODEL or ENDMDL lines; just create a model
          if(_model_ensemble.empty()) {
            _model_ensemble.emplace_back(std::map<char, model_chain_data>());
            _model_count = 1;
            _current_model_number = 1;
          } // if

          if(_current_model_complete) {
            LOG << "Found " << line_type << " line outside of model, after ENDMDL and before a new MODEL; ignoring.";
            return;
          } // if

          // convert strings to size_t and catch exceptions for ones which have a default or way the calculating
          size_t serial_number;
          try {
            serial_number = boost::lexical_cast<size_t>(serial_number_string);
            // check previously stored value
            if(serial_number != _last_serial_number + 1) {
              size_t const expected_serial_number(_last_serial_number + 1);
              DEBUG << "Found serial_number='" << serial_number << "' differing from expected serial_number='"
                    << expected_serial_number << "', using expected serial_number.";
              serial_number = expected_serial_number;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            serial_number = _model_ensemble.empty() ? 1 : _last_serial_number + 1;
            DEBUG << "Could not convert serial_number='" << serial_number_string
                  << "' to number, using calculated serial_number='" << serial_number << "'";
          } // catch

          if(_model_ensemble.back()[chain_id]._terminated) { // check for HETATM after chain TER, store them separately
            ++_model_ensemble.back()[chain_id]._hetatm_after_ter; // count how many HETATM were found
          } else {
            // add missing residues before the current position
            if(_model_ensemble.back()[chain_id]._last_residue_sequence_number + 1 < residue_sequence_number) {
              DEBUG << "Residues with residue_sequence_number="
                    << (_model_ensemble.back()[chain_id]._last_residue_sequence_number + 1) << ".."
                    << (residue_sequence_number - 1) << " missing, filling with unknown residues.";
              _model_ensemble.back()[chain_id]._cc_sequence.resize(residue_sequence_number - 1, cc('X'));
            } // if
            // add current residue; no else branch if it already exists, multiple atoms share the same residue
            if(_model_ensemble.back()[chain_id]._last_residue_sequence_number < residue_sequence_number) {
              _model_ensemble.back()[chain_id]._cc_sequence.emplace_back(__modres.get_cc(residue_name));
            } // if

            // convert position and add atom to cc
            double x, y, z;
            try {
              x = boost::lexical_cast<double>(x_string);
              y = boost::lexical_cast<double>(y_string);
              z = boost::lexical_cast<double>(z_string);
              _model_ensemble.back()[chain_id]._cc_sequence.back().add_determined_atom(
                  atom(atom_name, math::point({x, y, z})));
            } // try
            catch(boost::bad_lexical_cast const &__e) {
              LOG << "Could not convert spatial position '" << x_string << "|" << y_string << "|" << z_string
                  << "' to numbers, ignoring atom='" << serial_number << "|" << atom_name << "' of residue='"
                  << residue_name << "|" << chain_id << "|" << residue_sequence_number << "|" << residue_insertion_code
                  << "'.";
            } // catch
          } // else

          // save everything
          _model_ensemble.back()[chain_id]._last_residue_sequence_number = residue_sequence_number;
          _last_serial_number = serial_number;
        } // process_atom()
        // processes the input strings
        void process_terminate(std::vector<std::string> __line_values, modres_data const &__modres) {
          // check if contents are valid for conversion;
          // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
          // exception 2: always try to convert residue sequence number, there's no other way of knowing them
          // (if converting the residue sequence number does not work, lexical_cast will throw an exception)
          std::string const serial_number_string(boost::trim_copy(__line_values[0]));
          std::string const residue_name(boost::trim_copy(__line_values[1]));
          char const chain_id(boost::lexical_cast<char>(__line_values[2]));
          size_t const residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__line_values[3])));
          char const residue_insertion_code(boost::lexical_cast<char>(__line_values[4]));

          // print the submatches as strings before consistency checks and converting to size_t
          DEBUG << "Submatches: ter='" << serial_number_string << "' residue='" << residue_name << "|" << chain_id
                << "|" << residue_sequence_number << "|" << residue_insertion_code << "'";

          // convert strings to size_t and catch exceptions for ones which have a default or way the calculating
          size_t serial_number;
          try {
            serial_number = boost::lexical_cast<size_t>(serial_number_string);
            // check previously stored value
            if(serial_number != _last_serial_number + 1) {
              size_t const expected_serial_number(_last_serial_number + 1);
              DEBUG << "Found serial_number='" << serial_number << "' differing from expected serial_number='"
                    << expected_serial_number << "', using expected serial_number.";
              serial_number = expected_serial_number;
            } // if
          } // try
          catch(boost::bad_lexical_cast const &__e) {
            serial_number = _model_ensemble.empty() ? 1 : _last_serial_number + 1;
            DEBUG << "Could not convert serial_number='" << serial_number_string
                  << "' to number, using calculated serial_number='" << serial_number << "'";
          } // catch

          if(_model_ensemble.empty()) {
            DEBUG << "Found TER line without model, ignoring.";
            return;
          } // if

          // simply naming
          size_t const last_residue_sequence_number(_model_ensemble.back()[chain_id]._last_residue_sequence_number);
          std::string const last_residue_name(_model_ensemble.back()[chain_id]._cc_sequence.back().get_identifier());
          std::string const standard_residue_name(__modres.get_standard_name(residue_name));
          // verify last cc
          bool correct_residue_sequence_number(last_residue_sequence_number == residue_sequence_number);
          bool correct_residue_name(last_residue_name == residue_name);
          bool correct_standard_residue_name(last_residue_name == standard_residue_name);
          if(!correct_residue_sequence_number || !(correct_residue_name || correct_standard_residue_name)) {
            DEBUG << "Residue on TER line differs from last residue: correct_residue_sequence_number="
                  << correct_residue_sequence_number << " (" << last_residue_sequence_number << " and "
                  << residue_sequence_number << "), correct_residue_name=" << correct_residue_name << " ("
                  << last_residue_name << " and " << residue_name
                  << "), correct_standard_residue_name=" << correct_standard_residue_name << " (" << last_residue_name
                  << " and " << standard_residue_name << "), ignoring.";
            return;
          } // if

          // save everything
          _model_ensemble.back()[chain_id]._terminated = true;
          _last_serial_number = serial_number;
        } // process_terminate()
        // get model data for given chain id
        std::list<model_chain_data> get_models(char const &__chain_id) const {
          std::list<model_chain_data> models;
          for(auto chain_map : _model_ensemble) {
            models.push_back(chain_map.at(__chain_id));
          } // for

          return models;
        } // get_models()
      }; // class model_data

      // adjust seqres and ssdef data based on the given alignment
      void adjust(seqres_data &__seqres, ssdef_data &__ssdef, char const &__chain_id,
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
        // for serial number independent lines, maps line_keyword -> lines, each line consisting of a vector of values
        std::map<std::string, std::list<std::vector<std::string>>> non_model_lines;
        // for serial number dependent lines, vector of values, including the line_keyword
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
            non_model_lines[match[5]].push_back({match[6], match[7], match[8]});
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
            non_model_lines[match[25]].push_back({match[26], match[27], match[28], match[29], match[30], match[31],
                                                  match[32], match[33], match[34], match[35], match[36], match[37],
                                                  match[38], match[39], match[40], match[41]});
          } // else if
          else if(match[42].matched) {
            DEBUG << "Found MODRES line: " << match[42];
            non_model_lines[match[42]].push_back(
                {match[43], match[44], match[45], match[46], match[47], match[48], match[49]});
          } // else if
          else if(match[54].matched) {
            DEBUG << "Found HELIX line: " << match[54];
            non_model_lines[match[54]].push_back({match[55], match[56], match[57], match[58], match[59], match[60],
                                                  match[61], match[62], match[63], match[64], match[65], match[66],
                                                  match[67]});
          } // else if
          else if(match[68].matched) {
            DEBUG << "Found SHEET line: " << match[68];
            non_model_lines[match[68]].push_back({match[69], match[70], match[71], match[72], match[73], match[74],
                                                  match[75], match[76], match[77], match[78], match[79], match[80]});
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
            //             model.process_hetatm();
            model_lines.push_back({match[113], match[114], match[115], match[116], match[117], match[118], match[119],
                                   match[120], match[121], match[122], match[123], match[124], match[125], match[126],
                                   match[127]});
          } // else if
          else if(match[128].matched) {
            DEBUG << "Found ENDMDL line: " << match[128];
            model_lines.push_back({match[128]});
          } // else if
        } // for
        DEBUG << "Parsing found " << non_model_lines["CAVEAT"].size() << " CAVEAT lines, "
              << non_model_lines["SEQRES"].size() << " SEQRES lines, " << non_model_lines["MODRES"].size()
              << " MODRES lines, " << non_model_lines["HELIX"].size() << " HELIX lines, "
              << non_model_lines["SHEET"].size() << " SHEET lines, " << model_lines.size()
              << " model lines (NUMMDL, MODEL, ATOM, TER, HETATM, ENDMDL)";

        // to temporarily store processed data
        caveat_data caveat;
        modres_data modres;
        seqres_data seqres;
        ssdef_data ssdef;
        model_data model;

        DEBUG << "Processing pdb file";
        // process non-model lines
        std::for_each(non_model_lines["CAVEAT"].begin(), non_model_lines["CAVEAT"].end(),
                      [&](std::vector<std::string> const &__line) { caveat.process_caveat(); });
        // process MODRES first, even though it is sorted after SEQRES; needed to evaluate modified residues names
        std::for_each(non_model_lines["MODRES"].begin(), non_model_lines["MODRES"].end(),
                      [&](std::vector<std::string> const &__line) { modres.process(__line); });
        if(!modres.get_used_set().empty()) { // if any mappings were used, tell the user about it
          std::stringstream s;
          for(auto id : modres.get_used_set()) {
            s << id << "->" << modres.get_map().at(id) << "; ";
          } // for
          LOG << "The following " << modres.get_used_set().size()
              << " residue name to standard residue name mappings were using during processing: " << s.str();
        } // id
        std::for_each(non_model_lines["SEQRES"].begin(), non_model_lines["SEQRES"].end(),
                      [&](std::vector<std::string> const &__line) { seqres.process(__line, modres); });
        std::for_each(non_model_lines["HELIX"].begin(), non_model_lines["HELIX"].end(),
                      [&](std::vector<std::string> const &__line) { ssdef.process_helix(__line, modres); });
        std::for_each(non_model_lines["SHEET"].begin(), non_model_lines["SHEET"].end(),
                      [&](std::vector<std::string> const &__line) { ssdef.process_strand(__line, modres); });
        // process model lines
        std::for_each(model_lines.begin(), model_lines.end(), [&](std::vector<std::string> const &__line) {
          std::vector<std::string> line_without_keyword(__line.begin() + 1, __line.end());
          if(__line[0] == "NUMMDL") {
            model.process_number_models(line_without_keyword);
          } else if(__line[0] == "MODEL") {
            model.process_model_begin(line_without_keyword);
          } else if(__line[0] == "ATOM" || __line[0] == "HETATM") {
            model.process_atom(__line, modres);
          } else if(__line[0] == "TER") {
            model.process_terminate(line_without_keyword, modres);
          } else if(__line[0] == "ENDMDL") {
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
