#include "che/io/file_pdb_data.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace biosim {
  namespace che {
    namespace io {
      // processes the input strings
      void pdb_caveat_data::process() {}

      // processes the input strings to fill the map
      void pdb_modres_data::process(std::array<std::string, 7> __values) {
        std::string pdb_id(__values[0]);
        std::string residue_name(boost::trim_copy(__values[1])), chain_id_string(__values[2]);
        std::string residue_sequence_number_string(__values[3]), residue_insertion_code_string(__values[4]);
        std::string standard_residue_name(__values[5]), description(__values[6]);

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << "MODRES submatches: pdb_id='" << pdb_id << "' residue='" << residue_name << "|" << chain_id_string
              << "|" << residue_sequence_number_string << "' standard_residue='" << standard_residue_name
              << "' description='" << description << "'";

        if(_map.find(residue_name) != _map.end() && _map[residue_name] != standard_residue_name) {
          LOG << "Found two mappings to standard residue names for " << residue_name << "; using " << residue_name
              << "->" << _map[residue_name] << ", ignoring " << residue_name << "->" << standard_residue_name;
          return;
        } // if

        _map[residue_name] = standard_residue_name;
      } // process()
      // get the map
      std::map<std::string, std::string> const &pdb_modres_data::get_map() const { return _map; }
      // get the known cc type for the given residue_name id, taking modres mapping into account
      cc pdb_modres_data::get_cc(std::string const &__id) const {
        cc this_cc('X'); // default unknown
        try {
          this_cc = cc(__id); // try to find the cc with __id
        } // try
        catch(cc_data_not_found const &__e) { // if no cc with __id is stored, then try modres mapping
          try {
            this_cc = cc(_map.at(__id));
          } // try
          catch(std::out_of_range const &__e) { // if nothing is found in the map, use the unknown, and tell the user
            LOG << "Residue " << __id << " is undefined, and no MODRES data was found in the pdb file; using unknown.";
          } // catch
          catch(cc_data_not_found const &__e) { // if standard_residue_name is not found
            LOG << "Residue " << __id << " is undefined, and MODRES provided standard residue name " << _map.at(__id)
                << " is undefined; using unknown.";
          } // catch
        } // catch
        return this_cc;
      } // get_cc()
      // get the known residue_name for the given residue_name id; uses the standard_residue_name if in modres
      std::string pdb_modres_data::get_cc_name(std::string const &__id) const { return get_cc(__id).get_identifier(); }

      // processes the input strings
      void pdb_seqres_data::process(std::array<std::string, 16> __values, pdb_modres_data const &__modres) {
        // trim strings and convert chain id string to char
        std::string const trimmed_serial_number_string(boost::trim_copy(__values[0]));
        std::string const trimmed_residue_count_string(boost::trim_copy(std::string(__values[2])));
        std::vector<std::string> residues(__values.begin() + 3, __values.end());
        char const chain_id(boost::lexical_cast<char>(__values[1]));
        auto chain_itr(_chains.find(chain_id)); // if equal to _chains.end(), this chain_id is not in the map

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << "SEQRES submatches: seqres_serial_number='" << trimmed_serial_number_string << "' seqres_chain_id='"
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
          DEBUG << "Could not convert seqres_residue_count='" << trimmed_residue_count_string << "' to number, using "
                << "expected residue_count='" << residue_count << "'";
        } // catch

        // save everything
        _chains[chain_id]._last_serial_number = serial_number;
        _chains[chain_id]._last_residue_count = residue_count;
        // loop through all 13 residue name groups and check if there's a residue or space, and save the residues
        for(auto const &r : residues) {
          std::string const residue_name(boost::trim_copy(r));
          if(residue_name.empty()) {
            continue; // ignore empty names, likely at the end of the SEQRES definition
          } // if

          _chains[chain_id]._cc_sequence.emplace_back(__modres.get_cc(residue_name));
        } // for
      } // process()
      // get chain ids
      std::list<char> pdb_seqres_data::get_chain_ids() const {
        std::list<char> chain_ids;
        for(auto const &p : _chains) {
          chain_ids.emplace_back(p.first);
        } // for
        return chain_ids;
      } // get_chain_ids()
      // get the ps for the given chain id
      ps const &pdb_seqres_data::get_sequence(char const &__chain_id) const {
        return _chains.at(__chain_id)._cc_sequence;
      } // get_sequence()

      // processes the input strings
      void pdb_ssdef_data::process_helix(std::array<std::string, 13> __values, pdb_modres_data const &__modres) {
        // check if contents are valid for conversion;
        // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
        // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
        // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
        std::string const serial_number_string(boost::trim_copy(__values[0]));
        std::string const id(boost::trim_copy(__values[1]));
        std::string const initial_residue_name(boost::trim_copy(__values[2]));
        char const initial_residue_chain_id(boost::lexical_cast<char>(__values[3]));
        size_t const initial_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__values[4])));
        char const initial_residue_insertion_code(boost::lexical_cast<char>(__values[5]));
        std::string const terminal_residue_name(boost::trim_copy(__values[6]));
        char const terminal_residue_chain_id(boost::lexical_cast<char>(__values[7]));
        size_t const terminal_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__values[8])));
        char const terminal_residue_insertion_code(boost::lexical_cast<char>(__values[9]));
        std::string const type_string(boost::trim_copy(__values[10]));
        std::string const comment(__values[11]);
        std::string const length_string(boost::trim_copy(__values[12]));

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << "HELIX submatches: helix_serial_number='" << serial_number_string << "' helix_id='" << id
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
      // processes the input strings
      void pdb_ssdef_data::process_strand(std::array<std::string, 12> __values, pdb_modres_data const &__modres) {
        // check if contents are valid for conversion;
        // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
        // exception 2: always try to convert residue sequence numbers, if there's no other way of knowing them
        // (if converting residue sequence numbers does not work, lexical_cast will throw an exception)
        std::string const strand_number_in_sheet_string(boost::trim_copy(__values[0]));
        std::string const id(boost::trim_copy(__values[1]));
        std::string const number_strands_in_sheet_string(boost::trim_copy(__values[2]));
        std::string const initial_residue_name(boost::trim_copy(__values[3]));
        char const initial_residue_chain_id(boost::lexical_cast<char>(__values[4]));
        size_t const initial_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__values[5])));
        char const initial_residue_insertion_code(boost::lexical_cast<char>(__values[6]));
        std::string const terminal_residue_name(boost::trim_copy(__values[7]));
        char const terminal_residue_chain_id(boost::lexical_cast<char>(__values[8]));
        size_t const terminal_residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__values[9])));
        char const terminal_residue_insertion_code(boost::lexical_cast<char>(__values[10]));
        std::string const strand_sense_string(boost::trim_copy(__values[11]));

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << "SHEET submatches: strand_number_in_sheet='" << strand_number_in_sheet_string << "' sheet_id='" << id
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
      // get chain ids
      std::list<char> pdb_ssdef_data::get_chain_ids() const {
        std::list<char> chain_ids;
        for(auto const &p : _chains) {
          chain_ids.emplace_back(p.first);
        } // for
        return chain_ids;
      } // get_chain_ids()
      // returns the ps for the given chain id
      ps const &pdb_ssdef_data::get_sequence(char const &__chain_id) const {
        return _chains.at(__chain_id)._cc_sequence;
      } // get_sequence()
      // return the pool for the given chain id
      std::set<cchb_dssp_interval> const &pdb_ssdef_data::get_pool(char const &__chain_id) const {
        return _chains.at(__chain_id)._pool;
      } // get_pool()

      // default ctor
      pdb_model_data::model_chain_data::model_chain_data()
          : _cc_sequence(), _last_residue_sequence_number(0), _terminated(false), _hetatm_after_ter(0) {}
      // default ctor
      pdb_model_data::pdb_model_data()
          : _model_ensemble(),
            _model_count(0),
            _current_model_number(0),
            _current_model_complete(false),
            _last_serial_number(0) {}
      // processes the input strings
      void pdb_model_data::process_number_models(std::string __model_count_string) {
        std::string const model_count_string(boost::trim_copy(__model_count_string));
        DEBUG << "NUMMDL submatches: nummdl='" << model_count_string << "'"; // print the submatches

        try {
          _model_count = boost::lexical_cast<size_t>(model_count_string);
        } // try
        catch(boost::bad_lexical_cast const &__e) {
          _model_count = 1;
          DEBUG << "Could not convert number_models='" << model_count_string
                << "' to number, using default number_models=" << _model_count;
        } // catch
      } // process_number_models()
      // processes the input strings
      void pdb_model_data::process_model_begin(std::string __model_number_string) {
        std::string const model_number_string(boost::trim_copy(__model_number_string));
        DEBUG << "MODEL submatches: model='" << model_number_string << "'"; // print the submatches

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
          DEBUG << "Could not convert model_number='" << model_number_string
                << "' to number, using calculated model_number=" << model_number;
        } // catch

        _model_ensemble.emplace_back(std::map<char, model_chain_data>()); // add new model, initialize everything else
        _current_model_number = model_number;
        _current_model_complete = false;
        _last_serial_number = 0;
      } // process_model_begin()
      // processes the input strings
      void pdb_model_data::process_model_end() {
        DEBUG << "ENDMDL found, has no submatches.";
        _current_model_complete = true;
      } // process_model_end()
      // processes the input strings
      void pdb_model_data::process_atom(std::array<std::string, 15> __values, pdb_modres_data const &__modres) {
        // check if contents are valid for conversion;
        // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
        // exception 2: always try to convert residue sequence number, there's no other way of knowing them
        // (if converting the residue sequence number does not work, lexical_cast will throw an exception)
        std::string const line_type(__values[0]);
        std::string const serial_number_string(boost::trim_copy(__values[1]));
        std::string const atom_name(boost::trim_copy(__values[2]));
        std::string const alternate_location(boost::trim_copy(__values[3]));
        std::string const residue_name(boost::trim_copy(__values[4]));
        std::string const chain_id_string(__values[5]);
        char const chain_id(boost::lexical_cast<char>(chain_id_string));
        std::string residue_sequence_number_string(__values[6]);
        size_t const residue_sequence_number(
            boost::lexical_cast<size_t>(boost::trim_copy(residue_sequence_number_string)));
        std::string residue_insertion_code_string(__values[7]);
        char const residue_insertion_code(boost::lexical_cast<char>(residue_insertion_code_string));
        std::string const x_string(boost::trim_copy(__values[8]));
        std::string const y_string(boost::trim_copy(__values[9]));
        std::string const z_string(boost::trim_copy(__values[10]));
        std::string const occupancy_string(boost::trim_copy(__values[11]));
        std::string const temperature_factor_string(boost::trim_copy(__values[12]));
        std::string const element(boost::trim_copy(__values[13]));
        std::string const charge_string(boost::trim_copy(__values[14]));

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << line_type << " submatches: atom='" << serial_number_string << "|" << atom_name << "|"
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
      void pdb_model_data::process_terminate(std::array<std::string, 5> __values, pdb_modres_data const &__modres) {
        // check if contents are valid for conversion;
        // exception 1: lexical_cast on char should be fine, as the submatches are only single char strings;
        // exception 2: always try to convert residue sequence number, there's no other way of knowing them
        // (if converting the residue sequence number does not work, lexical_cast will throw an exception)
        std::string const serial_number_string(boost::trim_copy(__values[0]));
        std::string const residue_name(boost::trim_copy(__values[1]));
        char const chain_id(boost::lexical_cast<char>(__values[2]));
        size_t const residue_sequence_number(boost::lexical_cast<size_t>(boost::trim_copy(__values[3])));
        char const residue_insertion_code(boost::lexical_cast<char>(__values[4]));

        // print the submatches as strings before consistency checks and converting to size_t
        DEBUG << "TER submatches: ter='" << serial_number_string << "' residue='" << residue_name << "|" << chain_id
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
        std::string const standard_residue_name(__modres.get_cc_name(residue_name));
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
      // get chain ids
      std::list<char> pdb_model_data::get_chain_ids() const {
        std::list<char> chain_ids;
        for(auto const &m : _model_ensemble) {
          for(auto const &p : m) {
            chain_ids.emplace_back(p.first);
          } // for
        } // for
        return chain_ids;
      } // get_chain_ids()
      // get sequences for given chain id
      std::list<ps> pdb_model_data::get_sequences(char const &__chain_id) const {
        std::list<ps> sequences;
        for(auto const &chain_map : _model_ensemble) {
          sequences.push_back(chain_map.at(__chain_id)._cc_sequence);
        } // for

        return sequences;
      } // get_sequences()
    } // namespace io
  } // namespace che
} // namespace biosim
