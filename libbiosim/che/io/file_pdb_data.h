#ifndef che_io_file_pdb_modres_data_h
#define che_io_file_pdb_modres_data_h

#include <array>
#include <map>
#include "che/sequence_interval.h"

namespace biosim {
  namespace che {
    namespace io {
      // stores the caveat data of a pdb file
      class pdb_caveat_data {
      public:
        // processes the input strings
        void process();
      }; // class pdb_caveat_data

      // stores the modres data of a pdb file
      class pdb_modres_data {
      public:
        // processes the input strings to fill the map
        void process(std::array<std::string, 7> __values);
        // get the map
        std::map<std::string, std::string> const &get_map() const;
        // get the known cc type for the given residue_name id, taking modres mapping into account
        cc get_cc(std::string const &__id) const;
        // get the known residue_name for the given residue_name id; uses the standard_residue_name if in modres
        std::string get_cc_name(std::string const &__id) const;

      private:
        // mapping residue_name -> standard_residue_name; make it as simple as possible, use only residue_name as key;
        // chain_id, residue_sequence_number, residue_insertion_code are not always known
        std::map<std::string, std::string> _map;
      }; // class pdb_modres_data

      // stores the seqres data of a pdb file
      class pdb_seqres_data {
      public:
        // processes the input strings
        void process(std::array<std::string, 16> __values, pdb_modres_data const &__modres);
        // insert cc [begin..end) at the front of the cc_sequence with the given chain id
        template <class InputIterator>
        void insert_front(char const &__chain_id, InputIterator __begin, InputIterator __end) {
          _chains[__chain_id]._cc_sequence.insert(_chains[__chain_id]._cc_sequence.begin(), __begin, __end);
        } // insert_front()
        // insert cc [begin..end) at the back of the cc_sequence with the given chain id
        template <class InputIterator>
        void insert_back(char const &__chain_id, InputIterator __begin, InputIterator __end) {
          _chains[__chain_id]._cc_sequence.insert(_chains[__chain_id]._cc_sequence.end(), __begin, __end);
        } // insert_back()
        // get chain ids
        std::list<char> get_chain_ids() const;
        // get the ps for the given chain id
        ps const &get_sequence(char const &__chain_id) const;

      private:
        // stores the seqres data for a single chain
        struct seqres_chain_data {
          size_t _last_serial_number; // last seen serial number; starts at 1 for each chain, + 1 for each SEQRES line
          size_t _last_residue_count; // last seen residue count
          ps _cc_sequence; // converted sequence data
        }; // struct seqres_chain_data

        std::map<char, seqres_chain_data> _chains; // map chain_id->chain_data
      }; // class pdb_seqres_data

      // stores the ss definition data of a pdb file
      class pdb_ssdef_data {
      public:
        // processes the input strings
        void process_helix(std::array<std::string, 13> __values, pdb_modres_data const &__modres);
        // processes the input strings
        void process_strand(std::array<std::string, 12> __values, pdb_modres_data const &__modres);
        // shift the given chain by __insert_length
        void shift(char const &__chain_id, int __insert_length);
        // get chain ids
        std::list<char> get_chain_ids() const;
        // returns the ps for the given chain id
        ps const &get_sequence(char const &__chain_id) const;
        // return the pool for the given chain id
        std::set<cchb_dssp_interval> const &get_pool(char const &__chain_id) const;

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
      }; // class pdb_ssdef_data

      // stores the model data of a pdb file
      class pdb_model_data {
      public:
        // default ctor
        pdb_model_data();
        // processes the input strings
        void process_number_models(std::string __model_count_string);
        // processes the input strings
        void process_model_begin(std::string __model_number_string);
        // processes the input strings
        void process_model_end();
        // processes the input strings
        void process_atom(std::array<std::string, 15> __values, pdb_modres_data const &__modres);
        // processes the input strings
        void process_terminate(std::array<std::string, 5> __values, pdb_modres_data const &__modres);
        // get chain ids
        std::list<char> get_chain_ids() const;
        // get sequences for given chain id
        std::list<ps> get_sequences(char const &__chain_id) const;

      private:
        // stores the model data for a single chain
        struct model_chain_data {
          // default ctor
          model_chain_data();

          ps _cc_sequence; // atom/hetatm converted into a sequence
          size_t _last_residue_sequence_number; // last see residue sequence number
          bool _terminated; // if the chain was terminated with a TER line
          size_t _hetatm_after_ter; // how many atom/hetatm lines were processed after the TER line
        }; // struct model_chain_data

        std::vector<std::map<char, model_chain_data>> _model_ensemble; // ensemble of models, each model with chains
        size_t _model_count; // number of models, either default or as specified by NUMMDL
        size_t _current_model_number; // number of the current model
        bool _current_model_complete; // if the current model was completed with a ENDMDL line
        size_t _last_serial_number; // last seen serial number
      }; // class pdb_model_data
    } // namespace io
  } // namespace che
} // namespace biosim

#endif // che_io_file_pdb_modres_data_h
