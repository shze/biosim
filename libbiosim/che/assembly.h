#ifndef che_assembly_h
#define che_assembly_h

#include "che/structure_ensemble.h"

namespace biosim {
  namespace che {
    // assembly of multiple ensembles, one ensemble, i.e. a list of structure samples, for each chain id
    class assembly {
    public:
      // default ctor
      assembly();
      // ctor from structure
      explicit assembly(structure __s);
      // ctor from structure_ensemble
      explicit assembly(structure_ensemble __e);
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // get the size of the ensemble for the given chain id
      size_t get_ensemble_size(std::string const &__chain_id) const;
      // returns if an ensemble with the given chain id exists
      bool has_ensemble(std::string const &__chain_id) const;
      // get the ensemble for the given chain id
      structure_ensemble const &get_ensemble(std::string const &__chain_id) const;
      // get the first structure of the structure samples; throws out of range if no structure with this chain_id exists
      structure const &get_first_structure(std::string const &__chain_id) const;
      // get all first structures
      std::vector<structure> get_first_structures() const;
      // add structure to assembly, assign next available chain_id, and returns the chain_id
      std::string add(structure __s);
      // add ensemble to assembly, assign next available chain_id, and returns the chain_id
      std::string add(structure_ensemble __e);
      // set a structure to have a specific chain_id
      void set(std::string __chain_id, structure __s);
      // set an ensemble to have a specific chain_id
      void set(std::string __chain_id, structure_ensemble __e);

    private:
      std::map<std::string, structure_ensemble> _chains; // maps chain_id -> structure_ensemble
    }; // class assembly

    // output operator for assembly
    inline std::ostream &operator<<(std::ostream &__out, assembly const &__a) {
      __out << "Assembly: chains=" << __a.get_chain_id_list().size() << "\n";
      for(auto const &chain_id : __a.get_chain_id_list()) {
        // get ensemble for this chain id
        structure_ensemble const &e(__a.get_ensemble(chain_id));
        // output all occupancy values into the stringstream
        std::stringstream occupancy_stream;
        if(!e.get_samples().empty()) {
          std::for_each(e.get_samples().begin(), --e.get_samples().end(),
                        [&](structure_sample const &__s) { occupancy_stream << __s.occupancy << ", "; });
          occupancy_stream << e.get_samples().rbegin()->sample;
        } // if
        // output
        __out << "Chain " << chain_id << ": samples=" << e.get_samples().size() << ", occupancies={"
              << occupancy_stream.str() << "}\n" << e.get_samples().front().sample;
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_assembly_h
