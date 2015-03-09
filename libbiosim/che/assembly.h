#ifndef che_assembly_h
#define che_assembly_h

#include "che/molecule.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // assembly of multiple molecules with additional information
    class assembly {
    public:
      // default ctor
      assembly();
      // ctor from molecule and ss
      assembly(molecule __m, ss __ss);
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // if the molecule with the given chain_id exists
      bool has_molecule(std::string const &__chain_id) const;
      // if the ss with the given chain_id exists
      bool has_ss(std::string const &__chain_id) const;
      // get molecule with given chain_id, throws out of range if no ts with this chain_id exists
      molecule const &get_molecule(std::string const &__chain_id) const;
      // get ss with given chain_id, throws out of range if no ts with this chain_id exists
      ss const &get_ss(std::string const &__chain_id) const;
      // add molecule to assembly, assign next available chain_id, and returns the chain_id
      std::string add(molecule __m);
      // add corresponding molecule and ss to assembly; assign next available chain_id, and returns the chain_id
      std::string add(molecule __m, ss __ss);
      // set a molecule and ss to have a specific chain_id
      void set(std::string const &__chain_id, molecule __m, ss __ss);

    private:
      std::map<std::string, molecule> _molecules; // maps chain_id -> molecule
      std::map<std::string, ss> _ss; // maps chain_id -> ss
    }; // class assembly

    // output operator for assembly
    inline std::ostream &operator<<(std::ostream &__out, assembly const &__a) {
      __out << "Assembly structure: chains=" << __a.get_chain_id_list().size() << "\n";
      for(auto const &chain_id : __a.get_chain_id_list()) {
        __out << "Chain " << chain_id << "\n";
        if(__a.has_molecule(chain_id)) {
          __out << __a.get_molecule(chain_id);
        } // if
        if(__a.has_ss(chain_id)) {
          __out << __a.get_ss(chain_id);
        } // if
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_assembly_h
