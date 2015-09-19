#ifndef che_assembly_h
#define che_assembly_h

#include "che/structure.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // assembly of multiple structures
    class assembly {
    public:
      // default ctor
      assembly();
      // ctor from structure
      explicit assembly(structure __s);
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // get all structures
      std::vector<structure> get_structures() const;
      // if the structure with the given chain_id exists
      bool has_structure(std::string const &__chain_id) const;
      // get structure with given chain_id, throws out of range if no ts with this chain_id exists
      structure const &get_structure(std::string const &__chain_id) const;
      // add structure to assembly, assign next available chain_id, and returns the chain_id
      std::string add(structure __s);
      // set a structure to have a specific chain_id
      void set(std::string __chain_id, structure __s);

    private:
      std::map<std::string, structure> _chains; // maps chain_id -> structure
    }; // class assembly

    // output operator for assembly
    inline std::ostream &operator<<(std::ostream &__out, assembly const &__a) {
      __out << "Assembly structure: chains=" << __a.get_chain_id_list().size() << "\n";
      for(auto const &chain_id : __a.get_chain_id_list()) {
        __out << "Chain " << chain_id << "\n" << __a.get_structure(chain_id);
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_assembly_h
