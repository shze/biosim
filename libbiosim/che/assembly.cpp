#include "che/assembly.h"
#include "tools/mapper.h"

namespace biosim {
  namespace che {
    // default ctor
    assembly::assembly() : _chains() {}
    // ctor from structure
    assembly::assembly(structure __s) : _chains() { add(__s); }
    // get list of all chain_ids
    std::list<std::string> assembly::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto c : _chains) {
        chain_ids.push_back(c.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // get all structures
    std::vector<structure> assembly::get_structures() const {
      std::vector<structure> structures;
      structures.reserve(_chains.size());
      for(auto const &c : _chains) {
        structures.push_back(c.second);
      } // for
      return structures;
    } // get_structures()
    // if the structure with the given chain_id exists
    bool assembly::has_structure(std::string const &__chain_id) const {
      return _chains.find(__chain_id) != _chains.end();
    } // has_structure()
    // get structure with given chain_id, throws out of range if no ts with this chain_id exists
    structure const &assembly::get_structure(std::string const &__chain_id) const { return _chains.at(__chain_id); }
    // add structure to assembly, assign next available chain_id, and returns the chain_id
    std::string assembly::add(structure __s) {
      std::string chain_id_letters("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      std::vector<char> chain_id_alphabet(chain_id_letters.begin(), chain_id_letters.end());
      tools::mapper<std::string> m(std::vector<std::vector<char>>{chain_id_alphabet});
      std::string chain_id;
      for(size_t i(m.get_min()); has_structure(chain_id = m.encode(i)) && i <= m.get_max(); ++i) {
      } // for
      if(has_structure(chain_id)) {
        throw std::out_of_range("no letter was available to assign to this structure");
      } // if

      set(chain_id, __s);
      return chain_id;
    } // add()
    // set a structure to have a specific chain_id
    void assembly::set(std::string const &__chain_id, structure __s) { _chains[__chain_id] = __s; }
  } // namespace che
} // namespace biosim