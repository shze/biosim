#include "che/assembly.h"
#include "tools/mapper.h"

namespace biosim {
  namespace che {
    // default ctor
    assembly::assembly() : _chains() {}
    // ctor from structure
    assembly::assembly(structure __s) : _chains() { add(__s); }
    // ctor from structure_ensemble
    assembly::assembly(structure_ensemble __e) : _chains() { add(__e); }
    // get list of all chain_ids
    std::list<std::string> assembly::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto const &c : _chains) {
        chain_ids.push_back(c.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // get the size of the ensemble for the given chain id
    size_t assembly::get_ensemble_size(std::string const &__chain_id) const {
      return _chains.find(__chain_id) == _chains.end() ? 0 : _chains.at(__chain_id).get_samples().size();
    } // get_ensemble_size()
    // returns if an ensemble with the given chain id exists
    bool assembly::has_ensemble(std::string const &__chain_id) const { return get_ensemble_size(__chain_id) > 0; }
    // get the ensemble for the given chain id
    structure_ensemble const &assembly::get_ensemble(std::string const &__chain_id) const {
      return _chains.at(__chain_id);
    } // get_ensemble()
    // get the first structure of the ensemble; throws out of range if no structure with this chain_id exists
    structure const &assembly::get_first_structure(std::string const &__chain_id) const {
      return _chains.at(__chain_id).get_samples().front().sample;
    } // get_first_structure()
    // get all first structures
    std::vector<structure> assembly::get_first_structures() const {
      std::vector<structure> structures;
      structures.reserve(_chains.size());
      for(auto const &c : _chains) {
        structures.push_back(c.second.get_samples().front().sample);
      } // for
      return structures;
    } // get_first_structures()
    // add structure to assembly, assign next available chain_id, and returns the chain_id
    std::string assembly::add(structure __s) { return add(structure_ensemble(__s)); }
    // add ensemble to assembly, assign next available chain_id, and returns the chain_id
    std::string assembly::add(structure_ensemble __e) {
      std::string chain_id_letters("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
      std::vector<char> chain_id_alphabet(chain_id_letters.begin(), chain_id_letters.end());
      tools::mapper<std::string> m(std::vector<std::vector<char>>{chain_id_alphabet});
      std::string chain_id;
      for(size_t i(m.get_min()); has_ensemble(chain_id = m.encode(i)) && i <= m.get_max(); ++i) {
      } // for
      if(has_ensemble(chain_id)) {
        throw std::out_of_range("no letter was available to assign to this structure");
      } // if

      set(chain_id, __e);
      return chain_id;
    } // add()
    // set a structure to have a specific chain_id
    void assembly::set(std::string __chain_id, structure __s) { set(__chain_id, structure_ensemble(__s)); }
    // set an ensemble to have a specific chain_id
    void assembly::set(std::string __chain_id, structure_ensemble __e) {
      std::map<std::string, structure_ensemble>::iterator chain_id_itr(_chains.find(__chain_id));
      if(chain_id_itr != _chains.end()) { // see: http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3873.html
        chain_id_itr->second = __e; // update
      } // if
      else {
        _chains.emplace(__chain_id, __e); // insert
      } // else
    } // set()
  } // namespace che
} // namespace biosim