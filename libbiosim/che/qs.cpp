#include "che/qs.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // get list of all chain_ids
    std::list<std::string> qs::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto p : _ts) {
        chain_ids.push_back(p.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // get molecule with given chain_id, throws out of range if no molecule with this chain_id exists
    ts const &qs::get(std::string __chain_id) const { return _ts.at(__chain_id); }
    // add molecule to ts, assign next available chain_id, and returns the chain_id
    std::string qs::add(ts __ts) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor inc;
      while(_ts.find(chain_id) != _ts.end()) {
        chain_id = inc.next(chain_id);
      } // while

      _ts[chain_id] = __ts;
      return chain_id;
    } // add()
  } // namespace che
} // namespace biosim