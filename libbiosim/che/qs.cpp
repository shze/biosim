#include "che/qs.h"
#include "tools/incrementor.h"

namespace biosim {
  namespace che {
    // default ctor
    qs::qs() : _ts(), _ss() {}
    // ctor from ts and ss
    qs::qs(ts __ts, ss __ss) : _ts(), _ss() { add(__ts, __ss); }
    // get list of all chain_ids
    std::list<std::string> qs::get_chain_id_list() const {
      std::list<std::string> chain_ids;
      for(auto p : _ts) {
        chain_ids.push_back(p.first);
      } // for
      return chain_ids;
    } // get_chain_id_list()
    // if the ts with the given chain_id exists
    bool qs::has_ts(std::string const &__chain_id) const { return _ts.find(__chain_id) != _ts.end(); }
    // if the ss with the given chain_id exists
    bool qs::has_ss(std::string const &__chain_id) const { return _ss.find(__chain_id) != _ss.end(); }
    // get ts with given chain_id, throws out of range if no ts with this chain_id exists
    ts const &qs::get_ts(std::string const &__chain_id) const { return _ts.at(__chain_id); }
    // get ss with given chain_id, throws out of range if no ts with this chain_id exists
    ss const &qs::get_ss(std::string const &__chain_id) const { return _ss.at(__chain_id); }
    // add molecule to ts, assign next available chain_id, and returns the chain_id
    std::string qs::add(ts __ts) {
      std::string chain_id = "A"; // determine next available chain_id, start from A
      tools::incrementor inc;
      while(_ts.find(chain_id) != _ts.end()) {
        chain_id = inc.next(chain_id);
      } // while

      _ts.insert(std::pair<std::string, ts>(chain_id, __ts));
      return chain_id;
    } // add()
    // add corresponding ts and ss to qs; assign next available chain_id, and returns the chain_id
    std::string qs::add(ts __ts, ss __ss) {
      std::string chain_id = add(__ts);
      _ss.insert(std::pair<std::string, ss>(chain_id, __ss));
      return chain_id;
    } // add()
    // set a ts and ss to have a specific chain_id
    void qs::set(std::string const &__chain_id, ts __ts, ss __ss) {
      _ts[__chain_id] = __ts;
      // for ss use erase() and insert(), because operator[] requires a default ctor which ss does not have
      _ss.erase(__chain_id);
      _ss.insert(std::pair<std::string, ss>(__chain_id, __ss));
    } // set()
  } // namespace che
} // namespace biosim