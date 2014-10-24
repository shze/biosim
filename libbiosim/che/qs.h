#ifndef che_qs_h
#define che_qs_h

#include "che/ts.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // quarternary structure of multiple ts with additional information
    class qs {
    public:
      // default ctor
      qs();
      // ctor from ts and ss
      qs(ts __ts, ss __ss);
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // if the ts with the given chain_id exists
      bool has_ts(std::string __chain_id) const;
      // if the ss with the given chain_id exists
      bool has_ss(std::string __chain_id) const;
      // get ts with given chain_id, throws out of range if no ts with this chain_id exists
      ts const &get_ts(std::string __chain_id) const;
      // get ss with given chain_id, throws out of range if no ts with this chain_id exists
      ss const &get_ss(std::string __chain_id) const;
      // add ts to qs, assign next available chain_id, and returns the chain_id
      std::string add(ts __ts);
      // add corresponding ts and ss to qs; assign next available chain_id, and returns the chain_id
      std::string add(ts __ts, ss __ss);

    private:
      std::map<std::string, ts> _ts; // maps chain_id -> ts
      std::map<std::string, ss> _ss; // maps chain_id -> ss
    }; // class qs

    // output operator for qs
    inline std::ostream &operator<<(std::ostream &__out, qs const &__qs) {
      for(auto const &chain_id : __qs.get_chain_id_list()) {
        __out << "Chain: " << chain_id << "\n";
        if(__qs.has_ts(chain_id)) {
          __out << __qs.get_ts(chain_id);
        } // if
        if(__qs.has_ss(chain_id)) {
          __out << __qs.get_ss(chain_id);
        } // if
      } // for
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_qs_h
