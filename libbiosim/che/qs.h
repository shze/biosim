#ifndef che_qs_h
#define che_qs_h

#include "che/ts.h"
#include "che/ss.h"

namespace biosim {
  namespace che {
    // quarternary structure of multiple ts with additional information
    class qs {
    public:
      // get list of all chain_ids
      std::list<std::string> get_chain_id_list() const;
      // get ts with given chain_id, throws out of range if no ts with this chain_id exists
      ts const &get(std::string __chain_id) const;
      // add ts to qs, assign next available chain_id, and returns the chain_id
      std::string add(ts __ts);

    private:
      std::map<std::string, ts> _ts; // maps chain_id -> ts
      ss _ss; // secondary structure
    }; // class qs
  } // namespace che
} // namespace biosim

#endif // che_qs_h
