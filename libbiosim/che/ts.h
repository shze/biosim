#ifndef che_ts_h
#define che_ts_h

#include "che/sequence.h"
#include <string>
#include <algorithm>

namespace biosim {
  namespace che {
    // teriary structure contains the consitution (sequence) and spatial positions
    class ts {
    public:
      // ctor from storage and identifier (also default ctor)
      explicit ts(std::string const &__storage = std::string(), std::string const &__identifier = std::string());
      // get storage
      std::string const &get_storage() const;
      // get identifier
      std::string const &get_identifier() const;
      // get length
      size_t get_length() const;
      // get ps
      ps const &get_ps() const;
      // set ps
      void set_ps(ps __ps);

    private:
      std::string _storage; // persistent storage
      std::string _id; // identifier (eg. fasta identifier)
      ps _ps; // sequence
    }; // class ts

    // output operator for molecule
    inline std::ostream &operator<<(std::ostream &__out, ts const &__ts) {
      std::string id(__ts.get_identifier());
      id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
      __out << __ts.get_storage() << ": " << id << " length=" << __ts.get_length() << "\n" << __ts.get_ps();
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_molecule_h
