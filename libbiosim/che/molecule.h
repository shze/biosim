#ifndef che_molecule_h
#define che_molecule_h

#include "che/ps.h"
#include "che/ts.h"
#include <string>
#include <algorithm>

namespace biosim {
  namespace che {
    // molecule contains the consitution (ps) and spatial positions (ts)
    class molecule {
    public:
      // ctor from storage and identifier (also default ctor)
      explicit molecule(std::string const &__storage = std::string(), std::string const &__identifier = std::string());
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
      ts _ts; // coordinates
    }; // class molecule

    // output operator for molecule
    inline std::ostream &operator<<(std::ostream &__out, molecule const &__molecule) {
      std::string id(__molecule.get_identifier());
      id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
      __out << __molecule.get_storage() << ": " << id << " length=" << __molecule.get_length() << "\n"
            << __molecule.get_ps();
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_molecule_h
