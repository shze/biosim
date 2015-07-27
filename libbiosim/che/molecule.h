#ifndef che_molecule_h
#define che_molecule_h

#include "che/sequence.h"
#include "che/ss.h"
#include <string>
#include <algorithm>

namespace biosim {
  namespace che {
    // molecule contains the consitution (sequence) and spatial positions
    class molecule {
    public:
      // ctor from storage and identifier (also default ctor)
      explicit molecule(std::string __storage = std::string(), std::string __identifier = std::string());
      // ctor from storage, identifier, and ps
      molecule(std::string __storage, std::string __identifier, ps __ps);
      // ctor from storage, identifier, ps, and ss
      molecule(std::string __storage, std::string __identifier, ps __ps, ss __ss);
      // get storage
      std::string const &get_storage() const;
      // get identifier
      std::string const &get_identifier() const;
      // get length
      size_t get_length() const;
      // get ps
      ps const &get_ps() const;
      // get ss
      ss const &get_ss() const;

    private:
      std::string _storage; // persistent storage
      std::string _id; // identifier (eg. fasta identifier)
      ps _ps; // sequence
      ss _ss; // secondary structure
    }; // class molecule

    // output operator for molecule
    inline std::ostream &operator<<(std::ostream &__out, molecule const &__m) {
      std::string id(__m.get_identifier());
      id.erase(std::remove(id.begin(), id.end(), '\n'), id.end());
      __out << __m.get_storage() << ": " << id << " length=" << __m.get_length() << "\n" << __m.get_ps() << "\n";
      if(__m.get_ss().defined()) {
        __out << __m.get_ss(); // no newline here, b/c ss data consists of lines
      } // if
      return __out;
    } // operator<<()
  } // namespace che
} // namespace biosim

#endif // che_molecule_h
