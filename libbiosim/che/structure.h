#ifndef che_structure_h
#define che_structure_h

#include "che/sequence.h"
#include "che/ss.h"
#include <string>
#include <algorithm>

namespace biosim {
  namespace che {
    // structure contains the primary structure (sequence, positions) and the secondary structure
    class structure {
    public:
      // ctor from ps (also default ctor)
      explicit structure(ps __ps = ps());
      // ctor from storage and identifier
      structure(std::string __storage, std::string __identifier);
      // ctor from storage, identifier, and ps
      structure(std::string __storage, std::string __identifier, ps __ps);
      // ctor from storage, identifier, ps, and ss
      structure(std::string __storage, std::string __identifier, ps __ps, ss __ss);
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
      ps _ps; // sequence
      boost::optional<std::string> _storage; // optional persistent storage
      boost::optional<std::string> _id; // optional identifier (eg. fasta identifier)
      boost::optional<ss> _ss; // optional secondary structure
    }; // class structure

    // output operator for structure
    inline std::ostream &operator<<(std::ostream &__out, structure const &__m) {
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

#endif // che_structure_h
