#ifndef score_function_h
#define score_function_h

#include <string>
#include "che/qs.h"

namespace biosim {
  namespace score {
    // base class for energy functions; design idea based on deprecated binary_function,
    // http://www.cplusplus.com/reference/functional/binary_function/ and
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2010/n3145.html
    template <class T>
    class energy_function {
    public:
      // returns score identifier; do not automate this by calling typeid() b/c name mangling is compiler specific
      // see e.g. http://stackoverflow.com/questions/3649278/how-can-i-get-the-class-name-from-a-c-object
      virtual std::string const &get_identifier() const = 0;
      // computes the energy for the given T
      virtual double energy(T const &) const = 0;
      // computes the energy for the given T; ensures energies are <= 0
      double operator()(T const &__t) const {
        double e(energy(__t));
        return e <= 0 ? e : 0;
      } // operator()
    }; // class energy_function

    // base class for distance functions
    template <class T>
    class distance_function {
      // returns score identifier
      virtual std::string const &get_identifier() const = 0;
      // computes the distance between the two given T
      virtual double distance(T const &, T const &) const = 0;
      // computes the distance between the two given T; ensures distances are >= 0
      double operator()(T const &__t1, T const &__t2) const {
        double d(distance(__t1, __t2));
        return d >= 0 ? d : 0;
      } // operator()
    }; // class distance_function

    // defines distance types
    using ss_distance = distance_function<che::ss>;
    using qs_distance = distance_function<che::qs>;
  } // namespace score
} // namespace biosim

#endif // score_function_h
