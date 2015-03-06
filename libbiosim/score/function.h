#ifndef score_function_h
#define score_function_h

#include <string>

namespace biosim {
  namespace score {
    // base class for energy functions; design idea based on deprecated binary_function,
    // http://www.cplusplus.com/reference/functional/binary_function/ and
    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2010/n3145.html
    template <class T>
    class en_function {
    public:
      // returns identifier; do not automate this by calling typeid() b/c name mangling is compiler specific
      // see e.g. http://stackoverflow.com/questions/3649278/how-can-i-get-the-class-name-from-a-c-object
      virtual std::string get_identifier() const = 0;
      // computes the energy for the given T
      virtual double energy(T const &) const = 0;
      // computes the energy for the given T
      double operator()(T const &__t) const { return energy(__t); } // operator()
    }; // class en_function

    // base class for comparison functions (distance, similarity, dissimilarity)
    template <class T>
    class cm_function {
    public:
      // returns identifier
      virtual std::string get_identifier() const = 0;
      // compares the given two instances of T
      virtual double compare(T const &, T const &) const = 0;
      // compares the given two instances of T
      double operator()(T const &__t1, T const &__t2) const { return compare(__t1, __t2); }
    }; // cm_function
  } // namespace score
} // namespace biosim

#endif // score_function_h
