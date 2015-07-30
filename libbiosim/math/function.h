#ifndef math_function_h
#define math_function_h

#include <string>
#include <stdexcept>
#include <memory>

namespace biosim {
  namespace math {
    // base class for evaluation functions, e.g. energy functions, probability functions;
    // design idea based on deprecated binary_function, http://www.cplusplus.com/reference/functional/binary_function/
    // and http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2010/n3145.html
    template <class T>
    class ev_function {
    public:
      // returns identifier; do not automate this by calling typeid() b/c name mangling is compiler specific
      // see e.g. http://stackoverflow.com/questions/3649278/how-can-i-get-the-class-name-from-a-c-object
      virtual std::string get_identifier() const = 0;
      // computes the energy for the given T
      virtual double evaluate(T const &) const = 0;
      // computes the energy for the given T
      double operator()(T const &__t) const { return evaluate(__t); }
    }; // class ev_function

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
    }; // class cm_function

    // define identifier function for enumerate as template for all cm_function ptrs
    template <class T>
    std::string identifier(std::shared_ptr<cm_function<T>> const &__cm_function_ptr) {
      return __cm_function_ptr->get_identifier();
    } // identifier()

    // exception for compare function errors
    struct compare_error : std::runtime_error {
      explicit compare_error(std::string const &__msg)
          : std::runtime_error("Cannot calculate compare score: " + __msg) {}
    }; // struct compare_error
  } // namespace math
} // namespace biosim

#endif // math_function_h
