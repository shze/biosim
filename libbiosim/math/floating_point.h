#ifndef math_floating_point_h
#define math_floating_point_h

#include "tools/type_with_size.h"

namespace biosim {
  namespace math {
    // class for comparing and manipulating floating point numbers
    // design:
    // it could be designed around boost::math::float_distance(); however, this requires inf and nan to be handled
    // separately in almost_equal(), because float_distance() throws on either encounter;
    // since with this approach there is no way to manipulate the internal representation to allow for testing, this
    // class was based on https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    template <class T>
    class floating_point {
    public:
      // max number of ulp (units in last place) difference tolerated between two numbers to still be almost equal;
      // value taken from the google test framework, see:
      // https://code.google.com/p/googletest/source/browse/trunk/include/gtest/internal/gtest-internal.h
      // for ulps, see: https://en.wikipedia.org/wiki/Unit_in_the_last_place
      static size_t get_max_ulp_tolerance() { return 4; }

      using bits = typename tools::type_with_size<sizeof(T)>::uint; // simplify naming

      // get the number of bits of T (constexpr to allow its use in the definition of the bitfield)
      constexpr static size_t get_type_bit_count() { return 8 * sizeof(T); }
      // get the number of bits in T's mantissa
      constexpr static size_t get_mantissa_bit_count() { return std::numeric_limits<T>::digits - 1; }
      // get the number of bits in T's exponent
      constexpr static size_t get_exponent_bit_count() { return get_type_bit_count() - get_mantissa_bit_count() - 1; }
      // get the mask for the sign bit
      static bits get_sign_bit_mask() { return static_cast<bits>(1) << (get_type_bit_count() - 1); }

      // ctor from T; default ctor
      explicit floating_point(T __value = 0.0) : value(__value) {}

      // returns true iff both numbers are at most get_max_ulp_tolerance() ULP's away from each other
      bool almost_equal(floating_point const &__rhs) const {
        // make sure to return false if either compared value is nan
        return !std::isnan(value) && !std::isnan(__rhs.value) &&
               bits_distance(representation, __rhs.representation) <= get_max_ulp_tolerance();
      } // almost_equal()

      // bit fields (inside the struct) cannot be accessed by non-const ref; make members public to allow write access
      union {
        T value; // original floating point data type
        bits representation; // internal bit representation
        struct {
          bits mantissa : get_mantissa_bit_count(); // mantissa part
          bits exponent : get_exponent_bit_count(); // exponent part
          bits sign : 1; // sign part
        }; // struct
      }; // union

    private:
      // convert the bit representation into an ascending representation
      static bits to_ascending_representation(bits const &__bits) {
        // if branch: __bits represents negative number; else branch: __bits represents positive number
        return get_sign_bit_mask() & __bits ? ~__bits + 1 : (get_sign_bit_mask() | __bits);
      } // to_ascending_representation()

      // calculate the difference between two bit representations
      static bits bits_distance(bits const &__bits1, bits const &__bits2) {
        bits const asc1(to_ascending_representation(__bits1)), asc2(to_ascending_representation(__bits2));
        return (asc1 >= asc2) ? (asc1 - asc2) : (asc2 - asc1);
      } // bits_distance()
    }; // class floating_point
  } // namespace math
} // namespace biosim

#endif // math_floating_point_h
