#ifndef tools_compare_h
#define tools_compare_h

namespace biosim {
  namespace tools {
    // class to use any comparison function as functor; for the function template parameter see:
    // http://bytes.com/topic/c/answers/138989-class-template-parameter-function
    template <class T, bool F(T const &, T const &)>
    class compare {
    public:
      // function operator taking two instances of type T and using F for the comparison
      bool operator()(T const &__instance1, T const &__instance2) { return F(__instance1, __instance2); }
    }; // class compare
  } // namespace tools
} // namespace biosim

#endif // tools_compare_h
