#ifndef tools_less_h
#define tools_less_h

namespace biosim {
  namespace tools {
    // class to use any comparison function as functor; for the function template parameter see:
    // http://bytes.com/topic/c/answers/138989-class-template-parameter-function
    template <class T, bool F(T, T)>
    class less {
    public:
      // function operator taking two instances of type T and using F for the comparison
      bool operator()(T __instance1, T __instance2) { return F(__instance1, __instance2); }
    }; // class less
  } // namespace tools
} // namespace biosim

#endif // math_less_h
