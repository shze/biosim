namespace biosim {
  namespace tools {
    // enumerate is wrapper of T; it uses identifier() to get an identifier string for a T
    template <typename T>
    class enumerate {
    public:
      // constructor
      explicit enumerate(T __obj) : _object(__obj) {}
      // get object held by enumerate; use this explicit name instead of overloading operator * or operator T
      T const &get_object() const { return _object; }

      using instance_container = std::map<std::string, enumerate<T>>;

    private:
      // returns changeable ref to static instance_container
      static instance_container &get_instances_priv() {
        static instance_container instances;
        return instances;
      } // get_instances_priv()

    public:
      // returns const ref to static instance_container
      static instance_container const &get_instances() { return get_instances_priv(); }
      // adds new key identifier(__obj)->data enumerate to the instance_container
      static void add(T __obj) {
        get_instances_priv().insert(std::pair<std::string, enumerate<T>>(identifier(__obj), enumerate(__obj)));
      } // add()

    private:
      T _object; // the enumerated object
    };           // class enumerate

    // output operator for enumerate<T>
    template <typename T>
    inline std::ostream &operator<<(std::ostream &__out, enumerate<T> const &__rhs) {
      __out << __rhs.get_object();
      return __out;
    } // operator<<
  }
}
