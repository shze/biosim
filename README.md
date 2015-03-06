#Biosim

#Contributing code

##Formatting

* Format using `clang-format` and the style defined in `.clang-format`.  
  Use `clang-format-3.5` or newer because they can fit short enums on a single line.  
  See http://clang.llvm.org/docs/ClangFormat.html on how to integrate formatting for different editors.
  For KDevelop, add a custom source formatter script with `clang-format-3.4 -style=file -assume-filename=$FILE` as
  command and assign your favorite shortcut key to the _Reformat Source_ command.
* All identifiers are lowercase, words are separated by `_`. Member variables have one, passed function variables
  have two leading `_`. Template parameters are uppercase.
* Always place the `const` modifer after the part which constness it describes (so that types can be completely read
  backwards).
* Always use `{}`, even for single statements.
* Group methods by method/design pattern without empty line between methods.
* Use short, but descriptive names; comment at least every method and member.

##Programming

* Use exceptions, derive from std::runtime_exception; don't use boost::exceptions, they don't mix well, 
  see: http://stackoverflow.com/questions/25759293/
* Make single argument ctor explicit.
* Use enum class instead of enum for simple enumerations to create identifiers in separate namespace.

##Code Design

* Applications follow the Unix philosophy of small, focused tools.
* To reduce code size, the biosim lib is not linked tp applications with --whole-archive (gcc) or -force_load (clang). 
  Thus self-registering static instances cannot be used; they would be removed during linking optimizations, b/c they 
  are not explicitly called, 
  see: http://stackoverflow.com/questions/805555/; http://stackoverflow.com/questions/16294842/; 
  or http://stackoverflow.com/questions/17470350/.__
  Therefore, each enumerated class requires two changes: class and explicit enumeration.
