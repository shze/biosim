#Biosim

#Contributing code

##Formatting

* Format using `clang-format` and the style defined in `.clang-format`.  
  See http://clang.llvm.org/docs/ClangFormat.html on how to integrate formatting for different editors.  
  For KDevelop, add a custom source formatter script with `clang-format-3.4 -style=file -assume-filename=$FILE` as  
  command and assign your favorite shortcut key to the _Reformat Source_ command.
* All identifiers are lowercase, words are separated by `_`. Member variables have one, passed function variables  
  have two leading `_`. Template parameters are uppercase.
* Always place the `const` modifer after the part which constness it describes (so that types can be completely read  
  backwards).
* Group methods by method/design pattern without empty line between methods.
* Use short, but descriptive names; comment at least every method and member.

##Programming
