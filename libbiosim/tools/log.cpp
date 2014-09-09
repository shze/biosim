#include "tools/log.h"
#include <cstdio>
#include <cstring>

namespace biosim {
  namespace tools {
    // (static) get current output level
    size_t log_message::get_output_level() { return _current_output_level; }
    // (static) set current output level
    void log_message::set_output_level(size_t const &__level) {
      if(__level == level_error || __level == level_info || __level == level_debug) { // only set defined levels
        _current_output_level = __level;
      } // if
    } // set_output_level

    // ctor; used in macro; no default constructor needed
    log_message::log_message(char const *__file, int const __line, level_type const __level)
        : _level(__level), _string_stream() {
      _string_stream << "(" << _level_strings[_level];
      if(_current_output_level >= level_debug) { // print additional information
        _string_stream << ":" << basename(__file) << ":" << __line;
      } // if
      _string_stream << ") ";
    } // ctor
    // dtor
    log_message::~log_message() {
      // only print message if level of message is at least (i.e. equal or lower than) current output level
      if(_level <= _current_output_level) {
        _string_stream << "\n"; // add a final line break

        // avoid cerr, dtor may get called during exit code, and cerr may be partially or fully destroyed by then
        fwrite(_string_stream.str().c_str(), _string_stream.str().size(), 1, stderr);
      } // if
    } // dtor
    // returns the stream to which the message is written
    std::ostream &log_message::get_stream() { return _string_stream; }

    // (static) returns the basename of the given filename
    char const *log_message::get_basename(char const *__file_path) {
      char const *base(strrchr(__file_path, '/'));
      return base ? (base + 1) : __file_path;
    } // get_basename()

    char const *log_message::_level_strings[] = {"EE", "II", "DD"}; // initialize the message level strings
    size_t log_message::_current_output_level = log_message::level_info; // default message level
  } // namespace tools
} // namespace biosim
