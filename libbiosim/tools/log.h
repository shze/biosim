#ifndef tools_log_h
#define tools_log_h

#include <sstream>

// LOG macros: used like std::cout/cerr with the insert operator (<<)
#define LOG biosim::tools::log_message(__FILE__, __LINE__).get_stream()
#define ERROR biosim::tools::log_message(__FILE__, __LINE__, biosim::tools::log_message::level_type::error).get_stream()
#define DEBUG biosim::tools::log_message(__FILE__, __LINE__, biosim::tools::log_message::level_type::debug).get_stream()

namespace biosim {
  namespace tools {
    // log output to std::cout with additional information like file or line number; uses internal buffer to store the
    // message before outputting and fwrite to output it; a message is only printed if its priority is higher than the
    // current global level (i.e. value is equal or lower); very much simplified from Google Logging
    class log_message {
    public:
      enum class level_type { error, info, debug }; // message levels

      // ctor; used in macro; no default constructor needed
      log_message(char const *__file, int const __line, level_type const __level = level_type::info);
      // dtor
      ~log_message();
      // returns the stream to which the message is written
      std::ostream &get_stream();

      // get current output level
      static level_type get_output_level();
      // set current output level
      static void set_output_level(level_type const &__level);

    private:
      // returns the basename of the given filename
      static char const *get_basename(char const *__file_path);

      level_type _level; // message level for this message
      std::ostringstream _string_stream; // buffer to collect message string

      static char const *_level_strings[]; // const strings for the names of the different levels
      static level_type
          _current_output_level; // current level, of level_type; messages with level <= current are printed
    }; // class log_message
  } // namespace tools
} // namespace biosim

#endif // tools_log_h
