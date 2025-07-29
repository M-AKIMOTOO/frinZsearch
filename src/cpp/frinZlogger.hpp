#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip> // For std::setprecision

class Logger {
public:
    std::ostream* console_stream = nullptr;
    std::ofstream file_stream;
    bool file_output_enabled = false;
    bool console_output_enabled = true;

    Logger(); // Default constructor

    void setup(bool enable_console, bool enable_file, const std::string& file_path);

    template<typename T>
    Logger& operator<<(const T& val);

    // Typedefs for manipulators
    typedef std::ostream& (*StandardEndlManipulator)(std::ostream&);
    typedef std::ios_base& (*StandardFlagManipulator)(std::ios_base&);

    Logger& operator<<(StandardEndlManipulator manip); // For std::endl
    Logger& operator<<(StandardFlagManipulator manip); // For std::hex, std::dec, std::fixed

    // Method for std::setprecision - needs to be handled carefully or rethought
    // For simplicity, you might pass std::setprecision directly to the stream
    // or create a wrapper if you really need to chain it like `logger.setprecision(5) << ...`
    // This current approach in frinZsearch.cpp for setprecision is a bit unusual.
    // A more common way is: logger << std::setprecision(5) << my_float;
    // However, to keep your existing Logger interface:
    Logger& setprecision(int n);


    ~Logger();
};

// Template implementation must be in the header or an included .tpp/.ipp file
template<typename T>
Logger& Logger::operator<<(const T& val) {
    if (console_output_enabled && console_stream) {
        *console_stream << val;
    }
    if (file_output_enabled && file_stream.is_open()) {
        file_stream << val;
    }
    return *this;
}

#endif // LOGGER_HPP