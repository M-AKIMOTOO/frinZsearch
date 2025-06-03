#include "frinZlogger.hpp"

Logger::Logger() : console_stream(nullptr), file_output_enabled(false), console_output_enabled(true) {
    // Default constructor can be empty if members are initialized as above
}

void Logger::setup(bool enable_console, bool enable_file, const std::string& file_path) {
    console_output_enabled = enable_console;
    if (enable_console) {
        console_stream = &std::cout;
    } else {
        console_stream = nullptr; // Ensure it's null if console is disabled
    }

    if (enable_file && !file_path.empty()) {
        // Close if already open, to handle multiple setups (though not ideal)
        if (file_stream.is_open()) {
            file_stream.close();
        }
        file_stream.open(file_path);
        if (file_stream.is_open()) {
            file_output_enabled = true;
        } else {
            std::cerr << "Error: Could not open log file: " << file_path << std::endl;
            file_output_enabled = false;
        }
    } else {
        file_output_enabled = false;
    }
}

Logger& Logger::operator<<(StandardEndlManipulator manip) {
    if (console_output_enabled && console_stream) manip(*console_stream);
    if (file_output_enabled && file_stream.is_open()) manip(file_stream);
    return *this;
}

Logger& Logger::operator<<(StandardFlagManipulator manip) {
    if (console_output_enabled && console_stream) manip(*console_stream);
    if (file_output_enabled && file_stream.is_open()) manip(file_stream);
    return *this;
}

Logger& Logger::setprecision(int n) {
    if (console_output_enabled && console_stream) *console_stream << std::setprecision(n);
    if (file_output_enabled && file_stream.is_open()) file_stream << std::setprecision(n);
    return *this;
}

Logger::~Logger() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}