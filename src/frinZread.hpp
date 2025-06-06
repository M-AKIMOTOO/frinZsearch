#ifndef FRINZREAD_HPP
#define FRINZREAD_HPP

#include <string>
#include <vector>
#include <complex>
#include <cstdint>
#include <fstream> // For std::ifstream, already present
#include "frinZlogger.hpp" // Include the new Logger header

#pragma pack(push, 1)
struct HeaderRegion {
    uint32_t magic_word;
    uint32_t header_version;
    uint32_t software_version;
    uint32_t sampling_speed;
    double observed_sky_freq;
    uint32_t fft_point;
    uint32_t num_sector;
    char station1_name[16];
    double station1_pos_x;
    double station1_pos_y;
    double station1_pos_z;
    char station1_key[8];
    char station2_name[16];
    double station2_pos_x;
    double station2_pos_y;
    double station2_pos_z;
    char station2_key[8];
    char source_name[16];
    double source_ra;
    double source_dec;
};
#pragma pack(pop)

struct FileData {
    HeaderRegion header;
    std::vector<std::vector<std::complex<float>>> spectrum_data;
    std::vector<std::string> sector_start_times_utc;
    float first_effective_integration_length = 0.0f;
    bool success = false;
    std::string error_message;
};

// Function to read binary data from the file
FileData read_binary_file_data(const std::string& filepath, Logger& logger, bool output_header_info_flag);

#endif // FRINZREAD_HPP
