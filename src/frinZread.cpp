#include "frinZread.hpp"
#include <iostream> // For std::cerr in case logger is not fully set up
#include <iomanip>  // For std::hex, std::dec, std::fixed, std::setprecision
#include <cmath>    // For std::round
#include <ctime>    // For std::gmtime, std::strftime

// Note: Logger is currently included via frinZread.hpp -> frinZsearch.cpp
// This is not ideal and Logger should be in its own .hpp/.cpp pair.

FileData read_binary_file_data(const std::string& filepath, Logger& logger, bool output_header_info_flag) {
    FileData result;
    result.success = false; // Default to failure

    std::ifstream infile(filepath, std::ios::binary);
    if (!infile) {
        result.error_message = "Error opening file: " + filepath;
        // Use logger if available, otherwise cerr
        if (logger.console_output_enabled && logger.console_stream) {
             logger << result.error_message << std::endl;
        } else {
            std::cerr << result.error_message << std::endl;
        }
        return result;
    }

    infile.seekg(0);
    infile.read(reinterpret_cast<char*>(&result.header), sizeof(HeaderRegion));
    if (infile.gcount() != sizeof(HeaderRegion)) {
        result.error_message = "Error reading header from file: " + filepath;
        if (logger.console_output_enabled && logger.console_stream) {
            logger << result.error_message << std::endl;
        } else {
            std::cerr << result.error_message << std::endl;
        }
        infile.close();
        return result;
    }

    if (output_header_info_flag) {
        if (!logger.console_output_enabled || (logger.console_output_enabled && logger.console_stream)) { // Check if console is enabled
            logger << "=== Header Info (from frinZread) ===\n";
            logger << "Magic Word: 0x" << std::hex << result.header.magic_word << std::dec << "\n";
            logger << "Header Version: " << result.header.header_version << "\n";
            logger << "Software Version: " << result.header.software_version << "\n";
            logger << "Sampling Speed: " << result.header.sampling_speed << " Hz\n";
            logger << "Observed Sky Freq: " << result.header.observed_sky_freq << " Hz\n";
            logger << "FFT Point: " << result.header.fft_point << "\n";
            logger << "Number of Sectors in file: " << result.header.num_sector << "\n";
            logger << "Station1 Name: " << std::string(result.header.station1_name, 16).c_str() << "\n";
            logger << "Station1 Position: (" << result.header.station1_pos_x << ", " << result.header.station1_pos_y << ", " << result.header.station1_pos_z << ")\n";
            logger << "Station1 Key: " << std::string(result.header.station1_key, 8).c_str() << "\n";
            logger << "Station2 Name: " << std::string(result.header.station2_name, 16).c_str() << "\n";
            logger << "Station2 Position: (" << result.header.station2_pos_x << ", " << result.header.station2_pos_y << ", " << result.header.station2_pos_z << ")\n";
            logger << "Station2 Key: " << std::string(result.header.station2_key, 8).c_str() << "\n";
            logger << "Source Name: " << std::string(result.header.source_name, 16).c_str() << "\n";
            logger << "Source RA: " << result.header.source_ra << " rad\n";
            logger << "Source Dec: " << result.header.source_dec << " rad" << std::endl;
        }
    }

    if (result.header.fft_point == 0 || result.header.num_sector == 0) {
        result.error_message = "Error: Invalid header data (fft_point or num_sector is zero).";
         if (logger.console_output_enabled && logger.console_stream) {
            logger << result.error_message << std::endl;
        } else {
            std::cerr << result.error_message << std::endl;
        }
        infile.close();
        return result;
    }

    result.spectrum_data.reserve(result.header.num_sector);
    result.sector_start_times_utc.reserve(result.header.num_sector);
    bool first_sector_processed = false;

    for (uint32_t s = 0; s < result.header.num_sector; ++s) {
        std::streamoff sector_offset = 256 + s * (128 + 4 * result.header.fft_point);
        infile.seekg(sector_offset, std::ios::beg);

        uint32_t corr_start_sec_raw = 0;
        infile.read(reinterpret_cast<char*>(&corr_start_sec_raw), sizeof(uint32_t));
        std::time_t current_epoch_time = static_cast<std::time_t>(corr_start_sec_raw);
        
        std::tm* gmt_time = std::gmtime(&current_epoch_time);
        char time_buffer[80];
        char final_buffer[150];
        if (gmt_time) {
            std::strftime(time_buffer, sizeof(time_buffer), "%Y-%j %H:%M:%S UTC", gmt_time);
        }
    
        infile.seekg(108, std::ios::cur);

        float effective_integration_length = 0.0f;
        infile.read(reinterpret_cast<char*>(&effective_integration_length), sizeof(float));
        if (!first_sector_processed) {
            result.first_effective_integration_length = effective_integration_length;
            first_sector_processed = true;
        }

        if (gmt_time) {
            std::snprintf(final_buffer, sizeof(final_buffer), "%s (Epoch: %lu) (Int.Length: %.6f s)", 
                          time_buffer, static_cast<unsigned long>(current_epoch_time), effective_integration_length);
            result.sector_start_times_utc.push_back(final_buffer);
        } else {
            std::snprintf(final_buffer, sizeof(final_buffer), "Error converting epoch %lu (Int.Length: %.6f s)", 
                          static_cast<unsigned long>(current_epoch_time), effective_integration_length);
            result.sector_start_times_utc.push_back(final_buffer);
        }

        infile.seekg(12, std::ios::cur);

        std::vector<std::complex<float>> current_sector_spectrum;
        current_sector_spectrum.reserve(result.header.fft_point / 2);
        for (uint32_t i = 0; i < result.header.fft_point / 2; ++i) {
            float real = 0.0f, imag = 0.0f;
            infile.read(reinterpret_cast<char*>(&real), sizeof(float));
            infile.read(reinterpret_cast<char*>(&imag), sizeof(float));
            current_sector_spectrum.emplace_back(real, imag);
        }
        result.spectrum_data.push_back(current_sector_spectrum);
    }

    infile.close();
    result.success = true;
    return result;
}
