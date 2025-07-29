#include "frinZoutput.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio> // For popen, pclose

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void write_spectrum_to_text(
    const std::string& output_path,
    const HeaderRegion& header,
    const std::vector<std::vector<std::complex<float>>>& spectrum_data) {
    std::ofstream outfile(output_path);
    if (!outfile) {
        std::cerr << "Error: Cannot open file for writing: " << output_path << std::endl;
        return;
    }

    outfile << "# Header Information" << std::endl;
    outfile << "# Magic Word: 0x" << std::hex << header.magic_word << std::dec << std::endl;
    outfile << "# Header Version: " << header.header_version << std::endl;
    outfile << "# Software Version: " << header.software_version << std::endl;
    outfile << "# Sampling Speed: " << header.sampling_speed << " Hz" << std::endl;
    outfile << "# Observed Sky Freq: " << header.observed_sky_freq << " Hz" << std::endl;
    outfile << "# FFT Point: " << header.fft_point << std::endl;
    outfile << "# Number of Sectors in file: " << header.num_sector << std::endl;
    outfile << "# Station1 Name: " << std::string(header.station1_name, 16).c_str() << std::endl;
    outfile << "# Station2 Name: " << std::string(header.station2_name, 16).c_str() << std::endl;
    outfile << "# Source Name: " << std::string(header.source_name, 16).c_str() << std::endl;
    outfile << "# Source RA: " << header.source_ra << " rad" << std::endl;
    outfile << "# Source Dec: " << header.source_dec << " rad" << std::endl;
    outfile << "# --- Spectrum Data ---" << std::endl;
    outfile << "# Sector, FFT_Point, Real, Imag, Amplitude, Phase(deg)" << std::endl;

    for (size_t s = 0; s < spectrum_data.size(); ++s) {
        for (size_t i = 0; i < spectrum_data[s].size(); ++i) {
            const auto& data = spectrum_data[s][i];
            outfile << s << ", "
                    << i << ", "
                    << data.real() << ", "
                    << data.imag() << ", "
                    << std::abs(data) << ", "
                    << std::arg(data) * 180.0 / M_PI << std::endl;
        }
        outfile << std::endl; // Add a blank line after each sector
    }
}

void plot_spectrum_heatmaps_with_gnuplot(
    const std::string& output_path_amp,
    const std::string& output_path_phs,
    const std::vector<std::vector<std::complex<float>>>& spectrum_data,
    int fft_point,
    int num_sectors) {

    // Amplitude Heatmap
    FILE *gpipe_amp = popen("gnuplot -persist", "w");
    if (!gpipe_amp) {
        std::cerr << "Error: Could not open gnuplot pipe for amplitude plot." << std::endl;
        return;
    }

    fprintf(gpipe_amp, "set terminal pngcairo enhanced font 'sans,10' size 1024, 768\n");
    fprintf(gpipe_amp, "set output '%s'\n", output_path_amp.c_str());
    fprintf(gpipe_amp, "set title \"Complex Spectrum - Amplitude Heatmap\"\n");
    fprintf(gpipe_amp, "set xlabel \"FFT Point\"\n");
    fprintf(gpipe_amp, "set ylabel \"Sector\"\n");
    fprintf(gpipe_amp, "set zlabel \"Amplitude\"\n");
    fprintf(gpipe_amp, "set cblabel \"Amplitude\"\n");
    fprintf(gpipe_amp, "set pm3d map interpolate 1,1\n");
    fprintf(gpipe_amp, "set xrange [0:%d]\n", fft_point/2);
    fprintf(gpipe_amp, "set yrange [0:%d]\n", num_sectors);
    //fprintf(gpipe_amp, "set palette defined (0 \"blue\", 0.5 \"white\", 1 \"red\")\n");
    fprintf(gpipe_amp, "splot '-' using 2:1:3 with pm3d notitle\n"); // Read data from stdin

    for (int s = 0; s < num_sectors; ++s) {
        for (int i = 0; i < fft_point / 2; ++i) { // Assuming fft_point is total, and we use half
            const auto& data = spectrum_data[s][i];
            fprintf(gpipe_amp, "%d, %d, %f\n", s, i, std::abs(data));
        }
        fprintf(gpipe_amp, "\n"); // Blank line for pm3d
    }
    fflush(gpipe_amp);
    pclose(gpipe_amp);

    // Phase Heatmap
    FILE *gpipe_phs = popen("gnuplot -persist", "w");
    if (!gpipe_phs) {
        std::cerr << "Error: Could not open gnuplot pipe for phase plot." << std::endl;
        return;
    }

    fprintf(gpipe_phs, "set terminal pngcairo enhanced font 'sans,10' size 1024, 768\n");
    fprintf(gpipe_phs, "set output '%s'\n", output_path_phs.c_str());
    fprintf(gpipe_phs, "set title \"Complex Spectrum - Phase Heatmap\"\n");
    fprintf(gpipe_phs, "set xlabel \"FFT Point\"\n");
    fprintf(gpipe_phs, "set ylabel \"Sector\"\n");
    fprintf(gpipe_phs, "set zlabel \"Phase (degrees)\"\n");
    fprintf(gpipe_phs, "set cblabel \"Phase (degrees)\"\n");
    fprintf(gpipe_phs, "set cbrange [-180:180]\n");
    fprintf(gpipe_phs, "set pm3d map interpolate 1,1\n");
    fprintf(gpipe_phs, "set xrange [0:%d]\n", fft_point/2);
    fprintf(gpipe_phs, "set yrange [0:%d]\n", num_sectors);
    //fprintf(gpipe_phs, "set palette defined (0 \"blue\", 0.5 \"white\", 1 \"red\")\n");
    fprintf(gpipe_phs, "splot '-' using 2:1:3 with pm3d notitle\n"); // Read data from stdin

    for (int s = 0; s < num_sectors; ++s) {
        for (int i = 0; i < fft_point / 2; ++i) { // Assuming fft_point is total, and we use half
            const auto& data = spectrum_data[s][i];
            fprintf(gpipe_phs, "%d, %d, %f\n", s, i, std::arg(data) * 180.0 / M_PI);
        }
        fprintf(gpipe_phs, "\n"); // Blank line for pm3d
    }
    fflush(gpipe_phs);
    pclose(gpipe_phs);
}
