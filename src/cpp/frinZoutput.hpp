#ifndef FRINZOUTPUT_HPP
#define FRINZOUTPUT_HPP

#include <string>
#include <vector>
#include <complex>
#include "frinZread.hpp"

// Function to write spectrum data and header to a text file
void write_spectrum_to_text(
    const std::string& output_path,
    const HeaderRegion& header,
    const std::vector<std::vector<std::complex<float>>>& spectrum_data);

// Function to plot amplitude and phase heatmaps using gnuplot
void plot_spectrum_heatmaps_with_gnuplot(
    const std::string& output_path_amplitude,
    const std::string& output_path_phase,
    const std::vector<std::vector<std::complex<float>>>& spectrum_data,
    int fft_point,
    int num_sectors);

#endif // FRINZOUTPUT_HPP
