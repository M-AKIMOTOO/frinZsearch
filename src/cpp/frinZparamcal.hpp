#ifndef FRINZPARAMCAL_HPP
#define FRINZPARAMCAL_HPP

#include <vector>
#include <string> // For error messages or detailed results
#include <cmath>  // For std::sqrt, std::max
#include <utility> // For std::pair
#include <fftw3.h> // For fftwf_complex

// Structure to hold the calculated parameters
struct FftPeakParameters {
    float max_amplitude = -1.0f;
    int max_i = -1; // Row index of max amplitude (after fftshift)
    int max_j = -1; // Column index of max amplitude (after fftshift)
    float mean_amplitude = 0.0f;
    float stddev_amplitude = 0.0f;
    float snr = 0.0f;
    float physical_rate_hz = 0.0f;
    float physical_delay_samples = 0.0f;
    bool success = false;
    std::string message;
};

// Calculates mean and standard deviation of a 2D float vector
std::pair<float, float> calculate_mean_stddev_from_float_vec(const std::vector<std::vector<float>>& data);

// Calculates peak parameters from a 2D fft-shifted amplitude array
FftPeakParameters calculate_fft_peak_parameters(
    const std::vector<std::vector<float>>& fft_shifted_amplitude,
    const fftwf_complex* fft_out_complex_data, // FFT output (unshifted) for noise statistics
    float effective_integration_length, // For physical rate calculation
    int N_rows_padded,                  // Total rows in shifted_amplitude (for rate calc)
    int N_cols_padded_fft,              // Total cols in shifted_amplitude (for delay calc)
    // Added parameters for search range
    double delay_search_min_param,
    double delay_search_max_param,
    bool delay_search_range_specified_param,
    double rate_search_min_param,
    double rate_search_max_param,
    bool rate_search_range_specified_param
);

#endif // FRINZPARAMCAL_HPP
