
#ifndef FRINZFFTSHIFT_HPP
#define FRINZFFTSHIFT_HPP

#include <vector>
#include <complex> // For fftwf_complex, though it's a C-style struct
#include <fftw3.h> // For fftwf_complex type

// Performs FFT shift on 2D FFT output and calculates amplitude
std::vector<std::vector<float>> perform_fftshift_and_calc_amplitude(
    const fftwf_complex* fft_output_data, // Raw FFTW output
    int N_rows,                           // Number of rows in the FFT output
    int N_cols                            // Number of columns in the FFT output (padded for FFT)
);

#endif // FRINZFFTSHIFT_HPP
