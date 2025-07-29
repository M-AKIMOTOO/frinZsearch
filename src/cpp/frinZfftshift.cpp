#include "frinZfftshift.hpp"
#include <cmath> // For std::sqrt

std::vector<std::vector<float>> perform_fftshift_and_calc_amplitude(
    const fftwf_complex* fft_output_data,
    int N_rows,
    int N_cols) {

    std::vector<std::vector<float>> shifted_amplitude(N_rows, std::vector<float>(N_cols));

    if (!fft_output_data || N_rows <= 0 || N_cols <= 0) {
        // Return empty or appropriately sized zeroed vector if input is invalid
        return shifted_amplitude; // Or throw an exception
    }

    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            int shifted_i = (i + N_rows / 2) % N_rows;
            int shifted_j = (j + N_cols / 2) % N_cols;

            float real_part = fft_output_data[i * N_cols + j][0];
            float imag_part = fft_output_data[i * N_cols + j][1];
            
            shifted_amplitude[shifted_i][shifted_j] = std::sqrt(real_part * real_part + imag_part * imag_part);
        }
    }
    return shifted_amplitude;
}
