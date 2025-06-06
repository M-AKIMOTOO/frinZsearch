#include "frinZparamcal.hpp"
#include <numeric> // For std::accumulate (though not used in current mean/stddev)
#include <algorithm> // For std::max_element (if needed, not used here for 2D)
#include <fftw3.h> // For fftwf_complex
#include <cmath>   // For std::sqrt, std::max
#include <complex> // For std::complex, std::abs
#include <tuple>   // For std::tie


std::pair<float, float> calculate_mean_stddev_from_float_vec(const std::vector<std::vector<float>>& data) {
    double sum = 0.0;
    double sum_sq = 0.0;
    size_t count = 0;
    for (const auto& row : data) {
        for (float val : row) {
            sum += val;
            sum_sq += val * val;
            count++;
        }
    }
    if (count == 0) return {0.0f, 0.0f};
    double mean = sum / count;
    double variance = (sum_sq / count) - (mean * mean);
    return {static_cast<float>(mean), static_cast<float>(std::sqrt(std::max(0.0, variance)))};
}


struct ComplexStats {
    float mean_real;
    float stddev_real;
    float mean_imag;
    float stddev_imag;
};

ComplexStats calculate_complex_stats(
    const fftwf_complex* data, 
    int n_rows, 
    int n_cols, 
    bool exclude_dc) {

    double sum_real = 0.0;
    double sum_sq_real = 0.0;
    double sum_imag = 0.0;
    double sum_sq_imag = 0.0;
    size_t count = 0;

    for (int r = 0; r < n_rows; ++r) {
        for (int c = 0; c < n_cols; ++c) {
            if (exclude_dc && r == 0 && c == 0) { // DC component is at (0,0) in unshifted FFT data
                continue;
            }
            int idx = r * n_cols + c;
            sum_real += data[idx][0];
            sum_sq_real += data[idx][0] * data[idx][0];
            sum_imag += data[idx][1];
            sum_sq_imag += data[idx][1] * data[idx][1];
            count++;
        }
    }

    ComplexStats result = {0.0f, 0.0f, 0.0f, 0.0f};
    if (count > 0) {
        result.mean_real = static_cast<float>(sum_real / count);
        double variance_real = (sum_sq_real / count) - (result.mean_real * result.mean_real);
        result.stddev_real = static_cast<float>(std::sqrt(std::max(0.0, variance_real)));

        result.mean_imag = static_cast<float>(sum_imag / count);
        double variance_imag = (sum_sq_imag / count) - (result.mean_imag * result.mean_imag);
        result.stddev_imag = static_cast<float>(std::sqrt(std::max(0.0, variance_imag)));
    }
    return result;
}

FftPeakParameters calculate_fft_peak_parameters(
    const std::vector<std::vector<float>>& fft_shifted_amplitude,
    const fftwf_complex* fft_out_complex_data, // FFT output (unshifted) for noise statistics
    float effective_integration_length, // For physical rate calculation
    int N_rows_padded, // Dimension of fft_shifted_amplitude and fft_out_complex_data
    int N_cols_padded_fft, // Dimension of fft_shifted_amplitude and fft_out_complex_data
    double delay_search_min_param,
    double delay_search_max_param,
    bool delay_search_range_specified_param,
    double rate_search_min_param,
    double rate_search_max_param,
    bool rate_search_range_specified_param) {

    FftPeakParameters params_result; // Renamed to avoid conflict if ProgramOptions params were passed directly

    if (fft_shifted_amplitude.empty() || fft_shifted_amplitude[0].empty()) {
        params_result.message = "Input amplitude data is empty.";
        params_result.success = false;
        return params_result;
    }

    // Check consistency of input dimensions with provided N_rows_padded and N_cols_padded_fft
    if (static_cast<int>(fft_shifted_amplitude.size()) != N_rows_padded || static_cast<int>(fft_shifted_amplitude[0].size()) != N_cols_padded_fft) {
        params_result.message = "Dimension mismatch between fft_shifted_amplitude and N_rows/N_cols_padded.";
        // Optionally, could try to use fft_shifted_amplitude.size() and [0].size() directly if that's intended.
        // For now, strict check.
        params_result.success = false;
        return params_result;
    }

    bool peak_found_in_range = false;
    for (int i = 0; i < N_rows_padded; ++i) {
        for (int j = 0; j < N_cols_padded_fft; ++j) {
            // Calculate physical rate and delay for current pixel (i, j)
            float current_physical_rate_hz = 0.0f;
            if (std::abs(effective_integration_length) > 1e-9f && N_rows_padded > 0) {
                float max_fringe_freq_val = 1.0f / (2.0f * effective_integration_length);
                current_physical_rate_hz = (static_cast<float>(i) - static_cast<float>(N_rows_padded) / 2.0f) *
                                           (2.0f * max_fringe_freq_val / static_cast<float>(N_rows_padded));
            }
            float current_physical_delay_samples = 0.0f;
            if (N_cols_padded_fft > 0) {
                // Changed sign convention for delay:
                // j < N/2 now means positive delay
                // j > N/2 now means negative delay
                // Example: N=1024, j=0   -> +512
                //          N=1024, j=511 -> +1
                //          N=1024, j=512 -> 0
                //          N=1024, j=1023-> -511
                current_physical_delay_samples = (static_cast<float>(N_cols_padded_fft) / 2.0f) - static_cast<float>(j);
            }

            // Check if current pixel is within specified search ranges
            bool in_delay_range = true;
            if (delay_search_range_specified_param) {
                if (current_physical_delay_samples < static_cast<float>(delay_search_min_param) || current_physical_delay_samples > static_cast<float>(delay_search_max_param)) {
                    in_delay_range = false;
                }
            }

            bool in_rate_range = true;
            if (rate_search_range_specified_param) {
                if (current_physical_rate_hz < static_cast<float>(rate_search_min_param) || current_physical_rate_hz > static_cast<float>(rate_search_max_param)) {
                    in_rate_range = false;
                }
            }

            if (in_delay_range && in_rate_range) {
                if (fft_shifted_amplitude[static_cast<size_t>(i)][static_cast<size_t>(j)] > params_result.max_amplitude) {
                    params_result.max_amplitude = fft_shifted_amplitude[static_cast<size_t>(i)][static_cast<size_t>(j)];
                    params_result.max_i = i;
                    params_result.max_j = j;
                    // Store the physical values of the peak found so far
                    params_result.physical_rate_hz = current_physical_rate_hz;
                    params_result.physical_delay_samples = current_physical_delay_samples;
                    peak_found_in_range = true;
                }
            }
        }
    }

    // Calculate mean of the (shifted) amplitude data for the SNR numerator
    float original_mean_amplitude, temp_stddev_not_used_here;
    std::tie(original_mean_amplitude, temp_stddev_not_used_here) = calculate_mean_stddev_from_float_vec(fft_shifted_amplitude);
    params_result.mean_amplitude = original_mean_amplitude; // Store the original mean for reference

     // --- SNRの分母となるノイズレベルの計算 (ユーザー指示による変更) ---
    // fft_out_complex_data (unshifted complex FFT data) を使用する
    double complex_sum_real = 0.0;
    double complex_sum_imag = 0.0;
    size_t complex_data_count = 0;
    const int total_complex_elements = N_rows_padded * N_cols_padded_fft;

    // 1. データの平均値を計算 (複素数) - DC成分を含めて全体の平均を計算
    // fft_out_complex_data は unshifted なので、DCは (0,0) にある。
    // Pythonの np.mean(input_2D_data) の挙動に合わせ、ここではDCを除外せずに全体の平均を計算する。
    if (fft_out_complex_data && total_complex_elements > 0) {
        for (int i = 0; i < total_complex_elements; ++i) {
            complex_sum_real += fft_out_complex_data[i][0];
            complex_sum_imag += fft_out_complex_data[i][1];
            complex_data_count++;
        }
    }

    float noise_level_for_snr = 0.0f;
    if (complex_data_count > 0) {
        std::complex<double> mean_complex(
            complex_sum_real / complex_data_count,
            complex_sum_imag / complex_data_count
        );

        // 2. 元のデータから平均値を差し引いたものの絶対値の総和を計算
        double sum_abs_diff = 0.0;
        for (int i = 0; i < total_complex_elements; ++i) {
            std::complex<double> val_complex(fft_out_complex_data[i][0], fft_out_complex_data[i][1]);
            sum_abs_diff += std::abs(val_complex - mean_complex);
        }
        // 3. ノイズレベルを計算 (上記2の平均)
        noise_level_for_snr = static_cast<float>(sum_abs_diff / complex_data_count);
    }
    // params_result.stddev_amplitude に計算したノイズレベルを格納する。
    // 注意: この変数は元々標準偏差を格納していたが、SNR計算方法の変更に伴いノイズレベルを格納する。
    params_result.stddev_amplitude = noise_level_for_snr;
    // --- ノイズレベルの計算ここまで ---


    if (!peak_found_in_range) {
        params_result.message = "No peak found within the specified search ranges.";
        params_result.success = false;
        // Reset peak dependent values if no peak in range
        params_result.max_amplitude = -1.0f; // Already default, but explicit
        params_result.max_i = -1;
        params_result.max_j = -1;
        params_result.physical_rate_hz = 0.0f;
        params_result.physical_delay_samples = 0.0f;
        params_result.snr = 0.0f;
        return params_result;
    }

    // Calculate SNR: max_amplitude / noise_level
    // params_result.max_amplitude was found from the original fft_shifted_amplitude.
    // params_result.stddev_amplitude には上記で計算した noise_level が格納されている。
    if (params_result.stddev_amplitude > 1e-9f) { // Ensure noise_level is not zero or too small
        params_result.snr = (params_result.max_amplitude) / params_result.stddev_amplitude;
    } else {
        // ノイズレベルがほぼ0の場合の処理
        if (std::abs(params_result.max_amplitude) < 1e-9f) { // 最大振幅も0に近い場合はSNR=0
            params_result.snr = 0.0f;
        } else {
            // ノイズがほぼ0で信号がある場合、SNRは非常に大きな値になる (無限大として扱うことも可能)
            params_result.snr = std::numeric_limits<float>::infinity(); 
        }
    }
    
    if (params_result.max_i != -1) { // Double check if a peak was actually recorded
        params_result.success = true;
        params_result.message = "Parameters calculated successfully (peak within specified range).";
    } else {
        params_result.message = "Internal error: peak_found_in_range was true, but max_i is -1.";
        params_result.success = false; // Should not happen if logic is correct
    }
    return params_result;
}
