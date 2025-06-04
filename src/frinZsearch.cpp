//
// Made by M.AKIMOTO with Gemini 
// 2025/06/01
// 


#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <fftw3.h> 
#include <chrono> // For time conversion
#include <ctime>   // For std::gmtime, std::strftime
#include <cmath> // For M_PI, std::sqrt, std::sin, std::cos, std::atan2, fmod
#include <string> // For std::string
#include <vector>
#include <filesystem> // For directory creation (C++17)
#include <limits>     // For std::numeric_limits
//#include <thread>     // For std::thread::hardware_concurrency
#include "frinZargs.hpp" // ★ この行が正しく存在することを確認してください
#include "frinZread.hpp" // Include the new header for file reading functions and HeaderRegion
#include "frinZlogger.hpp"    // Include the new Logger header
#include "frinZparamcal.hpp" // For parameter calculation
#include "frinZfftshift.hpp" // For fftshift and amplitude calculation
#include "frinZfitting.hpp" 


// Generates the 2D shifted amplitude plane for given correction parameters
std::vector<std::vector<float>> generate_shifted_amplitude_plane(
    const std::vector<std::vector<std::complex<float>>>& base_spectrums,
    float rate_hz_for_correction,
    float delay_samples_for_correction,
    fftwf_complex* fft_in_buffer, // Output: Populated with phase-corrected data for vertical FFT
    fftwf_complex* fft_intermediate_buffer, // Added this parameter
    fftwf_complex* fft_out_buffer,
    fftwf_plan plan_vertical,
    fftwf_plan plan_horizontal,
    int N_rows_original_data,
    int N_cols_original_data,
    int N_rows_padded_fft,
    int N_cols_padded_fft_dim,
    const std::vector<std::pair<int, int>>& rfi_indices,
    const FileData& file_data_ref,
    float effective_integ_length
) {

    std::vector<std::vector<std::complex<float>>> temp_spectra_for_trial = base_spectrums;

    // Apply phase correction (Simplified: full logic from main phase_calibration_routine needed here)
    // Initialize fft_in_buffer to zero. This is crucial as not all padded elements might be filled.
    for (int idx = 0; idx < N_rows_padded_fft * N_cols_padded_fft_dim; ++idx) {
        fft_in_buffer[idx][0] = 0.0f;
        fft_in_buffer[idx][1] = 0.0f;
    }

    // Prepare for phase correction
    const float sampling_freq_hz_iter = static_cast<float>(file_data_ref.header.sampling_speed);
    const uint32_t total_fft_points_iter = file_data_ref.header.fft_point;
     // N_cols_original_data should be equivalent to total_fft_points_iter / 2

    float py_equiv_sampling_speed_mhz_iter = 0.0f;
    float stop_val_for_linspace_mhz_iter = 0.0f; // Used for calculating frequency for delay correction
    bool can_phase_correct = (sampling_freq_hz_iter > 1e-9f && total_fft_points_iter >= 2 && std::abs(effective_integ_length) > 1e-9f);

    if (can_phase_correct) {
        py_equiv_sampling_speed_mhz_iter = sampling_freq_hz_iter / 1.0e6f;
        // This calculation is for the frequency axis (delay correction)
        stop_val_for_linspace_mhz_iter = static_cast<float>(static_cast<int>(py_equiv_sampling_speed_mhz_iter / 2.0f)) - 1.0f;
    }

    // Copy base_spectrums to fft_in_buffer, applying RFI masking and phase correction
    const int rows_to_process_for_fft = std::min(N_rows_original_data, N_rows_padded_fft);
    for (int r_orig = 0; r_orig < rows_to_process_for_fft; ++r_orig) {
        std::complex<float> rate_corr_factor_iter = {1.0f, 0.0f};
        if (can_phase_correct) {
            float time_for_rate_corr_sec_iter = static_cast<float>(r_orig) * effective_integ_length;
            rate_corr_factor_iter = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * rate_hz_for_correction * time_for_rate_corr_sec_iter));
        }
        for (int c_orig = 0; c_orig < N_cols_original_data; ++c_orig) {
            bool is_rfi_channel_iter = false;
            for (const auto& rfi_idx_range : rfi_indices) {
                if (c_orig >= rfi_idx_range.first && c_orig <= rfi_idx_range.second) { is_rfi_channel_iter = true; break; }
            }
            int buffer_idx = r_orig * N_cols_padded_fft_dim + c_orig;

            if (is_rfi_channel_iter) {
                fft_in_buffer[buffer_idx][0] = 0.0f;
                fft_in_buffer[buffer_idx][1] = 0.0f;
            } else {
                // Bounds check for base_spectrums (should be safe if N_rows_original_data/N_cols_original_data are correct)
                std::complex<float> original_val = base_spectrums[static_cast<size_t>(r_orig)][static_cast<size_t>(c_orig)];
                std::complex<float> corrected_val = original_val;

                if (can_phase_correct) {
                    // Delay correction factor
                    float freq_k_hz_for_delay_corr_iter = (N_cols_original_data > 1) ? 
                        (static_cast<float>(c_orig) * stop_val_for_linspace_mhz_iter / (static_cast<float>(N_cols_original_data) - 1.0f)) * 1.0e6f : 0.0f;
                    float delay_seconds_iter = delay_samples_for_correction / sampling_freq_hz_iter;
                    std::complex<float> delay_corr_factor_iter = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * delay_seconds_iter * freq_k_hz_for_delay_corr_iter));
                    
                    corrected_val *= rate_corr_factor_iter * delay_corr_factor_iter;
                }

                fft_in_buffer[buffer_idx][0] = corrected_val.real();
                fft_in_buffer[buffer_idx][1] = corrected_val.imag();
            }
            
            // Explicitly set DC component (frequency bin 0) to 0+0j
            // This is applied after RFI masking and phase correction for c_orig = 0.
            if (c_orig == 0) {
                fft_in_buffer[buffer_idx][0] = 0.0f;
                fft_in_buffer[buffer_idx][1] = 0.0f;
            }
        }   
    }

    // Execute vertical FFT (output to fft_intermediate_buffer)
    fftwf_execute(plan_vertical); 

    // Apply normalization to fft_intermediate_buffer (consistent with main loop's initial FFT)
    // N_rows_padded_fft is header.fft_point. Factor simplifies to 1.0f / N_cols_padded_fft_dim.
    const float norm_factor =  1; //static_cast<float>(file_data_ref.header.fft_point) / static_cast<float>(N_rows_padded_fft);
    if (N_rows_original_data > 0 && std::abs(effective_integ_length) > 1e-9f) { // Condition from main loop
        for (int r = 0; r < N_rows_padded_fft; ++r) {
            for (int c = 0; c < N_cols_padded_fft_dim; ++c) {
                int index = r * N_cols_padded_fft_dim + c;
                fft_intermediate_buffer[index][0] *= norm_factor;
                fft_intermediate_buffer[index][1] *= norm_factor;
            }
        }
    }



    // Execute horizontal FFT (input from fft_intermediate_buffer, output to fft_out_buffer)
    fftwf_execute(plan_horizontal); // plan_horizontal is now FFTW_BACKWARD
    // Apply scaling for IFFT (1/N)
    if (N_cols_padded_fft_dim > 0) {
        float scale_factor_horizontal_ifft = 0.0f;
        // N_rows_original_data is the number of actual data rows (sectors) used for this FFT
        if (N_rows_original_data > 0 && N_cols_padded_fft_dim > 0) { 
            scale_factor_horizontal_ifft = (1/ static_cast<float>(N_cols_padded_fft_dim)) ; 
        } else if (N_cols_padded_fft_dim > 0) { // Fallback if N_rows_original_data is 0 (should not happen if data is present)
             scale_factor_horizontal_ifft =  static_cast<float>(N_cols_padded_fft_dim);
             // No logger here as this function is more generic and logger is not passed
        }

        if (std::abs(scale_factor_horizontal_ifft) > 1e-12f) { // Apply only if factor is non-zero
            for (int r = 0; r < N_rows_padded_fft; ++r) {
                for (int c = 0; c < N_cols_padded_fft_dim; ++c) {
                    int index = r * N_cols_padded_fft_dim + c;
                    fft_out_buffer[index][0] *= scale_factor_horizontal_ifft;
                    fft_out_buffer[index][1] *= scale_factor_horizontal_ifft;
                }
            }
        }
    }

    std::vector<std::vector<float>> trial_shifted_amplitude = perform_fftshift_and_calc_amplitude(fft_out_buffer, N_rows_padded_fft, N_cols_padded_fft_dim);

    return trial_shifted_amplitude;
}

// Helper function for iterative FFT and peak finding (now uses generate_shifted_amplitude_plane)
FftPeakParameters run_iterative_fft_and_find_peak(
    const std::vector<std::vector<std::complex<float>>>& base_spectrums,
    float rate_hz,
    float delay_samples,
    fftwf_complex* fft_in_buffer,
    fftwf_complex* fft_intermediate_buffer,
    fftwf_complex* fft_out_buffer,
    fftwf_plan plan_vertical,
    fftwf_plan plan_horizontal,
    int N_rows_original_data,
    int N_cols_original_data,
    int N_rows_padded_fft,
    int N_cols_padded_fft_dim,
    const std::vector<std::pair<int, int>>& rfi_indices,
    const FileData& file_data_ref, 
    float effective_integ_length,
    const ProgramOptions& prog_params) { // Added ProgramOptions to pass search range

    // Note: fft_out_buffer from generate_shifted_amplitude_plane is the unshifted complex data
    std::vector<std::vector<float>> trial_shifted_amplitude = generate_shifted_amplitude_plane(
        base_spectrums, rate_hz, delay_samples, fft_in_buffer, fft_intermediate_buffer, fft_out_buffer,
        plan_vertical, plan_horizontal, N_rows_original_data, N_cols_original_data,
        N_rows_padded_fft, N_cols_padded_fft_dim, rfi_indices, file_data_ref, effective_integ_length);
        
    FftPeakParameters trial_params;
    if (!trial_shifted_amplitude.empty()) {
        trial_params = calculate_fft_peak_parameters(
            trial_shifted_amplitude,
            fft_out_buffer, // Pass the complex FFT output for noise statistics 
            effective_integ_length, 
            N_rows_padded_fft, // This corresponds to N_rows_padded in calculate_fft_peak_parameters
            N_cols_padded_fft_dim, // This corresponds to N_cols_padded_fft in calculate_fft_peak_parameters
            prog_params.delay_search_min,
            prog_params.delay_search_max,
            prog_params.delay_search_range_specified,
            prog_params.rate_search_min,
            prog_params.rate_search_max,
            prog_params.rate_search_range_specified);
    } else {
        trial_params.success = false;
        trial_params.message = "Shifted amplitude data empty in iterative FFT helper.";
    }
    return trial_params;
}

// Helper function to write 2D float vector to a CSV file
void write_2d_float_vector_to_csv(const std::vector<std::vector<float>>& data, const std::string& filename, Logger& logger, bool noconsole) {
    if (data.empty()) {
        if (!noconsole) logger << "Data to write to " << filename << " is empty. Skipping." << std::endl;
        return;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        if (!noconsole) logger << "Error: Could not open file for writing: " << filename << std::endl;
        return;
    }

    if (!noconsole) logger << "Writing 2D data to " << filename << "..." << std::endl;
    outfile << std::fixed << std::setprecision(8); 

    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = 0; j < data[i].size(); ++j) {
            outfile << data[i][j] << (j < data[i].size() - 1 ? "," : "");
        }
        outfile << "\n";
    }
    outfile.close();
    if (!noconsole) logger << "Finished writing to " << filename << "." << std::endl;
}

int main(int argc, char* argv[]) {

    ProgramOptions params;
    Logger logger;
    
    // Parse arguments using the new function from frinZargs.cpp
    // This function handles --help, --version, and basic parsing errors.
    if (!parse_arguments(argc, argv, params)) {
        // parse_arguments will print help/version or error messages.
        // Return 1 indicates an error or that help/version was displayed.
        return 1;
    }

    // Post-process options (e.g., create output directory, set dependent options)
    post_process_options(params);

    /*
    // Initialize FFTW threads support
    if (!fftwf_init_threads()) {
        std::cerr << "Error: Failed to initialize FFTW threads support." << std::endl;
        // It might still be possible to run in single-threaded mode,
        // or you might choose to exit. For now, we'll just print an error.
    } else {
        // Determine the number of threads to use (e.g., number of hardware cores)
        int nthreads = std::thread::hardware_concurrency()/2;
        fftwf_plan_with_nthreads(nthreads > 0 ? nthreads : 1); // Use at least 1 thread
    }
        */

    // Setup logger (output_dir_final is now set by finalize_options)
    std::string text_log_file_full_path_str;
    if (params.enable_text_log_output) {
        // params.output_dir_final should be set correctly by finalize_options
        // If it's empty (e.g., due to an error during its creation), this block will be skipped.
        if (!params.output_dir_final.empty()) {
            namespace fs = std::filesystem;
            fs::path input_file_path_fs(params.input_filename); // input_filename is guaranteed by --input->required()
            std::string input_basename = input_file_path_fs.stem().string(); // Add this line
            text_log_file_full_path_str = (fs::path(params.output_dir_final) / (input_basename + "_frinZsearch_result.txt")).string();
        }
    }

    logger.setup(!params.noconsole, params.enable_text_log_output, text_log_file_full_path_str);

    if (params.enable_text_log_output && !params.output_dir_final.empty() && !params.noconsole) {
        logger << "Output directory for all files: " << params.output_dir_final << std::endl;
    }

    // === Step 1: ヘッダーとデータ読み込み (frinZread.cpp を使用) ===
    FileData loaded_data = read_binary_file_data(params.input_filename, logger, params.output_header_info);
    if (!loaded_data.success) {
        // Error message already printed by read_binary_file_data or logger
        return 1;
    }

    // Access loaded data
    const HeaderRegion& header = loaded_data.header; // Use const& for read-only access
    const std::vector<std::vector<std::complex<float>>>& raw_all_spectrums = loaded_data.spectrum_data;
    const std::vector<std::string>& raw_all_sector_start_times_utc = loaded_data.sector_start_times_utc;
    float first_effective_integration_length = loaded_data.first_effective_integration_length;

    // fft_point, PP, fft_half are now derived from loaded_data.header
    uint32_t fft_point = header.fft_point;
    // uint32_t PP = header.num_sector; // PP is not directly used later, raw_all_spectrums.size() is used
    // size_t fft_half = fft_point / 2; // fft_half is also not directly used later, N_cols_original is used

    if (raw_all_spectrums.empty()) { // This check is still valid
        std::cerr << "Error: No data read from file.\n";
        return 1;
    }
    if (loaded_data.first_effective_integration_length <= 1e-9f && (params.skip_seconds > 0 || params.specified_length_sec > 0)) {
        std::cerr << "Warning: Effective integration length is zero or too small. Cannot reliably use --skip or --length." << std::endl;
        // Proceed with caution, or exit if these options are critical
    }

    // Apply --skip
    std::vector<std::vector<std::complex<float>>> data_after_skip;
    // std::vector<std::string> times_after_skip; // If needed for logging per segment

    uint32_t num_sectors_to_skip = 0;
    if (params.skip_seconds > 0 && first_effective_integration_length > 1e-9) {
        num_sectors_to_skip = static_cast<uint32_t>(std::round(params.skip_seconds / loaded_data.first_effective_integration_length));
    }

    if (num_sectors_to_skip > 0 && num_sectors_to_skip < raw_all_spectrums.size()) {
        data_after_skip.assign(raw_all_spectrums.begin() + num_sectors_to_skip, raw_all_spectrums.end());
        // times_after_skip.assign(raw_all_sector_start_times_utc.begin() + num_sectors_to_skip, raw_all_sector_start_times_utc.end());
        if (!params.noconsole) {
            logger << "Skipped " << num_sectors_to_skip << " initial sectors (" << params.skip_seconds << "s)." << std::endl;
        }
    } else if (num_sectors_to_skip >= raw_all_spectrums.size() && !raw_all_spectrums.empty()) {
        std::cerr << "Error: --skip value (" << params.skip_seconds << "s) results in skipping all data (" << num_sectors_to_skip << " sectors)." << std::endl;
        return 1;
    } else {
        data_after_skip = raw_all_spectrums;
        // times_after_skip = raw_all_sector_start_times_utc;
    }

    if (data_after_skip.empty()) {
        std::cerr << "Error: No data remaining after applying --skip.\n";
        return 1;
    }

    // Hardcoded padding factors and FFTW plan flag
    //const double pad_factor_rows_val = 2.0;
    const double pad_factor_cols_val = 2.0; // No Edit this line
    const unsigned fftw_plan_flag_val = FFTW_ESTIMATE;

    uint32_t total_data_after_skip_size = data_after_skip.size();

    for (int loop_idx = 0; loop_idx < params.loop_count; ++loop_idx) {
        
        uint32_t current_segment_num_sectors;
        uint32_t segment_start_index_in_data_to_process;

        if (params.specified_length_sec > 0 && first_effective_integration_length > 1e-9) {
            current_segment_num_sectors = static_cast<uint32_t>(std::round(params.specified_length_sec / first_effective_integration_length));
            if (current_segment_num_sectors == 0 && params.specified_length_sec > 0) current_segment_num_sectors = 1; // Process at least one if length specified
            segment_start_index_in_data_to_process = loop_idx * current_segment_num_sectors;
        } else { // No length specified, or cannot calculate from length
            uint32_t base_segment_size = (total_data_after_skip_size + params.loop_count - 1) / params.loop_count; // sectors per loop if dividing evenly
            current_segment_num_sectors = base_segment_size;
            segment_start_index_in_data_to_process = loop_idx * base_segment_size;
        }
        

        if (segment_start_index_in_data_to_process >= total_data_after_skip_size) {
            if (!params.noconsole) logger << "No more data for loop " << loop_idx + 1 << ". Ending." << std::endl;
            break; 
        }

        uint32_t actual_segment_end_index = segment_start_index_in_data_to_process + current_segment_num_sectors;
        if (actual_segment_end_index > total_data_after_skip_size) {
            actual_segment_end_index = total_data_after_skip_size;
        }

        if (!params.noconsole) {
            logger << "\n--- Processing Loop " << loop_idx + 1 << "/" << params.loop_count << " ---" << std::endl;
            logger << "Input file: " << params.input_filename << std::endl;
            logger << "Target: " << std::string(header.source_name, 16).c_str() << std::endl;
            if (segment_start_index_in_data_to_process < raw_all_sector_start_times_utc.size()) {
                logger << "Segment Start UTC: " << raw_all_sector_start_times_utc[num_sectors_to_skip + segment_start_index_in_data_to_process] << std::endl;
            }
            if (actual_segment_end_index > 0 && (num_sectors_to_skip + actual_segment_end_index -1) < raw_all_sector_start_times_utc.size()) {
                 logger << "Segment End UTC (approx): " << raw_all_sector_start_times_utc[num_sectors_to_skip + actual_segment_end_index - 1] << std::endl;
            }
             logger << "Processing " << current_segment_num_sectors << " sectors for this loop." << std::endl;
        }
        
        current_segment_num_sectors = actual_segment_end_index - segment_start_index_in_data_to_process;

        if (current_segment_num_sectors == 0) {
            if (!params.noconsole) logger << "Segment for loop " << loop_idx + 1 << " is empty. Skipping." << std::endl;
            continue;
        }

        std::vector<std::vector<std::complex<float>>> current_all_spectrums(
            data_after_skip.begin() + segment_start_index_in_data_to_process,
            data_after_skip.begin() + actual_segment_end_index
        );

        if (current_all_spectrums.empty() || current_all_spectrums[0].empty()) {
            if (!params.noconsole) std::cerr << "Error: Segment for loop " << loop_idx + 1 << " is empty or invalid.\n";
            continue;
        }

        // === 初期位相較正 (コマンドラインオプション --delay, --rate が指定された場合) ===
        // current_all_spectrums に対して直接適用し、最初のFFTの入力とする
        
        if (params.initial_delay_specified || params.initial_rate_specified) {
            if (!params.noconsole) {
                logger << "Applying initial phase correction based on --delay and/or --rate flags before first FFT." << std::endl;
                if (params.initial_delay_specified) logger << "  Initial Delay for correction: " << params.initial_delay_samples << " samples." << std::endl;
                if (params.initial_rate_specified) logger << "  Initial Rate for correction: " << params.initial_rate_hz << " Hz." << std::endl;
            }

            const float rate_hz_for_initial_correction = params.initial_rate_specified ? params.initial_rate_hz : 0.0f;
            const float delay_samples_for_initial_correction = params.initial_delay_specified ? params.initial_delay_samples : 0.0f;

            const float sampling_freq_hz_corr = static_cast<float>(loaded_data.header.sampling_speed);
            const uint32_t total_fft_points_corr = loaded_data.header.fft_point;
            // N_cols_original は current_all_spectrums[0].size() から取得するべきだが、
            // この時点ではまだ N_cols_original 変数は定義されていない。
            // ただし、current_all_spectrums[0].size() は total_fft_points_corr / 2 と等価のはず。
            const int N_cols_original_for_corr = static_cast<int>(current_all_spectrums[0].size());

            bool can_initial_phase_correct = (sampling_freq_hz_corr > 1e-9f && total_fft_points_corr >= 2 && std::abs(first_effective_integration_length) > 1e-9f && N_cols_original_for_corr > 0);

            if (can_initial_phase_correct) {
                float py_equiv_sampling_speed_mhz_corr = sampling_freq_hz_corr / 1.0e6f;
                float stop_val_for_linspace_mhz_corr = static_cast<float>(static_cast<int>(py_equiv_sampling_speed_mhz_corr / 2.0f)) - 1.0f;

                for (size_t s_idx = 0; s_idx < current_all_spectrums.size(); ++s_idx) {
                    float time_for_rate_corr_sec = static_cast<float>(s_idx) * first_effective_integration_length;
                    std::complex<float> rate_corr_factor = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * rate_hz_for_initial_correction * time_for_rate_corr_sec));

                    for (int k_idx = 0; k_idx < N_cols_original_for_corr; ++k_idx) {
                        float freq_k_hz_for_delay_corr = (N_cols_original_for_corr > 1) ?
                            (static_cast<float>(k_idx) * stop_val_for_linspace_mhz_corr / (static_cast<float>(N_cols_original_for_corr) - 1.0f)) * 1.0e6f : 0.0f;
                        float delay_seconds_corr = delay_samples_for_initial_correction / sampling_freq_hz_corr;
                        std::complex<float> delay_corr_factor = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * delay_seconds_corr * freq_k_hz_for_delay_corr));
                        current_all_spectrums[s_idx][static_cast<size_t>(k_idx)] *= rate_corr_factor * delay_corr_factor;
                    }
                }
                 if (!params.noconsole) logger << "Initial phase correction applied to current segment data." << std::endl;
            } else {
                if (!params.noconsole) logger << "Skipping initial phase correction due to insufficient data (sampling_freq, fft_points, integ_length, or N_cols_original_for_corr=0)." << std::endl;
            }
        }
        

        int N_rows_original = static_cast<int>(current_all_spectrums.size());
        int N_cols_original = static_cast<int>(current_all_spectrums[0].size()); // Assuming consistent inner size


        // === ゼロパディング計算 (縦軸) === N_rows_original (current_segment_num_sectors) を使用
        // Rate direction (vertical axis) FFT size is set to the smallest power of 2 >= N_rows_original.
        int N_rows_padded;
        if (N_rows_original > 0) {
            double temp_rows = static_cast<double>(N_rows_original);
            int n_power_rows = 0;
            while (std::pow(2, n_power_rows) < temp_rows) {
                n_power_rows++;
            }
            N_rows_padded = static_cast<int>(std::pow(2, n_power_rows)) * pad_factor_cols_val;
        } else {
            N_rows_padded = 1; // Avoid 0 if N_rows_original is 0, FFTW can handle size 1.
        }

        // === ゼロパディング計算 (横軸) === header.fft_point を使用
        // Delay direction (horizontal axis) FFT size is set to header.fft_point.
        // N_cols_original is header.fft_point / 2.
        // The raw data (N_cols_original) will be zero-padded to N_cols_fft.
        int N_cols_padded_fft = static_cast<int>(fft_point); // fft_point is loaded_data.header.fft_point
        int N_cols_fft = N_cols_padded_fft; // N_cols_padded_fft is used in some places as an alias or for clarity, ensure it's consistent
        if (!params.noconsole) {
            logger << "Original Data Dimensions: Vertical (Sectors) = " << N_rows_original
                   << ", Horizontal (Channels) = " << N_cols_original << std::endl;
            logger << "Padded FFT Dimensions: Vertical (N_rows_padded) = " << N_rows_padded
                   << ", Horizontal (N_cols_fft) = " << N_cols_fft << std::endl;
        }


        fftwf_complex *fft_in, *fft_intermediate, *fft_out;
        fftwf_plan plan_vertical, plan_horizontal;

        fft_in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
        fft_intermediate = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
        fft_out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);

        if (!fft_in || !fft_intermediate || !fft_out) {
            std::cerr << "Error: Failed to allocate memory for FFT (Loop " << loop_idx + 1 << ")." << std::endl;
            if(fft_in) fftwf_free(fft_in);
            if(fft_intermediate) fftwf_free(fft_intermediate);
            if(fft_out) fftwf_free(fft_out);
            continue; // Skip to next loop iteration or break
        }

        // fftwf_complex* 配列をゼロで初期化
        for (int i = 0; i < N_rows_padded * N_cols_fft; ++i) {
            fft_in[i][0] = 0.0f;
            fft_in[i][1] = 0.0f;
        }

        // Calculate RFI index ranges once
        std::vector<std::pair<int, int>> rfi_index_ranges;
        if (!params.rfi_ranges_mhz.empty() && loaded_data.header.sampling_speed > 0 && loaded_data.header.fft_point > 0) {
            if (!params.noconsole) logger << "Calculating RFI cut indices..." << std::endl;
            double total_bandwidth_hz = static_cast<double>(loaded_data.header.sampling_speed);
            double nyquist_bandwidth_hz = total_bandwidth_hz / 2.0; // Nyquist is sampling_speed / 2
            uint32_t num_channels_half_fft = loaded_data.header.fft_point / 2; // Number of channels up to Nyquist
            double channel_resolution_hz = nyquist_bandwidth_hz / num_channels_half_fft;

            if (channel_resolution_hz > 1e-9) {
                for (const auto& range_mhz : params.rfi_ranges_mhz) { // range_mhz is std::pair<double, double> from params
                    double start_freq_hz = range_mhz.first * 1.0e6;  // Convert MHz to Hz
                    double end_freq_hz = range_mhz.second * 1.0e6; // Convert MHz to Hz

                    int start_idx = static_cast<int>(std::floor(start_freq_hz / channel_resolution_hz));
                    int end_idx = static_cast<int>(std::ceil(end_freq_hz / channel_resolution_hz)) -1;

                    start_idx = std::max(0, start_idx);
                    end_idx = std::min(static_cast<int>(num_channels_half_fft) - 1, end_idx);

                    if (start_idx <= end_idx) {
                        rfi_index_ranges.emplace_back(start_idx, end_idx);
                        if (!params.noconsole) logger << "  RFI Cut: " << range_mhz.first << "-" << range_mhz.second << " MHz maps to Indices: " << start_idx << "-" << end_idx << std::endl;
                    } else {
                        if (!params.noconsole) logger << "  RFI Cut: " << range_mhz.first << "-" << range_mhz.second << " MHz results in an invalid or empty index range after calculation (start_idx=" << start_idx << ", end_idx=" << end_idx <<"). Skipping this range." << std::endl;
                    }
                }
            } else {
                if (!params.noconsole) std::cerr << "Warning: Channel resolution is too small or zero. Cannot calculate RFI cut indices." << std::endl;
                }
        }
    
    // 元のデータを fft_in の左上にコピー (ゼロパディング)
    const int rows_to_copy_to_fft_in = std::min(N_rows_original, N_rows_padded);
    // N_cols_original は header.fft_point / 2。新しい N_cols_fft は header.fft_point。
    // よって、横方向には N_cols_original 列分コピーし、残りはゼロパディングされる。
    const int cols_to_copy_to_fft_in = N_cols_original; 
    for (int i = 0; i < rows_to_copy_to_fft_in; ++i) { // コピーする行数を N_rows_padded 以下に制限
        for (int j = 0; j < cols_to_copy_to_fft_in; ++j) { 
            bool is_rfi_channel = false;
            for (const auto& rfi_idx_range : rfi_index_ranges) {
                if (static_cast<int>(j) >= rfi_idx_range.first && static_cast<int>(j) <= rfi_idx_range.second) {
                    is_rfi_channel = true;
                    break;
                }
            }

            // fft_in のインデックスは i * N_cols_fft + j
            if (is_rfi_channel) {
                fft_in[i * N_cols_fft + j][0] = 0.0f;
                fft_in[i * N_cols_fft + j][1] = 0.0f;
            } else {
                fft_in[i * N_cols_fft + j][0] = current_all_spectrums[static_cast<size_t>(i)][static_cast<size_t>(j)].real();
                fft_in[i * N_cols_fft + j][1] = current_all_spectrums[static_cast<size_t>(i)][static_cast<size_t>(j)].imag();
            }
        }
    }

    // === Step 2a: 縦方向 (列ごと) の1D FFT ===
    int rank_vert = 1;
    int n_vert[] = {N_rows_padded};
    int howmany_vert = N_cols_fft;
    int istride_vert = N_cols_fft; // 縦方向に進むには N_cols_fft だけインデックスを増やす
    int idist_vert = 1;          // 次の列の先頭へは 1 だけインデックスを増やす
    int ostride_vert = N_cols_fft;
    int odist_vert = 1;

    plan_vertical = fftwf_plan_many_dft(rank_vert, n_vert, howmany_vert,
                                        fft_in, nullptr, istride_vert, idist_vert,
                                        fft_intermediate, nullptr, ostride_vert, odist_vert,
                                        FFTW_FORWARD, fftw_plan_flag_val);

    fftwf_execute(plan_vertical);
    //fftwf_destroy_plan(plan_vertical);

    // === Normalization after vertical FFT ===
    // Original: norm_factor = header.fft_point / (N_rows_padded * N_cols_fft)
    // Since new N_cols_fft is header.fft_point, this simplifies to:
    // norm_factor = N_cols_fft / (N_rows_padded * N_cols_fft) = 1.0f / N_rows_padded
    const float norm_factor = 1; //static_cast<float>(loaded_data.header.fft_point) / static_cast<float>(N_rows_padded * N_cols_fft);
    if (current_segment_num_sectors > 0 && std::abs(first_effective_integration_length) > 1e-9f) {
        if (!params.noconsole) {
            logger << "Applying normalization factor (" << norm_factor << ") after vertical FFT." << std::endl;
        }
        for (int r = 0; r < N_rows_padded; ++r) {
            for (int c = 0; c < N_cols_fft; ++c) { // N_cols_fft is N_cols_padded_fft
                int index = r * N_cols_fft + c;
                fft_intermediate[index][0] *= norm_factor;
                fft_intermediate[index][1] *= norm_factor;
            }
        }
    } else {
        if (!params.noconsole) {
            logger << "Skipping normalization after vertical FFT due to zero sectors or zero integration length." << std::endl;
        }
    }
    
    // === Step 2c: 横方向 (行ごと) の1D FFT ===
    int rank_horiz = 1;
    int n_horiz[] = {N_cols_fft};
    int howmany_horiz = N_rows_padded;
    int istride_horiz = 1;           // 横方向に進むには 1 だけインデックスを増やす
    int idist_horiz = N_cols_fft;  // 次の行の先頭へは N_cols_fft だけインデックスを増やす
    int ostride_horiz = 1;
    int odist_horiz = N_cols_fft;

    plan_horizontal = fftwf_plan_many_dft(rank_horiz, n_horiz, howmany_horiz,
                                          fft_intermediate, nullptr, istride_horiz, idist_horiz,
                                          fft_out, nullptr, ostride_horiz, odist_horiz,
                                          FFTW_BACKWARD, fftw_plan_flag_val); 

    fftwf_execute(plan_horizontal); // plan_horizontal is now FFTW_BACKWARD
    // Apply scaling for IFFT (1/N) after horizontal FFT
    // The existing scaling factor structure is N_cols_fft / N_rows_original (current_segment_num_sectors)
    // This is unusual for a standard IFFT (which would be 1/N_cols_fft). We'll maintain the structure.
    if (N_cols_fft > 0) {
        float scale_factor_horizontal_ifft = 0.0f;
        if (current_segment_num_sectors > 0 && N_cols_fft > 0) { // Avoid division by zero
            scale_factor_horizontal_ifft = (1 / static_cast<float>(N_cols_fft));
        } else if (N_cols_fft > 0) { // Fallback if sectors is 0, just scale by N_cols_fft
            scale_factor_horizontal_ifft = 1.0f; // Or perhaps 1.0f / N_cols_fft if N_rows_original is 0
            if (!params.noconsole) logger << "Warning: current_segment_num_sectors is 0 for horizontal IFFT scaling. Using 1/N_cols_fft only." << std::endl;
        }

        if (std::abs(scale_factor_horizontal_ifft) > 1e-12f) { // Apply only if factor is non-zero
            for (int r = 0; r < N_rows_padded; ++r) {
                for (int c = 0; c < N_cols_fft; ++c) {
                    int index = r * N_cols_fft + c;
                    fft_out[index][0] *= scale_factor_horizontal_ifft;
                    fft_out[index][1] *= scale_factor_horizontal_ifft;
                }
            }
        }
    }

    // === Step 3: fftshiftと振幅計算 ===
    std::vector<std::vector<float>> fft_shifted_amplitude = perform_fftshift_and_calc_amplitude(fft_out, N_rows_padded, N_cols_padded_fft);
    
    // Output initial FFT shifted amplitude for heatmap
     namespace fs = std::filesystem; // Declare fs alias here for wider scope
    if (params.enable_text_log_output && !params.output_dir_final.empty()) {
        std::string initial_fft_heatmap_filename = (fs::path(params.output_dir_final) / (fs::path(params.input_filename).stem().string() + "_loop" + std::to_string(loop_idx+1) + "_initial_fft_heatmap.csv")).string();
        write_2d_float_vector_to_csv(fft_shifted_amplitude, initial_fft_heatmap_filename, logger, params.noconsole);
    }


    // === Step 3.5: FFT後の振幅の最大値とその座標を検索 ===
    FftPeakParameters initial_peak_params;
    if (!fft_shifted_amplitude.empty()) {
        initial_peak_params = calculate_fft_peak_parameters(
            fft_shifted_amplitude, // Shifted amplitude for peak finding
            fft_out,               // Unshifted complex data for noise statistics
            first_effective_integration_length,
            N_rows_padded,
            N_cols_padded_fft,
            params.delay_search_min,
            params.delay_search_max,
            params.delay_search_range_specified,
            params.rate_search_min,
            params.rate_search_max,
            params.rate_search_range_specified
        );
    } else {
        initial_peak_params.success = false;
        initial_peak_params.message = "Initial FFT shifted amplitude data is empty.";
    }

    float refined_delay_from_first_fit_samples = initial_peak_params.physical_delay_samples;
    float refined_rate_from_first_fit_hz = initial_peak_params.physical_rate_hz;

    
    
    if (initial_peak_params.success) {
        // Declare refined parameters here to be accessible by Step 4
        // Initialize with values from initial peak parameters as a default.
        if (!params.noconsole) {
            logger << "Initial FFT - StdDev: " << initial_peak_params.stddev_amplitude 
                   << ", SNR (Max/StdDev): " << initial_peak_params.snr << std::endl;

            if (!params.noconsole) {
                logger << std::fixed;
                logger.setprecision(8) << "Maximum FFT amplitude (%): " << 100* initial_peak_params.max_amplitude
                                       << " (within specified range if any) found at" << std::endl
                                       << "  Fringe Frequency (Hz): " << initial_peak_params.physical_rate_hz
                                       << " (Range: " << (params.rate_search_range_specified ? "[" + std::to_string(params.rate_search_min) + "," + std::to_string(params.rate_search_max) + "]" : "unrestricted") << ")" << std::endl
                                       << "  Delay (Sample): " << initial_peak_params.physical_delay_samples 
                                       << " (Range: " << (params.delay_search_range_specified ? "[" + std::to_string(params.delay_search_min) + "," + std::to_string(params.delay_search_max) + "]" : "unrestricted") << ")" << std::endl;
                logger << std::dec; // Reset to decimal mode separately
            }
            
            // Note: refined_..._from_first_fit_samples/hz might have been clipped above by "Clip 1".
            // The quadratic fitting below uses the original initial_peak_params.max_i/max_j (unclipped indices).
            // Output surrounding points for delay and rate
            if (!params.noconsole && !fft_shifted_amplitude.empty()) {
                logger << "  Surrounding points (Initial FFT Peak at i=" << initial_peak_params.max_i << ", j=" << initial_peak_params.max_j << "):" << std::endl;
                
                // Delay dimension (fixed rate at max_i)
                logger << "    Delay (fixed rate bin " << initial_peak_params.max_i << "):" << std::endl;
                std::vector<double> x_coords_delay_fit;
                std::vector<double> y_values_delay_fit;
                bool all_points_valid_for_delay_fit = true;
                int delay_half_steps = (params.delay_fit_points - 1) / 2;
                for (int dj = -delay_half_steps; dj <= delay_half_steps; ++dj) {
                    int current_j = initial_peak_params.max_j + dj;
                     // Ensure max_i is within bounds of fft_shifted_amplitude's first dimension
                    if (current_j >= 0 && current_j < N_cols_padded_fft &&
                        static_cast<size_t>(initial_peak_params.max_i) < fft_shifted_amplitude.size() &&
                        static_cast<size_t>(current_j) < fft_shifted_amplitude[static_cast<size_t>(initial_peak_params.max_i)].size()) {
                        // Consistent with updated sign convention in frinZparamcal.cpp
                        float physical_delay = (static_cast<float>(N_cols_padded_fft) / 2.0f) - static_cast<float>(current_j);
                        // Ensure peak physical delay matches this calculation if current_j is initial_peak_params.max_j
                        // This is for display of surrounding points. The main peak's physical_delay comes from calculate_fft_peak_parameters.
                        float amplitude = fft_shifted_amplitude[static_cast<size_t>(initial_peak_params.max_i)][static_cast<size_t>(current_j)];
                        logger << "      Delay[j=" << current_j << "] (Phys: " << physical_delay << " samp): Amp = "
                               << amplitude << std::endl;
                        x_coords_delay_fit.push_back(static_cast<double>(physical_delay));
                        y_values_delay_fit.push_back(static_cast<double>(amplitude));
                    } else {
                        logger << "      Delay[j=" << current_j << "]: Out of bounds for fitting" << std::endl;
                        all_points_valid_for_delay_fit = false;
                    }
                    }

                if (all_points_valid_for_delay_fit && x_coords_delay_fit.size() == static_cast<size_t>(params.delay_fit_points)) {
                    QuadraticFitResult delay_fit_result = fit_quadratic_least_squares(x_coords_delay_fit, y_values_delay_fit);
                    if (delay_fit_result.success) {
                        logger << "      Quadratic Fit (Delay): Refined Peak Delay (samp) = " << delay_fit_result.peak_x
                                << ", Max Amp Est (%): " << 100* (delay_fit_result.a * delay_fit_result.peak_x * delay_fit_result.peak_x + delay_fit_result.b * delay_fit_result.peak_x + delay_fit_result.c)
                               << " (Coeffs a=" << delay_fit_result.a << ", b=" << delay_fit_result.b << ", c=" << delay_fit_result.c << ")" << std::endl;
                        refined_delay_from_first_fit_samples = static_cast<float>(delay_fit_result.peak_x);
                    } else {
                        logger << "      Quadratic Fit (Delay) failed: " << delay_fit_result.message << std::endl;
                    }
                } else {
                    logger << "      Quadratic Fit (Delay): Not enough valid/contiguous data points (" << x_coords_delay_fit.size() << " collected) for " << params.delay_fit_points << "-point fitting." << std::endl;
                }

                // Rate dimension (fixed delay at max_j)
                logger << "    Rate (fixed delay bin " << initial_peak_params.max_j << "):" << std::endl;
                std::vector<double> x_coords_rate_fit;
                std::vector<double> y_values_rate_fit;
                bool all_points_valid_for_rate_fit = true;
                
                int rate_half_steps = (params.rate_fit_points - 1) / 2;
                for (int di = -rate_half_steps; di <= rate_half_steps; ++di) {
                    int current_i = initial_peak_params.max_i + di;
                     // Ensure max_j is within bounds of fft_shifted_amplitude's second dimension
                    if (current_i >= 0 && current_i < N_rows_padded &&
                        static_cast<size_t>(current_i) < fft_shifted_amplitude.size() &&
                        static_cast<size_t>(initial_peak_params.max_j) < fft_shifted_amplitude[static_cast<size_t>(current_i)].size()) {
                         float physical_rate = (static_cast<float>(current_i) - static_cast<float>(N_rows_padded) / 2.0f) *
                                              ( (std::abs(first_effective_integration_length) > 1e-9f && N_rows_padded > 0) ? (2.0f * (1.0f / (2.0f * first_effective_integration_length)) / static_cast<float>(N_rows_padded)) : 0.0f );
                        float amplitude = fft_shifted_amplitude[static_cast<size_t>(current_i)][static_cast<size_t>(initial_peak_params.max_j)];
                        logger << "      Rate[i=" << current_i << "] (Phys: " << physical_rate << " Hz): Amp = "
                               << amplitude << std::endl;
                        x_coords_rate_fit.push_back(static_cast<double>(physical_rate));
                        y_values_rate_fit.push_back(static_cast<double>(amplitude));
                    } else {
                        logger << "      Rate[i=" << current_i << "]: Out of bounds for fitting" << std::endl;
                        all_points_valid_for_rate_fit = false;
                    }
                }
                if (all_points_valid_for_rate_fit && x_coords_rate_fit.size() == static_cast<size_t>(params.rate_fit_points)) {
                    // --- Rateフィッティングのためのスケーリング ---
                    std::vector<double> scaled_x_coords_rate_fit = x_coords_rate_fit;
                    const double rate_scale_factor = 100.0;
                    for (double& x_val : scaled_x_coords_rate_fit) {
                        x_val *= rate_scale_factor;
                    }

                    QuadraticFitResult rate_fit_result = fit_quadratic_least_squares(scaled_x_coords_rate_fit, y_values_rate_fit);
                    
                    if (rate_fit_result.success) {

                        // peak_x と係数を元のスケールに戻す
                        rate_fit_result.peak_x /= rate_scale_factor;
                        // y = a_fit * (S*x)^2 + b_fit * (S*x) + c_fit
                        // y = (a_fit*S^2) * x^2 + (b_fit*S) * x + c_fit
                        // なので、元のスケールの係数は a_orig = a_fit*S^2, b_orig = b_fit*S
                        rate_fit_result.a *= (rate_scale_factor * rate_scale_factor);
                        rate_fit_result.b *= rate_scale_factor;
                        // c_fit は不変

                        logger << "      Quadratic Fit (Rate): Refined Peak Rate (Hz) = " << rate_fit_result.peak_x 
                               << ", Max Amp Est (%): " << 100* (rate_fit_result.a * rate_fit_result.peak_x * rate_fit_result.peak_x + rate_fit_result.b * rate_fit_result.peak_x + rate_fit_result.c)
                               << " (Coeffs a=" << rate_fit_result.a << ", b=" << rate_fit_result.b << ", c=" << rate_fit_result.c << ")" << std::endl;
                               refined_rate_from_first_fit_hz = static_cast<float>(rate_fit_result.peak_x);
                    } else {
                        logger << "      Quadratic Fit (Rate) failed: " << rate_fit_result.message << std::endl;
                    }
                } else {
                    logger << "      Quadratic Fit (Rate): Not enough valid/contiguous data points (" << x_coords_rate_fit.size() << " collected) for " << params.rate_fit_points << "-point fitting." << std::endl;
                }
            }
        }
    } else {
        if (!params.noconsole) logger << "Could not calculate initial FFT peak parameters: " << initial_peak_params.message << std::endl;
    }
    
    // Apply drange/rrange to the refined values AFTER quadratic fitting from initial FFT,
    // as fitting might push the peak outside the original search box if the true peak is near the edge.
    // This ensures the starting point for iteration respects the user-defined search boundaries.
    if (initial_peak_params.success) { // Only clip if initial peak finding and fitting were attempted
        if (params.delay_search_range_specified) {
            if (refined_delay_from_first_fit_samples < params.delay_search_min) {
                if (!params.noconsole) logger << "  Refined Initial Delay from fit " << refined_delay_from_first_fit_samples << " samp (post-fit) clipped to min search range: " << params.delay_search_min << " samp." << std::endl;
                refined_delay_from_first_fit_samples = static_cast<float>(params.delay_search_min);
            } else if (refined_delay_from_first_fit_samples > params.delay_search_max) {
                if (!params.noconsole) logger << "  Refined Initial Delay from fit " << refined_delay_from_first_fit_samples << " samp (post-fit) clipped to max search range: " << params.delay_search_max << " samp." << std::endl;
                refined_delay_from_first_fit_samples = static_cast<float>(params.delay_search_max);
            }
        }
        if (params.rate_search_range_specified) {
            if (refined_rate_from_first_fit_hz < params.rate_search_min) {
                if (!params.noconsole) logger << "  Refined Initial Rate from fit " << refined_rate_from_first_fit_hz << " Hz (post-fit) clipped to min search range: " << params.rate_search_min << " Hz." << std::endl;
                refined_rate_from_first_fit_hz = static_cast<float>(params.rate_search_min);
            } else if (refined_rate_from_first_fit_hz > params.rate_search_max) {
                if (!params.noconsole) logger << "  Refined Initial Rate from fit " << refined_rate_from_first_fit_hz << " Hz (post-fit) clipped to max search range: " << params.rate_search_max << " Hz." << std::endl;
                refined_rate_from_first_fit_hz = static_cast<float>(params.rate_search_max);
            }
        }
        if (!params.noconsole && (params.delay_search_range_specified || params.rate_search_range_specified)) {
             logger << "After post-fit clipping by drange/rrange, initial values for iteration are:" << std::endl;
             logger << "  Rate for Iteration Start: " << refined_rate_from_first_fit_hz << " Hz, Delay for Iteration Start: " << refined_delay_from_first_fit_samples << " samples" << std::endl;
        }
    }
        

        // === Step 4: Iterative Calibration Loop ===
        float accumulated_rate_correction_hz = refined_rate_from_first_fit_hz;
        float accumulated_delay_correction_samples = refined_delay_from_first_fit_samples;
        

        if (params.initial_rate_specified) {
            accumulated_rate_correction_hz = params.initial_rate_hz;
            if (!params.noconsole) logger << "Overriding initial rate for iteration with command-line specified rate: " << params.initial_rate_hz << " Hz." << std::endl;
        }
        if (params.initial_delay_specified) {
            accumulated_delay_correction_samples = params.initial_delay_samples;
            if (!params.noconsole) logger << "Overriding initial delay for iteration with command-line specified delay: " << params.initial_delay_samples << " samples." << std::endl;
        }

        if (initial_peak_params.success && params.iterations > 0 && !params.frequency_mode && !params.quick_mode) {
            if (!params.noconsole) {
                logger << "\n--- Starting Iterative Calibration (" << params.iterations << " iterations) ---" << std::endl;
            }

            for (int iter = 0; iter < params.iterations; ++iter) {
                if (!params.noconsole) {
                    logger << "\n-- Iteration " << iter + 1 << "/" << params.iterations << " --" << std::endl;
                    logger << "Applying Total Rate: " << accumulated_rate_correction_hz 
                           << " Hz, Total Delay: " << accumulated_delay_correction_samples 
                           << " samples for this iteration." << std::endl;
                }
                
                std::vector<std::vector<float>> iter_fft_shifted_amplitude = generate_shifted_amplitude_plane(
                    current_all_spectrums,
                    accumulated_rate_correction_hz,
                    accumulated_delay_correction_samples,
                    fft_in, fft_intermediate, fft_out,
                    plan_vertical, plan_horizontal,
                    N_rows_original, N_cols_original,
                    N_rows_padded, N_cols_padded_fft,
                    rfi_index_ranges,
                    loaded_data, first_effective_integration_length
                );
            if (params.enable_text_log_output && !params.output_dir_final.empty()) {
                    std::string iter_fft_heatmap_filename = (fs::path(params.output_dir_final) / (fs::path(params.input_filename).stem().string() + "_loop" + std::to_string(loop_idx+1) + "_iter" + std::to_string(iter+1) + "_fft_heatmap.csv")).string();
                    write_2d_float_vector_to_csv(iter_fft_shifted_amplitude, iter_fft_heatmap_filename, logger, params.noconsole);
                }

                FftPeakParameters iter_peak_params;
                if (!iter_fft_shifted_amplitude.empty()) {
                    iter_peak_params = calculate_fft_peak_parameters(
                        iter_fft_shifted_amplitude, // Shifted amplitude for peak finding
                        fft_out,                    // Unshifted complex data from generate_shifted_a
                        first_effective_integration_length,
                        N_rows_padded,
                        N_cols_padded_fft,
                        params.delay_search_min - initial_peak_params.physical_delay_samples,
                        params.delay_search_max - initial_peak_params.physical_delay_samples,
                        true,
                        params.rate_search_min - initial_peak_params.physical_rate_hz,
                        params.rate_search_max - initial_peak_params.physical_rate_hz,
                        true
                    );
                } else {
                    iter_peak_params.success = false;
                    iter_peak_params.message = "Iterative FFT shifted amplitude data is empty.";
                }

                if (iter_peak_params.success) {
                    if (!params.noconsole) {
                        logger << "  Iter " << iter + 1 << " FFT - Max Amp (%): " <<  100* iter_peak_params.max_amplitude
                               << ", SNR: " << iter_peak_params.snr
                               << ", Residual Rate (Hz): " << iter_peak_params.physical_rate_hz
                               << ", Residual Delay (samp): " << iter_peak_params.physical_delay_samples << std::endl;
                    }

                    float residual_delay_from_fit_samples = iter_peak_params.physical_delay_samples;
                    float residual_rate_from_fit_hz = iter_peak_params.physical_rate_hz;

                    // --- Quadratic fit for Delay residual ---
                    std::vector<double> x_coords_iter_delay_fit;
                    std::vector<double> y_values_iter_delay_fit;
                    bool all_points_valid_iter_delay = true;
                    int iter_delay_half_steps = (params.delay_fit_points - 1) / 2;

                    if (!params.noconsole) logger << "    Fitting Delay residual (around peak i=" << iter_peak_params.max_i << ", j=" << iter_peak_params.max_j << "):" << std::endl;
                    for (int dj = -iter_delay_half_steps; dj <= iter_delay_half_steps; ++dj) {
                        int current_j = iter_peak_params.max_j + dj;
                        if (current_j >= 0 && current_j < N_cols_padded_fft && static_cast<size_t>(iter_peak_params.max_i) < iter_fft_shifted_amplitude.size() && static_cast<size_t>(current_j) < iter_fft_shifted_amplitude[static_cast<size_t>(iter_peak_params.max_i)].size()) {
                            // Consistent with updated sign convention in frinZparamcal.cpp
                            float physical_delay = (static_cast<float>(N_cols_padded_fft) / 2.0f) - static_cast<float>(current_j);
                            float amplitude = iter_fft_shifted_amplitude[static_cast<size_t>(iter_peak_params.max_i)][static_cast<size_t>(current_j)];
                            if (!params.noconsole) logger << "      Delay[j=" << current_j << "] (PhysRes: " << physical_delay << " samp): Amp = " << amplitude << std::endl;
                            x_coords_iter_delay_fit.push_back(static_cast<double>(physical_delay));
                            y_values_iter_delay_fit.push_back(static_cast<double>(amplitude));
                        } else { all_points_valid_iter_delay = false; if (!params.noconsole) logger << "      Delay[j=" << current_j << "]: Out of bounds" << std::endl;}
                    }
                    if (all_points_valid_iter_delay && x_coords_iter_delay_fit.size() == static_cast<size_t>(params.delay_fit_points)) {
                        QuadraticFitResult iter_delay_fit_result = fit_quadratic_least_squares(x_coords_iter_delay_fit, y_values_iter_delay_fit);
                        if (iter_delay_fit_result.success) {
                            residual_delay_from_fit_samples = static_cast<float>(iter_delay_fit_result.peak_x);
                            if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Delay): Refined Residual Delay (samp) = " << residual_delay_from_fit_samples << std::endl;
                        } else { if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Delay) failed: " << iter_delay_fit_result.message << ". Using non-fitted residual." << std::endl; }
                    } else { if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Delay): Not enough valid points for " << params.delay_fit_points << "-point fit. Using non-fitted residual." << std::endl; }
                    
                    accumulated_delay_correction_samples += residual_delay_from_fit_samples;


                    // Clip accumulated delay to user-defined search range
                    if (params.delay_search_range_specified) {
                        if (accumulated_delay_correction_samples < params.delay_search_min) {
                            if (!params.noconsole) {
                                logger << "  Accumulated Delay " << accumulated_delay_correction_samples 
                                       << " samp clipped to min search range: " << params.delay_search_min << " samp." << std::endl;
                            }
                            accumulated_delay_correction_samples = static_cast<float>(params.delay_search_min);
                        } else if (accumulated_delay_correction_samples > params.delay_search_max) {
                            if (!params.noconsole) {
                                logger << "  Accumulated Delay " << accumulated_delay_correction_samples 
                                       << " samp clipped to max search range: " << params.delay_search_max << " samp." << std::endl;
                            }
                            accumulated_delay_correction_samples = static_cast<float>(params.delay_search_max);
                        }
                    }
                    

                    // --- Quadratic fit for Rate residual ---
                    std::vector<double> x_coords_iter_rate_fit;
                    std::vector<double> y_values_iter_rate_fit;
                    bool all_points_valid_iter_rate = true;
                    int iter_rate_half_steps = (params.rate_fit_points - 1) / 2;
                    if (!params.noconsole) logger << "    Fitting Rate residual (around peak i=" << iter_peak_params.max_i << ", j=" << iter_peak_params.max_j << "):" << std::endl;
                    for (int di = -iter_rate_half_steps; di <= iter_rate_half_steps; ++di) {
                        int current_i = iter_peak_params.max_i + di;
                        if (current_i >= 0 && current_i < N_rows_padded && static_cast<size_t>(current_i) < iter_fft_shifted_amplitude.size() && static_cast<size_t>(iter_peak_params.max_j) < iter_fft_shifted_amplitude[static_cast<size_t>(current_i)].size()) {
                            float physical_rate = (static_cast<float>(current_i) - static_cast<float>(N_rows_padded) / 2.0f) * ((std::abs(first_effective_integration_length) > 1e-9f && N_rows_padded > 0) ? (2.0f * (1.0f / (2.0f * first_effective_integration_length)) / static_cast<float>(N_rows_padded)) : 0.0f);
                            float amplitude = iter_fft_shifted_amplitude[static_cast<size_t>(current_i)][static_cast<size_t>(iter_peak_params.max_j)];
                            if (!params.noconsole) logger << "      Rate[i=" << current_i << "] (PhysRes: " << physical_rate << " Hz): Amp = " << amplitude << std::endl;
                            x_coords_iter_rate_fit.push_back(static_cast<double>(physical_rate));
                            y_values_iter_rate_fit.push_back(static_cast<double>(amplitude));
                        } else { all_points_valid_iter_rate = false; if (!params.noconsole) logger << "      Rate[i=" << current_i << "]: Out of bounds" << std::endl; }
                    }
                    if (all_points_valid_iter_rate && x_coords_iter_rate_fit.size() == static_cast<size_t>(params.rate_fit_points)) {
                        // --- Rateフィッティングのためのスケーリング (イテレーション中) ---
                        std::vector<double> scaled_x_coords_iter_rate_fit = x_coords_iter_rate_fit;
                        const double iter_rate_scale_factor = 100.0; // 初期フィッティングと同じファクターを使用
                        for (double& x_val : scaled_x_coords_iter_rate_fit) {
                            x_val *= iter_rate_scale_factor;
                        }
                        
                        QuadraticFitResult iter_rate_fit_result = fit_quadratic_least_squares(scaled_x_coords_iter_rate_fit, y_values_iter_rate_fit);

                        if (iter_rate_fit_result.success) {
                            // peak_x と係数を元のスケールに戻す
                            iter_rate_fit_result.peak_x /= iter_rate_scale_factor;
                            iter_rate_fit_result.a *= (iter_rate_scale_factor * iter_rate_scale_factor);
                            iter_rate_fit_result.b *= iter_rate_scale_factor;
                            residual_rate_from_fit_hz = static_cast<float>(iter_rate_fit_result.peak_x);
                            if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Rate): Refined Residual Rate (Hz) = " << residual_rate_from_fit_hz 
                                                          << " (Coeffs a="<< iter_rate_fit_result.a << ", b=" << iter_rate_fit_result.b << ", c=" << iter_rate_fit_result.c << ")" << std::endl;
                        } else { if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Rate) failed: " << iter_rate_fit_result.message << ". Using non-fitted residual." << std::endl; }
                    } else { if (!params.noconsole) logger << "      Quadratic Fit (Iter " << iter + 1 << " Rate): Not enough valid points for " << params.rate_fit_points << "-point fit. Using non-fitted residual." << std::endl; }

                    accumulated_rate_correction_hz += residual_rate_from_fit_hz;

                    // Clip accumulated rate to user-defined search range
                    if (params.rate_search_range_specified) {
                        if (accumulated_rate_correction_hz < params.rate_search_min) {
                            if (!params.noconsole) {
                                logger << "  Accumulated Rate " << accumulated_rate_correction_hz 
                                       << " Hz clipped to min search range: " << params.rate_search_min << " Hz." << std::endl;
                            }
                            accumulated_rate_correction_hz = static_cast<float>(params.rate_search_min);
                        } else if (accumulated_rate_correction_hz > params.rate_search_max) {
                            if (!params.noconsole) {
                                logger << "  Accumulated Rate " << accumulated_rate_correction_hz 
                                       << " Hz clipped to max search range: " << params.rate_search_max << " Hz." << std::endl;
                            }
                            accumulated_rate_correction_hz = static_cast<float>(params.rate_search_max);
                        }
                    }

                    if (!params.noconsole && params.rate_search_range_specified) { // DEBUG log after rate clipping in iter
                        logger << "DEBUG Iter Post-Clip Rate: accumulated_rate_correction_hz=" << accumulated_rate_correction_hz << std::endl;
                    }

                    if (!params.noconsole) {
                        logger << "  Iter " << iter + 1 << " - Updated Total Applied Rate: " << accumulated_rate_correction_hz
                               << " Hz, Updated Total Applied Delay: " << accumulated_delay_correction_samples << " samples." << std::endl;
                    }

                } else {
                    if (!params.noconsole) logger << "  Iter " << iter + 1 << " - Could not calculate peak parameters: " << iter_peak_params.message << ". Stopping iteration." << std::endl;
                    break; // Stop iterating if peak finding fails
                }
            } // End of iteration loop

            if (!params.noconsole) {
                 logger << "\n--- Iterative Calibration Finished ---" << std::endl;
                 logger << "Final Accumulated Rate Correction: " << accumulated_rate_correction_hz << " Hz" << std::endl;
                 logger << "Final Accumulated Delay Correction: " << accumulated_delay_correction_samples << " samples" << std::endl;
            }
        }
        // End of Iterative Calibration
        fftwf_destroy_plan(plan_vertical);
        fftwf_destroy_plan(plan_horizontal);
        fftwf_free(fft_in);
        fftwf_free(fft_intermediate);
        fftwf_free(fft_out);
    
        } 
    
    
    //fftwf_cleanup_threads();
    

    return 0;
}
