//
// Made by M.AKIMOTO with Gemini 
// 2025/06/01
// g++ -std=c++17 -Wall -O2  -c frinZsearch.cpp -o frinZsearch.o
// g++ frinZsearch.o -o frinZsearch  -lfftw3f -lfftw3f_threads -lm -lpthread -lstdc++fs
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
#include <thread>     // For std::thread::hardware_concurrency
#include "frinZargs.hpp" // Include the new header for ProgramOptions and arg functions
#include "frinZread.hpp" // Include the new header for file reading functions and HeaderRegion
#include "frinZlogger.hpp"    // Include the new Logger header
#include "frinZparamcal.hpp" // For parameter calculation
#include "frinZfftshift.hpp" // For fftshift and amplitude calculation


int main(int argc, char* argv[]) {

    ProgramOptions params;
    Logger logger;
    

    CLI::App app{"frinZsearch: Fringe Z Search Program"};
    // Set up command line options using the function from frinZargs.cpp
    setup_cli_options(app, params);

    // If no arguments (other than program name) are given, print help and exit.
    if (argc == 1) {
        std::cout << app.help() << std::endl;
        return 0;
    }


    // Parse the arguments. CLI11 will handle --help, --version, and errors.
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e); // Prints error message and exits
    }

    // Finalize options based on parsed values (e.g., dependencies between options)
    try {
        finalize_options(params, app);
    } catch (const CLI::ValidationError &e) { // Catch validation errors from finalize_options
        std::cerr << "Validation Error: " << e.what() << std::endl;
        return 1;
    }

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

        int N_rows_original = static_cast<int>(current_all_spectrums.size());
        int N_cols_original = static_cast<int>(current_all_spectrums[0].size()); // Assuming consistent inner size

        // === ゼロパディング計算 (縦軸) === fft_point を使用
        int N_rows_padded = static_cast<int>(fft_point); // ヘッダーから読み込んだ fft_point を縦軸のFFT点数とする


        // === ゼロパディング計算 (横軸) ===
        int N_cols_padded_fft = N_cols_original*2; 
        if (N_cols_original > 0) {
            double temp_cols = static_cast<double>(N_cols_original);
            int n_power_cols = 0;
            while (std::pow(2, n_power_cols) < temp_cols) { // Find smallest power of 2 >= temp_cols
                n_power_cols++;
            }
            N_cols_padded_fft = static_cast<int>(std::pow(2, n_power_cols)); // Round up to power of 2
            N_cols_padded_fft = static_cast<int>(std::round(N_cols_padded_fft * pad_factor_cols_val));
        }
        int N_cols_fft = N_cols_padded_fft;

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
    for (int i = 0; i < rows_to_copy_to_fft_in; ++i) { // コピーする行数を N_rows_padded 以下に制限
        for (int j = 0; j < N_cols_original; ++j) { // N_cols_original は fft_half と同じはず
            bool is_rfi_channel = false;
            for (const auto& rfi_idx_range : rfi_index_ranges) {
                if (static_cast<int>(j) >= rfi_idx_range.first && static_cast<int>(j) <= rfi_idx_range.second) {
                    is_rfi_channel = true;
                    break;
                }
            }

            if (is_rfi_channel) {
                fft_in[i * N_cols_padded_fft + j][0] = 0.0f;
                fft_in[i * N_cols_padded_fft + j][1] = 0.0f;
            } else {
                fft_in[i * N_cols_padded_fft + j][0] = current_all_spectrums[static_cast<size_t>(i)][static_cast<size_t>(j)].real();
                fft_in[i * N_cols_padded_fft + j][1] = current_all_spectrums[static_cast<size_t>(i)][static_cast<size_t>(j)].imag();
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
    fftwf_destroy_plan(plan_vertical);

    // === Normalization after vertical FFT ===
    const float norm_factor = static_cast<float>(loaded_data.header.fft_point) / static_cast<float>(N_rows_padded * N_cols_fft);
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
                                          FFTW_FORWARD, fftw_plan_flag_val);

    fftwf_execute(plan_horizontal);
    fftwf_destroy_plan(plan_horizontal);

    // === Step 3: fftshiftと振幅計算 ===
    std::vector<std::vector<float>> fft_shifted_amplitude = perform_fftshift_and_calc_amplitude(fft_out, N_rows_padded, N_cols_padded_fft);
    

    // === Step 3.5: FFT後の振幅の最大値とその座標を検索 ===
    FftPeakParameters initial_peak_params;
    if (!fft_shifted_amplitude.empty()) {
        initial_peak_params = calculate_fft_peak_parameters(
            fft_shifted_amplitude,
            first_effective_integration_length,
            N_rows_padded,
            N_cols_padded_fft
        );
    } else {
        initial_peak_params.success = false;
        initial_peak_params.message = "Initial FFT shifted amplitude data is empty.";
    }

    if (initial_peak_params.success) {
        if (!params.noconsole) {
            logger << "Initial FFT - StdDev: " << initial_peak_params.stddev_amplitude 
                   << ", SNR (Max/StdDev): " << initial_peak_params.snr << std::endl;

            if (!params.noconsole) {
                logger << std::fixed;
                logger.setprecision(8) << "Maximum FFT amplitude: " << initial_peak_params.max_amplitude
                                       << " found at Fringe Frequency (Hz): " << initial_peak_params.physical_rate_hz
                                       << ", Delay (Sample): " << initial_peak_params.physical_delay_samples << std::endl;
                logger << std::dec; // Reset to decimal mode separately
            }
            
        }
    } else {
        if (!params.noconsole) {
            logger << "Could not calculate initial FFT peak parameters: " << initial_peak_params.message << std::endl;
        }
    }

        // === Iterative Search for Optimal Calibration Parameters ===
        float best_rate_for_calibration = params.initial_rate_specified ? params.initial_rate_hz : initial_peak_params.physical_rate_hz;
        float best_delay_for_calibration = params.initial_delay_specified ? params.initial_delay_samples : initial_peak_params.physical_delay_samples;
        float highest_peak_amplitude_after_calibration = -1.0f; // Initialize with a value lower than any possible amplitude

         if (params.quick_mode || params.frequency_mode) {
            if (!params.noconsole) logger << "Skipping iterative search due to --quick or --frequency flag." << std::endl;
        } else if (params.iterations == 0) {
            if (!params.noconsole) logger << "Skipping iterative search as iterations set to 0." << std::endl;
        } else if (initial_peak_params.success && initial_peak_params.max_i != -1 && !current_all_spectrums.empty()) { // Only proceed if an initial peak was found and data exists
            if (!params.noconsole) {
                logger << "Starting iterative search for optimal calibration parameters..." << std::endl;
                logger << "Initial parameters for search: Rate=" << best_rate_for_calibration << " Hz, Delay=" << best_delay_for_calibration << " Samples" << std::endl;
            }
            std::vector<std::vector<std::complex<float>>> all_spectrums_original_copy = current_all_spectrums; // Keep original data for iteration

            float rate_resolution_hz = 0.0f;
            if (std::abs(first_effective_integration_length) > 1e-9f && N_rows_padded > 0) {
                float max_fringe_freq_val_iter = 1.0f / (2.0f * first_effective_integration_length);
                rate_resolution_hz = (2.0f * max_fringe_freq_val_iter / static_cast<float>(N_rows_padded));
            }
            if (!params.noconsole) logger << "Base Rate resolution for search: " << rate_resolution_hz << " Hz" << std::endl;

            float current_delay_search_half_range = 0.5f; // Initial delay search half-width
            float current_delay_search_step = 0.1f;       // Initial delay search step
            float current_rate_search_range_abs = (rate_resolution_hz > 1e-9f) ? rate_resolution_hz : 0.0f; // Initial rate search width
            float current_rate_search_step = (rate_resolution_hz > 1e-9f) ? rate_resolution_hz / 10.0f : 0.0f; // Initial rate search step

            fftwf_complex *fft_in_iter, *fft_intermediate_iter, *fft_out_iter;
            fft_in_iter = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
            fft_intermediate_iter = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
            fft_out_iter = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);

            if (!fft_in_iter || !fft_intermediate_iter || !fft_out_iter) {
                std::cerr << "Error: Failed to allocate memory for iterative FFT." << std::endl;
            } else {
                // Create plans for iterative FFT (similar to calibrated FFT plans)
                fftwf_plan plan_vertical_iter = fftwf_plan_many_dft(rank_vert, n_vert, howmany_vert, fft_in_iter, nullptr, istride_vert, idist_vert, fft_intermediate_iter, nullptr, ostride_vert, odist_vert, FFTW_FORWARD, fftw_plan_flag_val);
                fftwf_plan plan_horizontal_iter = fftwf_plan_many_dft(rank_horiz, n_horiz, howmany_horiz, fft_intermediate_iter, nullptr, istride_horiz, idist_horiz, fft_out_iter, nullptr, ostride_horiz, odist_horiz, FFTW_FORWARD, fftw_plan_flag_val);

                for (int iter_count = 0; iter_count < params.iterations; ++iter_count) {
                    if (!params.noconsole) {
                        logger << "--- Iteration " << iter_count + 1 << "/" << params.iterations << " ---" << std::endl;
                        logger << "Searching around Rate: " << best_rate_for_calibration << " Hz, Delay: " << best_delay_for_calibration << " Samples" << std::endl;
                        logger << "Rate search range: +/-" << current_rate_search_range_abs << " Hz (step: " << current_rate_search_step << " Hz)" << std::endl;
                        logger << "Delay search range: +/-" << current_delay_search_half_range << " Samples (step: " << current_delay_search_step << " Samples)" << std::endl;
                    }
                    // The center for this iteration's search is the best found so far
                    float center_delay_for_iter = best_delay_for_calibration;
                    float center_rate_for_iter = best_rate_for_calibration;

                    for (float delay_offset = -current_delay_search_half_range; delay_offset <= current_delay_search_half_range + current_delay_search_step * 0.5f; delay_offset += current_delay_search_step) {
                        float current_delay_candidate = center_delay_for_iter + delay_offset;
                        if (params.delay_search_range_specified && (current_delay_candidate < params.delay_search_min || current_delay_candidate > params.delay_search_max)) {
                            continue; // Skip if outside user-defined delay range
                        }
                        float rate_loop_start = (current_rate_search_step > 1e-9f) ? -current_rate_search_range_abs : 0.0f;
                        float rate_loop_end = (current_rate_search_step > 1e-9f) ? current_rate_search_range_abs + current_rate_search_step * 0.5f : 0.0f;
                        float rate_loop_step = (current_rate_search_step > 1e-9f) ? current_rate_search_step : 1.0f;
                        if (rate_loop_step < 1e-9f && rate_loop_start == 0.0f && rate_loop_end == 0.0f) rate_loop_step = 1.0f; // Ensure loop runs once if range is zero

                        for (float rate_offset = rate_loop_start; rate_offset <= rate_loop_end; rate_offset += rate_loop_step) {
                            float current_rate_candidate = center_rate_for_iter + rate_offset;

                            if (params.rate_search_range_specified && (current_rate_candidate < params.rate_search_min || current_rate_candidate > params.rate_search_max)) {
                                if (rate_loop_step < 1e-9f) break; // if step is effectively zero and we are out of range, break inner loop
                                continue; // Skip if outside user-defined rate range
                            }
                            if (rate_loop_step < 1e-9f && (rate_offset > rate_loop_end + 1e-9f || rate_offset < rate_loop_start - 1e-9f)) { // Safety for zero step
                                 break; // Avoid infinite loop if step is zero and we are somehow outside the single point check
                            }

                            std::vector<std::vector<std::complex<float>>> temp_spectrums = all_spectrums_original_copy;

                            const float eff_integ_length_sec_iter = first_effective_integration_length;
                            const float sampling_freq_hz_iter = static_cast<float>(loaded_data.header.sampling_speed);
                            const uint32_t total_fft_points_iter = loaded_data.header.fft_point;
                            const uint32_t num_freq_channels_half_iter = total_fft_points_iter / 2;
                            float py_equiv_sampling_speed_mhz_iter = sampling_freq_hz_iter / 1.0e6f;
                            float stop_val_for_linspace_mhz_iter = static_cast<float>(static_cast<int>(py_equiv_sampling_speed_mhz_iter / 2.0f)) - 1.0f;

                            for (size_t s = 0; s < temp_spectrums.size(); ++s) {
                                float time_for_rate_corr_sec = static_cast<float>(s) * eff_integ_length_sec_iter;
                                std::complex<float> rate_corr_factor = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * current_rate_candidate * time_for_rate_corr_sec));
                                for (size_t k = 0; k < num_freq_channels_half_iter; ++k) {
                                    float freq_k_hz_for_delay_corr;
                                    if (num_freq_channels_half_iter > 1) {
                                        freq_k_hz_for_delay_corr = (static_cast<float>(k) * stop_val_for_linspace_mhz_iter / (static_cast<float>(num_freq_channels_half_iter) - 1.0f)) * 1.0e6f;
                                    } else { freq_k_hz_for_delay_corr = 0.0f; }
                                    float delay_seconds = current_delay_candidate / sampling_freq_hz_iter;
                                    std::complex<float> delay_corr_factor = std::exp(std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * delay_seconds * freq_k_hz_for_delay_corr));
                                    if (s < temp_spectrums.size() && k < temp_spectrums[s].size()) {
                                       temp_spectrums[s][k] *= rate_corr_factor * delay_corr_factor;
                                    }
                            }
                        }

                        // Perform 2D FFT on calibrated temp_spectrums
                        for (int r_idx = 0; r_idx < N_rows_padded * N_cols_fft; ++r_idx) { fft_in_iter[r_idx][0] = 0.0f; fft_in_iter[r_idx][1] = 0.0f; }
                        for (size_t r_orig = 0; r_orig < static_cast<size_t>(N_rows_original); ++r_orig) {
                            for (size_t c_orig = 0; c_orig < static_cast<size_t>(N_cols_original); ++c_orig) {
                                if (r_orig < temp_spectrums.size() && c_orig < temp_spectrums[r_orig].size()) { // Comparisons are now between size_t
                                    fft_in_iter[r_orig * N_cols_padded_fft + c_orig][0] = temp_spectrums[r_orig][c_orig].real();
                                    fft_in_iter[r_orig * N_cols_padded_fft + c_orig][1] = temp_spectrums[r_orig][c_orig].imag();
                                }
                            }
                        }
                        fftwf_execute(plan_vertical_iter);
                        fftwf_execute(plan_horizontal_iter);


                        // FFTShift and find max amplitude
                        std::vector<std::vector<float>> trial_shifted_amplitude = perform_fftshift_and_calc_amplitude(fft_out_iter, N_rows_padded, N_cols_padded_fft);
                        
                        FftPeakParameters trial_peak_params;
                        if (!trial_shifted_amplitude.empty()) {
                             trial_peak_params = calculate_fft_peak_parameters(trial_shifted_amplitude, eff_integ_length_sec_iter, N_rows_padded, N_cols_padded_fft);
                        } else {
                            trial_peak_params.success = false;
                        }

                        if (trial_peak_params.success && trial_peak_params.max_amplitude > highest_peak_amplitude_after_calibration) {
                                highest_peak_amplitude_after_calibration = trial_peak_params.max_amplitude;
                                best_rate_for_calibration = current_rate_candidate;
                                best_delay_for_calibration = current_delay_candidate;
                            }
                            if (rate_loop_step < 1e-9f) break; 
                        } // end rate_offset loop
                    } // end delay_offset loop

                    // Reduce search range and step for the next iteration
                    current_delay_search_half_range /= 10.0f;
                    current_delay_search_step /= 10.0f;
                    current_rate_search_range_abs /= 10.0f;
                    current_rate_search_step /= 10.0f;

                    // Ensure steps do not become too small to be useful or cause infinite loops with float precision
                    if (current_delay_search_step < 1e-6f) current_delay_search_step = 1e-6f;
                    if (current_rate_search_step < 1e-9f && current_rate_search_range_abs > 1e-8f) current_rate_search_step = 1e-9f; // Avoid zero step if range is still non-zero
                    else if (current_rate_search_range_abs <= 1e-8f) current_rate_search_step = 0.0f; // Effectively stop rate search if range is tiny

                }
                fftwf_destroy_plan(plan_vertical_iter);
                fftwf_destroy_plan(plan_horizontal_iter);
                fftwf_free(fft_in_iter); fftwf_free(fft_intermediate_iter); fftwf_free(fft_out_iter);

                if (!params.noconsole) {
                    logger << "Iterative search complete." << std::endl;
                    logger << "Best parameters found: Rate=" << best_rate_for_calibration
                              << " Hz, Delay=" << best_delay_for_calibration
                              << " Samples, Max Calibrated Amp (from iter search): " << highest_peak_amplitude_after_calibration << std::endl;
                }
            }
        } else if (initial_peak_params.success && initial_peak_params.max_i != -1) { // Peak found, but all_spectrums might be empty (should not happen if max_i != -1)
            if (!params.noconsole) logger << "Skipping iterative search as data is empty, though an initial peak was found." << std::endl;
        } else { // No initial peak found
            if (!params.noconsole) logger << "Skipping iterative search as no initial peak was found." << std::endl;
        }

        // === Step 4 (Final): 位相較正の実行 (using best parameters) ===
        // This step is skipped if --quick or --frequency is set.
        // If iterations was 0 due to --delay/--rate, we still apply this one-time calibration with those values.
        if (params.quick_mode || params.frequency_mode) {
            if (!params.noconsole) logger << "Skipping final phase calibration due to --quick or --frequency flag." << std::endl;
        } else if (initial_peak_params.success && initial_peak_params.max_i != -1 && !current_all_spectrums.empty()) { // 最大値が見つかり、データが存在する場合
            if (!params.noconsole) logger << "Applying final phase calibration to all_spectrums data..." << std::endl;

            // current_all_spectrums is modified IN-PLACE here
            // Use the best parameters found from iteration (or initial if iteration was skipped/failed)
            const float rate_to_correct_hz = best_rate_for_calibration; // This will be initial if iter=0, or refined if iter > 0
            const float delay_samples_to_correct = best_delay_for_calibration; // Same as above


            const float eff_integ_length_sec = first_effective_integration_length;
            const float sampling_freq_hz = static_cast<float>(loaded_data.header.sampling_speed);
            const uint32_t total_fft_points = loaded_data.header.fft_point;
            const uint32_t num_freq_channels_half = total_fft_points / 2;

            if (sampling_freq_hz <= 0.0f) {
                std::cerr << "Warning: Sampling frequency is not positive. Cannot perform delay calibration." << std::endl;
            } else if (total_fft_points < 2) {
                std::cerr << "Warning: Not enough FFT points. Cannot perform delay calibration." << std::endl;
            } else {
                float py_equiv_sampling_speed_mhz = sampling_freq_hz / 1.0e6f;
                // Python's int() truncates towards zero. For positive numbers, same as static_cast<int>.
                float stop_val_for_linspace_mhz = static_cast<float>(static_cast<int>(py_equiv_sampling_speed_mhz / 2.0f)) - 1.0f;
                
                // Calibrate all_spectrums IN-PLACE for the final FFT output
                for (size_t s = 0; s < current_all_spectrums.size(); ++s) {
                    // Rate correction factor for this sector
                    float time_for_rate_corr_sec = static_cast<float>(s) * eff_integ_length_sec;
                    std::complex<float> rate_corr_factor = std::exp(
                        std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * rate_to_correct_hz * time_for_rate_corr_sec)
                    );

                    if (s >= current_all_spectrums.size() || current_all_spectrums[s].empty() || current_all_spectrums[s].size() != num_freq_channels_half) {
                        std::cerr << "Warning: Mismatch in expected number of frequency channels for sector " << s << std::endl;
                        continue;
                    }

                    for (size_t k = 0; k < num_freq_channels_half; ++k) {
                        // Delay correction factor for this frequency channel
                        float freq_k_hz_for_delay_corr;
                        if (num_freq_channels_half > 1) {
                            freq_k_hz_for_delay_corr = (static_cast<float>(k) * stop_val_for_linspace_mhz / (static_cast<float>(num_freq_channels_half) - 1.0f)) * 1.0e6f;
                        } else { // num_freq_channels_half == 1 (so k must be 0)
                            freq_k_hz_for_delay_corr = 0.0f; // Mimicking np.linspace(0, X, 1) -> [0.0]
                        }

                        float delay_seconds = delay_samples_to_correct / sampling_freq_hz;
                        std::complex<float> delay_corr_factor = std::exp(
                            std::complex<float>(0.0f, -2.0f * static_cast<float>(M_PI) * delay_seconds * freq_k_hz_for_delay_corr)
                        );

                        if (k < current_all_spectrums[s].size()){ // Bounds check
                            current_all_spectrums[s][k] *= rate_corr_factor * delay_corr_factor;
                        }
                    }
                }
                if (!params.noconsole) logger << "Phase calibration applied." << std::endl;
                // ここで較正された all_spectrums をファイルに書き出すなどの処理を追加できます。
                // 例:
                // std::cout << "First calibrated spectrum value: " << all_spectrums[0][0] << std::endl;
                // === Step 5: 較正後データの2D FFTと振幅計算 ===
                if (!params.noconsole) logger << "Performing 2D FFT on calibrated data..." << std::endl;

                fftwf_complex *fft_in_calibrated, *fft_intermediate_calibrated, *fft_out_calibrated;
                fft_in_calibrated = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
                fft_intermediate_calibrated = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);
                fft_out_calibrated = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N_rows_padded * N_cols_fft);

                if (!fft_in_calibrated || !fft_intermediate_calibrated || !fft_out_calibrated) {
                    std::cerr << "Error: Failed to allocate memory for calibrated FFT." << std::endl;
                    if (fft_in_calibrated) fftwf_free(fft_in_calibrated);
                    if (fft_intermediate_calibrated) fftwf_free(fft_intermediate_calibrated);
                    if (fft_out_calibrated) fftwf_free(fft_out_calibrated);
                } else {
                    // fft_in_calibrated をゼロで初期化
                    for (int r_idx = 0; r_idx < N_rows_padded * N_cols_fft; ++r_idx) {
                        fft_in_calibrated[r_idx][0] = 0.0f;
                        fft_in_calibrated[r_idx][1] = 0.0f;
                    }

                     // *最終的に較正された* all_spectrums データを fft_in_calibrated にコピー (ゼロパディング)
                    for (size_t r_orig = 0; r_orig < static_cast<size_t>(N_rows_original); ++r_orig) {
                        for (size_t c_orig = 0; c_orig < static_cast<size_t>(N_cols_original); ++c_orig) { // N_cols_original is from current_all_spectrums[0].size()
                            // Ensure r_orig and c_orig are within bounds of all_spectrums
                            if (r_orig < current_all_spectrums.size() && !current_all_spectrums[r_orig].empty() && c_orig < current_all_spectrums[r_orig].size()) { // Comparisons are now between size_t
                                fft_in_calibrated[r_orig * N_cols_padded_fft + c_orig][0] = current_all_spectrums[r_orig][c_orig].real() * norm_factor;
                                fft_in_calibrated[r_orig * N_cols_padded_fft + c_orig][1] = current_all_spectrums[r_orig][c_orig].imag() * norm_factor;
                            }
                        }
                    }

                    // 較正データに対する縦方向FFT (plan_vertical のパラメータを再利用)
                    fftwf_plan plan_vertical_calibrated = fftwf_plan_many_dft(
                        rank_vert, n_vert, howmany_vert,
                        fft_in_calibrated, nullptr, istride_vert, idist_vert,
                        fft_intermediate_calibrated, nullptr, ostride_vert, odist_vert,
                        FFTW_FORWARD, fftw_plan_flag_val);
                    fftwf_execute(plan_vertical_calibrated);
                    fftwf_destroy_plan(plan_vertical_calibrated);

                    // 較正データに対する横方向FFT (plan_horizontal のパラメータを再利用)
                    fftwf_plan plan_horizontal_calibrated = fftwf_plan_many_dft(
                        rank_horiz, n_horiz, howmany_horiz,
                        fft_intermediate_calibrated, nullptr, istride_horiz, idist_horiz,
                        fft_out_calibrated, nullptr, ostride_horiz, odist_horiz,
                        FFTW_FORWARD, fftw_plan_flag_val);
                    fftwf_execute(plan_horizontal_calibrated);
                    fftwf_destroy_plan(plan_horizontal_calibrated);

                    // 較正データのFFTShiftと振幅計算
                    std::vector<std::vector<float>> fft_shifted_amplitude_calibrated = perform_fftshift_and_calc_amplitude(fft_out_calibrated, N_rows_padded, N_cols_padded_fft);
                    
                    if (!params.noconsole) logger << "Calibrated data 2D FFT and amplitude calculation complete." << std::endl;

                    // Calculate and display SNR for the final calibrated data
                    FftPeakParameters final_calibrated_params;
                    if (!fft_shifted_amplitude_calibrated.empty()) {
                        final_calibrated_params = calculate_fft_peak_parameters(
                            fft_shifted_amplitude_calibrated,
                            first_effective_integration_length,
                            N_rows_padded,
                            N_cols_padded_fft
                        );
                    } else {
                        final_calibrated_params.success = false;
                        final_calibrated_params.message = "Final calibrated FFT shifted amplitude data is empty.";
                    }
                    if (final_calibrated_params.success && !params.noconsole) {
                        logger << "Final Calibrated FFT - Max Amp: " << final_calibrated_params.max_amplitude 
                               << ", StdDev: " << final_calibrated_params.stddev_amplitude 
                               << ", SNR: " << final_calibrated_params.snr << std::endl;
                    } else if (!params.noconsole) {
                        logger << "Could not calculate final calibrated FFT parameters: " << final_calibrated_params.message << std::endl;
                    }
                        
                    fftwf_free(fft_in_calibrated);
                    fftwf_free(fft_intermediate_calibrated);
                    fftwf_free(fft_out_calibrated);
                }
            }
        } else {
            if (!params.noconsole) logger << "Skipping final phase calibration and subsequent FFT as no initial peak was found, data is empty, or mode prevents it." << std::endl;
        }
    
        // Free FFTW arrays allocated in this loop iteration
        fftwf_free(fft_in);
        fftwf_free(fft_intermediate);
        fftwf_free(fft_out);
    
        } // End of main processing loop (for params.loop_count)
    
    
    // Cleanup FFTW threads support
    fftwf_cleanup_threads();


    return 0;
}
