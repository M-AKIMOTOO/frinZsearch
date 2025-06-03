#include "frinZargs.hpp"
#include <iostream> // For std::cerr
#include <filesystem> // For directory creation

void setup_cli_options(CLI::App& app, ProgramOptions& params) {
    app.add_option("--input", params.input_filename, "Specify the input binary file")
        ->required()
        ->check(CLI::ExistingFile);
    app.add_option("--length", params.specified_length_sec, "Specify integration time for FFT segment (seconds, default: all data in segment)");
    app.add_option("--skip", params.skip_seconds, "Specify time to skip from the beginning of the data (seconds, default: 0.0)")
        ->check(CLI::NonNegativeNumber);
    app.add_option("--loop", params.loop_count, "Loop FFT process for segments (default: 1)")
        ->check(CLI::PositiveNumber);
    app.add_option("--iter", params.iterations, "Number of phase calibration iterations (default: 5)") // Option name "iter"
        ->check(CLI::NonNegativeNumber);
    app.add_flag("--output", params.enable_text_log_output, "Enable output of console messages and data files to a 'frinZ/frinZsearch' subdirectory relative to the input file. Log filename: <input_basename>_frinZsearch_result.txt");
    app.add_flag("--frequency", params.frequency_mode, "Perform only the first FFT, no calibration/iteration");
    app.add_option_function<std::vector<std::string>>("--rfi", [&](const std::vector<std::string> &ranges) {
        for (const auto& range_str : ranges) {
            size_t comma_pos = range_str.find(',');
            if (comma_pos != std::string::npos) {
                try {
                    double start_mhz = std::stod(range_str.substr(0, comma_pos));
                    double end_mhz = std::stod(range_str.substr(comma_pos + 1));
                    if (start_mhz < 0 || end_mhz < 0 || start_mhz > end_mhz) {
                        throw CLI::ValidationError("--rfi", "Invalid RFI range: " + range_str + ". Frequencies must be non-negative and start <= end.");
                    }
                    params.rfi_ranges_mhz.emplace_back(start_mhz, end_mhz);
                } catch (const std::exception& e) {
                    throw CLI::ValidationError("--rfi", "Invalid RFI range format: " + range_str + ". Expected format: start_mhz,end_mhz. Error: " + e.what());
                }
            } else {
                throw CLI::ValidationError("--rfi", "Invalid RFI range format: " + range_str + ". Expected comma-separated pair: start_mhz,end_mhz");
            }
        }
    }, "Specify RFI frequency ranges to cut (MHz). E.g., --rfi 10,50 100,120")
        ->type_name("RFI_RANGE (MHz)");
    app.add_option("--delay", params.initial_delay_samples, "Initial delay correction in samples. Sets --iter to 0 if --iter not specified (default: 0.0)"); // Option name "delay"
    app.add_option("--rate", params.initial_rate_hz, "Initial rate correction in Hz. Sets --iter to 0 if --iter not specified (default: 0.0)");    // Option name "rate"
    
    // For --drange and --rrange, CLI11 will store results in these local vectors.
    // We'll retrieve them in finalize_options using app.get_option("...").
    // These dummy vectors are just to satisfy add_option's signature here.
    static std::vector<double> drange_vals_dummy; // static to have a stable address if CLI::App stores a pointer internally, though not ideal.
    static std::vector<double> rrange_vals_dummy; // A better way would be callbacks or direct binding if ProgramOptions members were public and simple.

    app.add_option("--drange", drange_vals_dummy, "Limit delay search range for calibration [min_samp max_samp]") // Option name "drange"
        ->expected(2);
    app.add_option("--rrange", rrange_vals_dummy, "Limit rate search range for calibration [min_Hz max_Hz]") // Option name "rrange"
        ->expected(2);
    app.add_option("--dstep", params.delay_fit_points, "Number of points for delay quadratic fit (e.g., 3, 5, 7). Default: 3")
        ->check([](const std::string& str) {
            int val = std::stoi(str);
            if (val < 3 || val % 2 == 0) return std::string("Value must be an odd integer >= 3.");
            return std::string();
        });
    app.add_option("--rstep", params.rate_fit_points, "Number of points for rate quadratic fit (e.g., 3, 5, 7). Default: 3")
        ->check([](const std::string& str) {
            int val = std::stoi(str);
            if (val < 3 || val % 2 == 0) return std::string("Value must be an odd integer >= 3.");
            return std::string();
        });
    app.add_flag("--header", params.output_header_info, "Output header information from the input file");
    app.add_flag("--quick", params.quick_mode, "Skip phase calibration, output initial FFT results only");
    app.add_flag("--noconsole", params.noconsole, "Suppress console output (still writes to files if --output is set)");
    app.set_version_flag("--version", "frinZsearch v0.1 in 2025/05/31 (Made by M.AKIMOTO with Gemini)\nfrinZsearch v1.0 in 2025/06/01\nfrinZsearch v2.0 in 2025/05/03 ");
}

void finalize_options(ProgramOptions& params, CLI::App& app) {
    // Retrieve option counts or values using app.get_option()
    if (app.get_option("--iter")->count() > 0) {
        params.iter_explicitly_set = true;
    }
    if (app.get_option("--delay")->count() > 0) {
        params.initial_delay_specified = true;
    }
    if (app.get_option("--rate")->count() > 0) {
        params.initial_rate_specified = true;
    }

    CLI::Option* drange_opt = app.get_option("--drange");
    if (drange_opt->count() > 0) {
        const std::vector<double>& drange_vals = drange_opt->as<std::vector<double>>();
        // Validation for size 2 is implicitly handled by ->expected(2) during parsing.
        // If parsing succeeded, drange_vals will have 2 elements.
        params.delay_search_min = drange_vals[0];
        params.delay_search_max = drange_vals[1];
        if (params.delay_search_min > params.delay_search_max) {
             throw CLI::ValidationError("--drange", "min_samp (" + std::to_string(params.delay_search_min) + ") > max_samp (" + std::to_string(params.delay_search_max) + ") for --drange.");
        }
        params.delay_search_range_specified = true;
    }

    CLI::Option* rrange_opt = app.get_option("--rrange");
    if (rrange_opt->count() > 0) {
        const std::vector<double>& rrange_vals = rrange_opt->as<std::vector<double>>();
        params.rate_search_min = rrange_vals[0];
        params.rate_search_max = rrange_vals[1];
         if (params.rate_search_min > params.rate_search_max) {
            throw CLI::ValidationError("--rrange", "min_Hz (" + std::to_string(params.rate_search_min) + ") > max_Hz (" + std::to_string(params.rate_search_max) + ") for --rrange.");
        }
        params.rate_search_range_specified = true;
    }

    // Adjust iterations based on other flags.
    // If --delay or --rate is specified, iterations are NOT automatically set to 0.
    // They will use the default (e.g., 5) unless --iter is also specified.
    // if (params.initial_delay_specified || params.initial_rate_specified) {
    //     if (!params.iter_explicitly_set) {
    //         // params.iterations = 0; // Keep default iterations unless explicitly set by --iter
    //     }
    // }
    if (params.frequency_mode || params.quick_mode) {
        params.iterations = 0;
    }

    // Setup output directory path
    if (params.enable_text_log_output) {
        if (params.input_filename.empty()) { 
            // This should have been caught by --input->required() during app.parse()
            std::cerr << "Critical Error: --output enabled, but input filename is empty after parsing. This indicates an issue." << std::endl;
            params.enable_text_log_output = false; // Disable to prevent further errors
        } else {
            namespace fs = std::filesystem;
            fs::path input_file_path_fs(params.input_filename);
            fs::path absolute_input_path = fs::absolute(input_file_path_fs);
            fs::path base_output_dir = absolute_input_path.parent_path();
            
            params.output_dir_final = (base_output_dir / "frinZ" / "frinZsearch").string();
            try {
                fs::create_directories(params.output_dir_final);
            } catch (const fs::filesystem_error& e) {
                std::cerr << "Error creating output directory '" << params.output_dir_final << "': " << e.what() << std::endl;
                params.enable_text_log_output = false; // Disable if dir creation fails
                params.output_dir_final.clear();
            }
        }
    }
}