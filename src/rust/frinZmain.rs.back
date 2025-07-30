use clap::{Parser, CommandFactory};
use std::path::PathBuf;
use std::f64::INFINITY;
use std::fs;


mod frinz_read;
mod frinz_fft;
mod frinz_fitting;
mod frinz_logger;
mod frinz_error;

/// Program to perform fringe search
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None, after_help = "(c) M.AKIMOTO with Gemini in 2025/07/17
This program is licensed under the MIT License
see https://opensource.org/license/mit")]
struct Args {
    /// Specify input binary file (required)
    #[arg(short, long, aliases = ["in", "inp", "i", "inpu"])]
    input: String,

    /// Specify integration time for FFT segment (seconds, default: all data)
    #[arg(long, alias = "len", default_value_t = -1.0)]
    length: f64,

    /// Specify time to skip from the beginning of the data (seconds, default: 0.0)
    #[arg(long, default_value_t = 0.0)]
    skip: f64,

    /// Number of times to loop FFT processing per segment (default: 1)
    #[arg(long, alias = "loop", default_value_t = 1)]
    loop_count: usize,

    /// Number of phase calibration iterations (default: 3)
    #[arg(long, default_value_t = 3)]
    iter: usize,

    /// Output console messages and data files to 'frinZ/frinZsearch' subdirectory relative to input file
    #[arg(long)]
    output: bool,

    /// Perform only the first FFT, no calibration/iteration
    #[arg(long)]
    frequency: bool,

    /// Specify RFI frequency range in MHz (e.g.: --rfi 10,50 --rfi 100,120)
    #[arg(long, value_parser = parse_rfi_range)]
    rfi: Vec<(f64, f64)>,

    /// Initial delay correction (in samples, default: 0.0)
    #[arg(long, default_value_t = 0.0)]
    delay: f64,

    /// Initial rate correction (in Hz, default: 0.0)
    #[arg(long, default_value_t = 0.0)]
    rate: f64,

    /// Limit delay search range during calibration [min_samp max_samp]
    #[arg(long, num_args = 2, value_names = ["MIN_SAMP", "MAX_SAMP"], default_values_t = [-INFINITY, INFINITY], value_delimiter = ' ', allow_negative_numbers = true)]
    drange: Vec<f64>,

    /// Limit rate search range during calibration [min_Hz max_Hz]
    #[arg(long, num_args = 2, value_names = ["MIN_HZ", "MAX_HZ"], default_values_t = [-INFINITY, INFINITY], value_delimiter = ' ', allow_negative_numbers = true)]
    rrange: Vec<f64>,

    /// Number of points for delay quadratic fitting (3, 5 or 7, default: 3)
    #[arg(long, default_value_t = 3)]
    dstep: usize,

    /// Number of points for rate quadratic fitting (3, 5 or 7, default: 3)
    #[arg(long, default_value_t = 3)]
    rstep: usize,

    /// Output header information from input file
    #[arg(long)]
    header: bool,

    /// Skip phase calibration and output only initial FFT results
    #[arg(long)]
    quick: bool,

    /// Suppress console output (file output if --output is specified)
    #[arg(long)]
    noconsole: bool,

        /// Output histogram and cumulative distribution function of time-domain fringes to CSV (filename.rayleigh.csv)
    #[arg(long)]
    pub rayleigh: bool,

    // Fields for internal use
    #[arg(skip)]
    pub output_dir_final: Option<PathBuf>,
    #[arg(skip)]
    delay_search_range_specified: bool,
    #[arg(skip)]
    rate_search_range_specified: bool,
    #[arg(skip)]
    initial_delay_specified: bool,
    #[arg(skip)]
    initial_rate_specified: bool,
}

// --rfi オプションのパーサーヘルパー
fn parse_rfi_range(s: &str) -> Result<(f64, f64), String> {
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 2 {
        return Err(format!("Invalid RFI range format: {}. Expected format: start_mhz, end_mhz", s));
    }
    let start = parts[0].parse::<f64>().map_err(|e| format!("Failed to parse RFI start frequency: {}", e))?;
    let end = parts[1].parse::<f64>().map_err(|e| format!("Failed to parse RFI end frequency: {}", e))?;
    if start < 0.0 || end < 0.0 || start > end {
        return Err(format!("Invalid RFI range: {}. Frequencies must be non-negative and start <= end.", s));
    }
    Ok((start, end))
}

// 引数の後処理を行う関数
fn post_process_args(args: &mut Args) {
    // If frequency_mode or quick_mode is specified, set iterations to 0
    if args.frequency || args.quick {
        args.iter = 0;
    }

    // If drange/rrange are not default values, mark as specified
    if args.drange[0] != -INFINITY || args.drange[1] != INFINITY {
        args.delay_search_range_specified = true;
    }
    if args.rrange[0] != -INFINITY || args.rrange[1] != INFINITY {
        args.rate_search_range_specified = true;
    }

    // If delay/rate are not default values, mark as specified
    if args.delay != 0.0 {
        args.initial_delay_specified = true;
    }
    if args.rate != 0.0 {
        args.initial_rate_specified = true;
    }

    // If --output or --Rayleigh is specified, set output directory
    if args.output || args.rayleigh {
        let input_path = PathBuf::from(&args.input);
        if let Some(parent) = input_path.parent() {
            let base_output_dir = parent.to_path_buf();
            let output_path = base_output_dir.join("frinZ").join("frinZsearch");
            args.output_dir_final = Some(output_path);
        } else {
            eprintln!("Error: Could not determine parent directory from input file path.");
        }
    }
}

fn main() {
    let mut args = match Args::try_parse() {
        Ok(args) => args,
        Err(e) => {
            // If no arguments are provided, print help and exit
            if std::env::args().len() <= 1 {
                let mut cmd = Args::command();
                cmd.print_help().expect("Failed to print help");
                std::process::exit(0);
            } else {
                // For other parsing errors, let clap handle it (prints error and usage)
                e.exit();
            }
        }
    };
    post_process_args(&mut args);

    let mut logger = frinz_logger::Logger::new();
    let mut log_filepath: Option<String> = None;

    if let Some(output_dir) = &args.output_dir_final {
        if args.output {
            let input_path = PathBuf::from(&args.input);
            if let Some(file_stem) = input_path.file_stem() {
                let log_filename = format!("{}_frinZsearch_result.txt", file_stem.to_string_lossy());
                log_filepath = Some(output_dir.join(log_filename).to_string_lossy().into_owned());
            }
        }
        // ディレクトリが存在しない場合は作成
        if let Err(e) = fs::create_dir_all(output_dir) {
            logger.log_fmt(format_args!("Error: Failed to create output directory '{}': {}", output_dir.display(), e));
            // Disable file logging if directory creation fails
            args.output = false;
            args.rayleigh = false;
            log_filepath = None;
        }
    }

    if let Err(e) = logger.setup(!args.noconsole, args.output, log_filepath.as_deref()) {
        logger.log_fmt(format_args!("Logger setup error: {}", e));
        return;
    }

    // If header info output mode
    if args.header {
        match frinz_read::read_binary_file_data(&args.input) {
            Ok(file_data) => {
                logger.log_fmt(format_args!("=== Header Information ==="));
                let magic_word = file_data.header.magic_word;
                let header_version = file_data.header.header_version;
                let software_version = file_data.header.software_version;
                let sampling_speed = file_data.header.sampling_speed;
                let observed_sky_freq = file_data.header.observed_sky_freq;
                let fft_point = file_data.header.fft_point;
                let num_sector = file_data.header.num_sector;
                let station1_name = String::from_utf8_lossy(&file_data.header.station1_name).trim_end_matches('\0').to_string();
                let station1_pos_x = file_data.header.station1_pos_x;
                let station1_pos_y = file_data.header.station1_pos_y;
                let station1_pos_z = file_data.header.station1_pos_z;
                let station1_key = String::from_utf8_lossy(&file_data.header.station1_key).trim_end_matches('\0').to_string();
                let station2_name = String::from_utf8_lossy(&file_data.header.station2_name).trim_end_matches('\0').to_string();
                let station2_pos_x = file_data.header.station2_pos_x;
                let station2_pos_y = file_data.header.station2_pos_y;
                let station2_pos_z = file_data.header.station2_pos_z;
                let station2_key = String::from_utf8_lossy(&file_data.header.station2_key).trim_end_matches('\0').to_string();
                let source_name = String::from_utf8_lossy(&file_data.header.source_name).trim_end_matches('\0').to_string();
                let source_ra = file_data.header.source_ra;
                let source_dec = file_data.header.source_dec;

                logger.log_fmt(format_args!("Magic Word: {:#X}", magic_word));
                logger.log_fmt(format_args!("Header Version: {}", header_version));
                logger.log_fmt(format_args!("Software Version: {}", software_version));
                logger.log_fmt(format_args!("Sampling Speed: {} Hz", sampling_speed));
                logger.log_fmt(format_args!("Observed Sky Freq: {} Hz", observed_sky_freq));
                logger.log_fmt(format_args!("FFT Point: {}", fft_point));
                logger.log_fmt(format_args!("Number of Sectors in file: {}", num_sector));
                logger.log_fmt(format_args!("Station1 Name: {}", station1_name));
                logger.log_fmt(format_args!("Station1 Position: ({}, {}, {})", station1_pos_x, station1_pos_y, station1_pos_z));
                logger.log_fmt(format_args!("Station1 Key: {}", station1_key));
                logger.log_fmt(format_args!("Station2 Name: {}", station2_name));
                logger.log_fmt(format_args!("Station2 Position: ({}, {}, {})", station2_pos_x, station2_pos_y, station2_pos_z));
                logger.log_fmt(format_args!("Station2 Key: {}", station2_key));
                logger.log_fmt(format_args!("Source Name: {}", source_name));
                logger.log_fmt(format_args!("Source RA: {} rad", source_ra));
                logger.log_fmt(format_args!("Source Dec: {} rad", source_dec));
                logger.log_fmt(format_args!("--------------------"));
            }
            Err(e) => {
                logger.log_fmt(format_args!("File reading error: {}", e));
            }
        }
        return; // Exit here if in header display mode
    }

    // Normal processing
    logger.log_fmt(format_args!("--- Input Parameters ---"));
    logger.log_fmt(format_args!("Input file: {}", args.input));
    logger.log_fmt(format_args!("Number of iterations: {}", args.iter));
    logger.log_fmt(format_args!("--------------------"));

    // File reading process
    match frinz_read::read_binary_file_data(&args.input) {
        Ok(file_data) => {
            logger.log_fmt(format_args!("File reading successful!"));

            let first_effective_integration_length = file_data.first_effective_integration_length as f64;

            // --skip の適用
            let num_sectors_to_skip = if args.skip > 0.0 && first_effective_integration_length > 1e-9 {
                (args.skip / first_effective_integration_length).round() as usize
            } else {
                0
            };

            if num_sectors_to_skip > 0 && num_sectors_to_skip < file_data.spectrum_data.len() {
                logger.log_fmt(format_args!("Skipped {} initial sectors ({}s).", num_sectors_to_skip, args.skip));
            } else if num_sectors_to_skip >= file_data.spectrum_data.len() && !file_data.spectrum_data.is_empty() {
                logger.log_fmt(format_args!("Error: --skip value ({}s) results in skipping all data ({} sectors).", args.skip, num_sectors_to_skip));
                return;
            }

            let data_after_skip: Vec<_> = file_data.spectrum_data.iter().skip(num_sectors_to_skip).cloned().collect();
            let times_after_skip: Vec<_> = file_data.sector_start_times_utc.iter().skip(num_sectors_to_skip).cloned().collect();

            if data_after_skip.is_empty() {
                logger.log_fmt(format_args!("Error: No data remaining after applying --skip."));
                return;
            }

            let total_data_after_skip_size = data_after_skip.len();

            for loop_idx in 0..args.loop_count {
                let (current_segment_num_sectors, segment_start_index) = if args.length > 0.0 && first_effective_integration_length > 1e-9 {
                    let mut num_sectors = (args.length / first_effective_integration_length).round() as usize;
                    if num_sectors == 0 { num_sectors = 1; }
                    let start_index = loop_idx * num_sectors;
                    (num_sectors, start_index)
                } else {
                    let base_segment_size = (total_data_after_skip_size + args.loop_count - 1) / args.loop_count;
                    let start_index = loop_idx * base_segment_size;
                    (base_segment_size, start_index)
                };

                if segment_start_index >= total_data_after_skip_size {
                    logger.log_fmt(format_args!("No more data for loop {}. Ending.", loop_idx + 1));
                    break;
                }

                let segment_end_index = (segment_start_index + current_segment_num_sectors).min(total_data_after_skip_size);
                let actual_segment_num_sectors = segment_end_index - segment_start_index;

                if actual_segment_num_sectors == 0 {
                    logger.log_fmt(format_args!("Segment for loop {} is empty. Skipping.", loop_idx + 1));
                    continue;
                }

                let current_all_spectrums = &data_after_skip[segment_start_index..segment_end_index];

                let n_rows_original = current_all_spectrums.len();
                let n_cols_original = if n_rows_original > 0 { current_all_spectrums[0].len() } else { 0 };

                if n_rows_original == 0 || n_cols_original == 0 {
                    logger.log_fmt(format_args!("Error: Loaded data is empty."));
                    return;
                }

                let n_rows_padded = n_rows_original.next_power_of_two() * 2;
                let n_cols_padded = file_data.header.fft_point as usize;

                // C++版の出力に合わせるための情報
                let input_file_path_fs = PathBuf::from(&args.input);
                let input_basename = input_file_path_fs.file_name().unwrap().to_string_lossy();
                let source_name = String::from_utf8_lossy(&file_data.header.source_name).trim_end_matches('\0').to_string();
                
                logger.log_fmt(format_args!("
--- Processing Loop {}/{} ---", loop_idx + 1, args.loop_count));
                logger.log_fmt(format_args!("Input file: {}", input_basename));
                logger.log_fmt(format_args!("Target: {}", source_name));
                if let Some(start_time) = times_after_skip.get(segment_start_index) {
                    logger.log_fmt(format_args!("Segment Start UTC: {}", start_time));
                }
                if let Some(end_time) = times_after_skip.get(segment_end_index - 1) {
                    logger.log_fmt(format_args!("Segment End UTC (approx): {}", end_time));
                }
                logger.log_fmt(format_args!("Processing {} sectors for this loop.", actual_segment_num_sectors));
                logger.log_fmt(format_args!("Original Data Dimensions: Vertical (Sectors) = {}, Horizontal (Channels) = {}", n_rows_original, n_cols_original));
                logger.log_fmt(format_args!("Padded FFT Dimensions: Vertical (N_rows_padded) = {}, Horizontal (N_cols_fft) = {}", n_rows_padded, n_cols_padded));

                let effective_integration_length = first_effective_integration_length;

                // Initial FFT
                match frinz_fft::perform_2d_fft(
                    current_all_spectrums,
                    n_rows_padded,
                    n_cols_padded,
                    &args.rfi,
                    file_data.header.sampling_speed,
                    file_data.header.fft_point,
                ) {
                    Ok((amplitude, fft_out_complex_data)) => {
                        logger.log_fmt(format_args!("Initial FFT calculation successful!"));

                        match frinz_fft::calculate_fft_peak_parameters(
                            &amplitude,
                            &fft_out_complex_data,
                            effective_integration_length as f32,
                            args.drange[0],
                            args.drange[1],
                            args.delay_search_range_specified,
                            args.rrange[0],
                            args.rrange[1],
                            args.rate_search_range_specified,
                        ) {
                            Ok(initial_peak_params) => {
                                logger.log_fmt(format_args!("Initial FFT - StdDev: {:.5e}, SNR (Max/StdDev): {:.3}", 
                                                            initial_peak_params.stddev_amplitude, initial_peak_params.snr));
                                logger.log_fmt(format_args!("Maximum FFT amplitude (%): {:.8} (within specified range if any) found at", 
                                                            100.0 * initial_peak_params.max_amplitude));
                                logger.log_fmt(format_args!("  Fringe Frequency (Hz): {:.8} (Range: {})", 
                                                            initial_peak_params.physical_rate_hz, 
                                                            if args.rate_search_range_specified { format!("[{:.8},{:.8}]", args.rrange[0], args.rrange[1]) } else { "unrestricted".to_string() }));
                                logger.log_fmt(format_args!("  Delay (Sample): {:.8} (Range: {})", 
                                                            initial_peak_params.physical_delay_samples, 
                                                            if args.delay_search_range_specified { format!("[{:.8},{:.8}]", args.drange[0], args.drange[1]) } else { "unrestricted".to_string() }));
                                
                                // Surrounding points and quadratic fit for initial FFT
                                logger.log_fmt(format_args!("  Surrounding points (Initial FFT Peak at i={}, j={}):", initial_peak_params.max_i, initial_peak_params.max_j));
                                logger.log_fmt(format_args!("    Delay (fixed rate bin {}):", initial_peak_params.max_i));
                                let fit_points = args.dstep; // フィッティングに使う点の数 (奇数)
                                let half_fit_points = (fit_points - 1) / 2;

                                // ディレイ方向のフィッティング
                                let mut x_coords_delay_fit = Vec::new();
                                let mut y_values_delay_fit = Vec::new();
                                for dj in - (half_fit_points as isize) ..= (half_fit_points as isize) {
                                    let current_j = (initial_peak_params.max_j as isize + dj) as usize;
                                    if current_j < n_cols_padded {
                                        let physical_delay = (n_cols_padded as f64 / 2.0) - current_j as f64;
                                        x_coords_delay_fit.push(physical_delay);
                                        y_values_delay_fit.push(amplitude[initial_peak_params.max_i][current_j] as f64);
                                        logger.log_fmt(format_args!("      Delay[j={}] (Phys: {:.8} samp): Amp (%): {:.8}", current_j, physical_delay, 100.0 * amplitude[initial_peak_params.max_i][current_j]));
                                    } else {
                                        logger.log_fmt(format_args!("      Delay[j={}]: Out of bounds for fitting", current_j));
                                    }
                                }

                                let mut refined_delay_from_first_fit_samples = initial_peak_params.physical_delay_samples;
                                if x_coords_delay_fit.len() == fit_points {
                                    match frinz_fitting::fit_quadratic_least_squares(&x_coords_delay_fit, &y_values_delay_fit) {
                                        Ok(fit_result) => {
                                            refined_delay_from_first_fit_samples = fit_result.peak_x as f32;
                                            logger.log_fmt(format_args!("      Quadratic Fit (Delay): Refined Peak Delay (samp) = {:.8}, Max Amp Est (%): {:.8} (Coeffs a={:.8}, b={:.8}, c={:.8})",
                                                                        fit_result.peak_x, 100.0 * (fit_result.a * fit_result.peak_x * fit_result.peak_x + fit_result.b * fit_result.peak_x + fit_result.c), 
                                                                        fit_result.a, fit_result.b, fit_result.c));
                                        }
                                        Err(e) => {
                                            logger.log_fmt(format_args!("      Quadratic Fit (Delay) failed: {}. Using non-fitted residual.", e));
                                        }
                                    }
                                } else {
                                    logger.log_fmt(format_args!("      Quadratic Fit (Delay): Not enough valid/contiguous data points ({} collected) for {}-point fitting.", x_coords_delay_fit.len(), fit_points));
                                }

                                // Rate direction fitting
                                logger.log_fmt(format_args!("    Rate (fixed delay bin {}):", initial_peak_params.max_j));
                                let fit_points = args.rstep;
                                let half_fit_points = (fit_points - 1) / 2;
                                let mut x_coords_rate_fit = Vec::new();
                                let mut y_values_rate_fit = Vec::new();
                                for di in - (half_fit_points as isize) ..= (half_fit_points as isize) {
                                    let current_i = (initial_peak_params.max_i as isize + di) as usize;
                                    if current_i < n_rows_padded {
                                        let physical_rate = if effective_integration_length > 1e-9 && n_rows_padded > 0 {
                                            let max_fringe_freq = 1.0 / (2.0 * effective_integration_length as f64);
                                            (current_i as f64 - n_rows_padded as f64 / 2.0)
                                                * (2.0 * max_fringe_freq / n_rows_padded as f64)
                                        } else { 0.0 };
                                        x_coords_rate_fit.push(physical_rate);
                                        y_values_rate_fit.push(amplitude[current_i][initial_peak_params.max_j] as f64);
                                        logger.log_fmt(format_args!("      Rate[i={}] (Phys: {:.8} Hz): Amp (%): {:.8}", current_i, physical_rate, 100.0 * amplitude[current_i][initial_peak_params.max_j]));
                                    } else {
                                        logger.log_fmt(format_args!("      Rate[i={}]: Out of bounds for fitting", current_i));
                                    }
                                }

                                let mut refined_rate_from_first_fit_hz = initial_peak_params.physical_rate_hz;
                                if x_coords_rate_fit.len() == fit_points {
                                    let rate_scale_factor = (10.0 * n_rows_padded as f64) * effective_integration_length as f64;
                                    let scaled_x_coords_rate_fit: Vec<f64> = x_coords_rate_fit.iter().map(|&x| x * rate_scale_factor).collect();

                                    match frinz_fitting::fit_quadratic_least_squares(&scaled_x_coords_rate_fit, &y_values_rate_fit) {
                                        Ok(mut fit_result) => {
                                            refined_rate_from_first_fit_hz = (fit_result.peak_x / rate_scale_factor) as f32;
                                            // 係数を元のスケールに戻す
                                            fit_result.a *= rate_scale_factor * rate_scale_factor;
                                            fit_result.b *= rate_scale_factor;
                                            logger.log_fmt(format_args!("      Quadratic Fit (Rate): Refined Peak Rate (Hz) = {:.8}, Max Amp Est (%): {:.8} (Coeffs a={:.8}, b={:.8}, c={:.8})",
                                                                        refined_rate_from_first_fit_hz, 100.0 * (fit_result.a * fit_result.peak_x * fit_result.peak_x + fit_result.b * fit_result.peak_x + fit_result.c), 
                                                                        fit_result.a, fit_result.b, fit_result.c));
                                        }
                                        Err(e) => {
                                            logger.log_fmt(format_args!("      Quadratic Fit (Rate) failed: {}. Using non-fitted residual.", e));
                                        }
                                    }
                                } else {
                                    logger.log_fmt(format_args!("      Quadratic Fit (Rate): Not enough valid/contiguous data points ({} collected) for {}-point fitting.", x_coords_rate_fit.len(), fit_points));
                                }

                                // --delay または --rate が指定されている場合、初期値を上書き
                                let mut accumulated_rate_correction_hz = refined_rate_from_first_fit_hz;
                                let mut accumulated_delay_correction_samples = refined_delay_from_first_fit_samples;

                                if args.initial_rate_specified {
                                    accumulated_rate_correction_hz = args.rate as f32;
                                    logger.log_fmt(format_args!("Overriding initial rate with command-line specified rate: {} Hz", accumulated_rate_correction_hz));
                                }
                                if args.initial_delay_specified {
                                    accumulated_delay_correction_samples = args.delay as f32;
                                    logger.log_fmt(format_args!("Overriding initial delay with command-line specified delay: {} samples", accumulated_delay_correction_samples));
                                }

                                // Skip iterations if in quick mode or frequency mode
                                if args.quick || args.frequency {
                                    logger.log_fmt(format_args!("Skipping iterations due to quick mode or frequency mode."));
                                } else if args.iter > 0 {
                                    logger.log_fmt(format_args!("
--- Starting Iterative Correction ({} iterations) ---", args.iter));

                                    for iter in 0..args.iter {
                                        logger.log_fmt(format_args!("
-- イテレーション {}/{} --", iter + 1, args.iter));
                                        logger.log_fmt(format_args!("Applied Rate: {:.8} Hz, Applied Delay: {:.8} samples",
                                                                    accumulated_rate_correction_hz, accumulated_delay_correction_samples));

                                        // 位相補正を適用
                                        let corrected_data = frinz_fft::apply_phase_correction(
                                            current_all_spectrums, // 現在のセグメントに補正を適用
                                            accumulated_rate_correction_hz,
                                            accumulated_delay_correction_samples,
                                            effective_integration_length as f32,
                                            file_data.header.sampling_speed,
                                            file_data.header.fft_point,
                                        );

                                        // 補正後のデータでFFT
                                        match frinz_fft::perform_2d_fft(
                                            &corrected_data,
                                            n_rows_padded,
                                            n_cols_padded,
                                            &args.rfi,
                                            file_data.header.sampling_speed,
                                            file_data.header.fft_point,
                                        ) {
                                            Ok((iter_amplitude, iter_fft_out_complex_data)) => {
                                                // 残差ピークパラメータを計算
                                                match frinz_fft::calculate_fft_peak_parameters(
                                                    &iter_amplitude,
                                                    &iter_fft_out_complex_data,
                                                    effective_integration_length as f32,
                                                    args.drange[0],
                                                    args.drange[1],
                                                    args.delay_search_range_specified,
                                                    args.rrange[0],
                                                    args.rrange[1],
                                                    args.rate_search_range_specified,
                                                ) {
                                                    Ok(iter_peak_params) => {
                                                        logger.log_fmt(format_args!("  Iter {} FFT - Max Amp (%): {:.8}, SNR: {:.3}, Residual Rate (Hz): {:.8}, Residual Delay (samp): {:.8}",
                                                                                    iter + 1, 100.0 * iter_peak_params.max_amplitude, iter_peak_params.snr, iter_peak_params.physical_rate_hz, iter_peak_params.physical_delay_samples));

                                                        // --- 2次関数フィッティング (残差) ---
                                                        let fit_points = args.dstep; // フィッティングに使う点の数 (奇数)
                                                        let half_fit_points = (fit_points - 1) / 2;

                                                        // ディレイ方向のフィッティング
                                                        let mut x_coords_delay_fit = Vec::new();
                                                        let mut y_values_delay_fit = Vec::new();
                                                        for dj in - (half_fit_points as isize) ..= (half_fit_points as isize) {
                                                            let current_j = (iter_peak_params.max_j as isize + dj) as usize;
                                                            if current_j < n_cols_padded {
                                                                let physical_delay = (n_cols_padded as f64 / 2.0) - current_j as f64;
                                                                x_coords_delay_fit.push(physical_delay);
                                                                y_values_delay_fit.push(iter_amplitude[iter_peak_params.max_i][current_j] as f64);
                                                                logger.log_fmt(format_args!("      Delay[j={}] (PhysRes: {:.8} samp): Amp (%): {:.8}", current_j, physical_delay, 100.0 * iter_amplitude[iter_peak_params.max_i][current_j]));
                                                            }
                                                        }

                                                        let mut refined_residual_delay = iter_peak_params.physical_delay_samples;
                                                        if x_coords_delay_fit.len() == fit_points {
                                                            match frinz_fitting::fit_quadratic_least_squares(&x_coords_delay_fit, &y_values_delay_fit) {
                                                                Ok(fit_result) => {
                                                                    refined_residual_delay = fit_result.peak_x as f32;
                                                                    logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Delay): Refined Residual Delay (samp) = {:.8}", iter + 1, refined_residual_delay));
                                                                }
                                                                Err(e) => {
                                                                    logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Delay) failed: {}. Using non-fitted residual.", iter + 1, e));
                                                                }
                                                            }
                                                        } else {
                                                            logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Delay): Not enough valid points for {}-point fit. Using non-fitted residual.", iter + 1, fit_points));
                                                        }

                                                        // レート方向のフィッティング
                                                        let fit_points = args.rstep;
                                                        let half_fit_points = (fit_points - 1) / 2;
                                                        let mut x_coords_rate_fit = Vec::new();
                                                        let mut y_values_rate_fit = Vec::new();
                                                        for di in - (half_fit_points as isize) ..= (half_fit_points as isize) {
                                                            let current_i = (iter_peak_params.max_i as isize + di) as usize;
                                                            if current_i < n_rows_padded {
                                                                let physical_rate = if effective_integration_length > 1e-9 && n_rows_padded > 0 {
                                                                    let max_fringe_freq = 1.0 / (2.0 * effective_integration_length as f64);
                                                                    (current_i as f64 - n_rows_padded as f64 / 2.0)
                                                                        * (2.0 * max_fringe_freq / n_rows_padded as f64)
                                                                } else { 0.0 };
                                                                x_coords_rate_fit.push(physical_rate);
                                                                y_values_rate_fit.push(iter_amplitude[current_i][iter_peak_params.max_j] as f64);
                                                                logger.log_fmt(format_args!("      Rate[i={}] (PhysRes: {:.8} Hz): Amp (%): {:.8}", current_i, physical_rate, 100.0 * iter_amplitude[current_i][iter_peak_params.max_j]));
                                                            }
                                                        }

                                                        let mut refined_residual_rate = iter_peak_params.physical_rate_hz;
                                                        if x_coords_rate_fit.len() == fit_points {
                                                            let rate_scale_factor = (10.0 * n_rows_padded as f64) * effective_integration_length as f64;
                                                            let scaled_x_coords_rate_fit: Vec<f64> = x_coords_rate_fit.iter().map(|&x| x * rate_scale_factor).collect();

                                                            match frinz_fitting::fit_quadratic_least_squares(&scaled_x_coords_rate_fit, &y_values_rate_fit) {
                                                                Ok(mut fit_result) => {
                                                                    refined_residual_rate = (fit_result.peak_x / rate_scale_factor) as f32;
                                                                    // 係数を元のスケールに戻す
                                                                    fit_result.a *= rate_scale_factor * rate_scale_factor;
                                                                    fit_result.b *= rate_scale_factor;
                                                                    logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Rate): Refined Residual Rate (Hz) = {:.8} (Coeffs a={:.8}, b={:.8}, c={:.8})",
                                                                                                iter + 1, refined_residual_rate, fit_result.a, fit_result.b, fit_result.c));
                                                                }
                                                                Err(e) => {
                                                                    logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Rate) failed: {}. Using non-fitted residual.", iter + 1, e));
                                                                }
                                                            }
                                                        } else {
                                                            logger.log_fmt(format_args!("      Quadratic Fit (Iter {} Rate): Not enough valid points for {}-point fit. Using non-fitted residual.", iter + 1, fit_points));
                                                        }

                                                        // 累積補正量を更新
                                                        accumulated_rate_correction_hz += refined_residual_rate;
                                                        accumulated_delay_correction_samples += refined_residual_delay;

                                                        logger.log_fmt(format_args!("  Iter {} - Updated Total Applied Rate: {:.8} Hz, Updated Total Applied Delay: {:.8} samples.",
                                                                                    iter + 1, accumulated_rate_correction_hz, accumulated_delay_correction_samples));

                                                    }
                                                    Err(e) => {
                                                        logger.log_fmt(format_args!("Error calculating peak parameters during iteration: {}", e));
                                                        break; // Exit loop on error
                                                    }
                                                }
                                            }
                                            Err(e) => {
                                                logger.log_fmt(format_args!("Error calculating FFT during iteration: {}", e));
                                                break; // Exit loop on error
                                            }
                                        }
                                    }

                                    logger.log_fmt(format_args!("
--- Iterative Correction Complete ---"));
                                    logger.log_fmt(format_args!("Final Rate Correction: {:.8} Hz", accumulated_rate_correction_hz));
                                                                    logger.log_fmt(format_args!("Final Delay Correction: {:.8} samples", accumulated_delay_correction_samples));
                                }

                                // Rayleigh CSV出力
                                if args.rayleigh {
                                    if let Some(output_dir) = &args.output_dir_final {
                                        if let Err(e) = write_rayleigh_csv(
                                            output_dir,
                                            &input_basename,
                                            &initial_peak_params,
                                            &amplitude, // 初期FFTの振幅データ
                                            &logger,
                                        ) {
                                            logger.log_fmt(format_args!("Failed to create Rayleigh CSV: {}", e));
                                        }
                                    } else {
                                        logger.log_fmt(format_args!("Output directory not set, cannot create Rayleigh CSV."));
                                    }
                                }

                            }
                            Err(e) => {
                                logger.log_fmt(format_args!("Error calculating initial peak parameters: {}", e));
                            }
                        }
                    }
                    Err(e) => {
                        logger.log_fmt(format_args!("Initial FFT calculation error: {}", e));
                    }
                }
            }
        }
        Err(e) => {
            logger.log_fmt(format_args!("File reading error: {}", e));
        }
    }
}

// Function to output Rayleigh distribution histogram and CDF to CSV file
fn write_rayleigh_csv(
    output_dir: &PathBuf,
    input_basename: &str,
    peak_params: &frinz_fft::FftPeakParameters,
    amplitude_data: &Vec<Vec<f32>>,
    logger: &frinz_logger::Logger,
) -> Result<(), frinz_error::FrinZError> {
    use std::fs::File;
    use std::io::Write;

    // Construct file path
    let output_filename = format!("{}.rayleigh.csv", input_basename);
    let output_filepath = output_dir.join(output_filename);

    let mut file = File::create(&output_filepath)
        .map_err(|e| frinz_error::FrinZError::Io(e))?;

    // Calculate histogram
    let mut all_amplitudes: Vec<f32> = amplitude_data.iter().flatten().cloned().collect();
    if all_amplitudes.is_empty() {
        return Err(frinz_error::FrinZError::Logic("Amplitude data is empty, cannot create Rayleigh CSV.".to_string()));
    }
    all_amplitudes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let min_amplitude = *all_amplitudes.first().unwrap_or(&0.0);
    let max_amplitude = *all_amplitudes.last().unwrap_or(&0.0);

    const NUM_BINS: usize = 1000;
    let bin_width = if max_amplitude > min_amplitude {
        (max_amplitude - min_amplitude) / NUM_BINS as f32
    } else {
        // If all amplitudes are the same, or only one data point
        // Setting bin_width to 0 would cause division by zero, so return an error or a suitable small value
        // Here, histogram output is meaningless, so return an error
        return Err(frinz_error::FrinZError::Logic("Amplitude data range is zero, cannot create Rayleigh CSV.".to_string()));
    };

    let mut bins: Vec<usize> = vec![0; NUM_BINS];
    for &amp in &all_amplitudes {
        if amp >= min_amplitude && amp <= max_amplitude && bin_width > 0.0 {
            let mut bin_index = ((amp - min_amplitude) / bin_width) as usize;
            if bin_index >= NUM_BINS {
                bin_index = NUM_BINS - 1; // Adjustment for max_amplitude exactly at max_amplitude
            }
            bins[bin_index] += 1;
        }
    }

    // Write header information
    writeln!(file, "# Amplitude Histogram and CDF Data").map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Initial FFT Peak Parameters:").map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Max Amplitude: {:.8}", peak_params.max_amplitude).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Peak Delay (samples): {:.8}", peak_params.physical_delay_samples).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Peak Rate (Hz): {:.8}", peak_params.physical_rate_hz).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# SNR (Max/NoiseLevel): {:.3}", peak_params.snr).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Histogram Parameters:").map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Number of Bins: {}", NUM_BINS).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Min Amplitude for Binning: {:.8}", min_amplitude).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Max Amplitude for Binning: {:.8}", max_amplitude).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "# Bin Width: {:.8}", bin_width).map_err(|e| frinz_error::FrinZError::Io(e))?;
    writeln!(file, "Bin,Frequency,CDF").map_err(|e| frinz_error::FrinZError::Io(e))?;

    let total_count = all_amplitudes.len() as f32;
    let mut cumulative_count = 0;

    for i in 0..NUM_BINS {
        let bin_center = min_amplitude + (i as f32 + 0.5) * bin_width;
        let frequency_count = bins[i];
        cumulative_count += frequency_count;
        let cdf = cumulative_count as f32 / total_count;

        writeln!(file, "{:.8},{:.8},{:.8}", bin_center, frequency_count, cdf)
            .map_err(|e| frinz_error::FrinZError::Io(e))?;
    }

    logger.log_fmt(format_args!("Rayleigh CSV file '{}' created.", output_filepath.display()));
    Ok(())
}

