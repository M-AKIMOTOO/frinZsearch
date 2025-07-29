use rustfft::{FftPlanner, num_complex::Complex};
use std::f32::consts::PI;

use crate::frinz_error::FrinZError;

// Corresponds to C++ FftPeakParameters struct
#[derive(Debug, Clone, Copy)]
pub struct FftPeakParameters {
    pub max_amplitude: f32,
    pub max_i: usize, // Row index of max amplitude
    pub max_j: usize, // Column index of max amplitude
    // pub mean_amplitude: f32, // Added
    pub stddev_amplitude: f32, // Added
    pub physical_rate_hz: f32,
    pub physical_delay_samples: f32,
    pub snr: f32,
    pub success: bool,
}

impl Default for FftPeakParameters {
    fn default() -> Self {
        Self {
            max_amplitude: -1.0,
            max_i: 0,
            max_j: 0,
            // mean_amplitude: 0.0, // Initialization
            stddev_amplitude: 0.0, // Initialization
            physical_rate_hz: 0.0,
            physical_delay_samples: 0.0,
            snr: 0.0,
            success: false,
        }
    }
}


/// Calculates 2D FFT and computes shifted amplitude and post-FFT complex data.
pub fn perform_2d_fft(
    input_data: &[Vec<Complex<f32>>],
    n_rows_padded: usize, // Vertical FFT size (after padding)
    n_cols_padded: usize, // Horizontal FFT size (after padding)
    rfi_ranges_mhz: &[(f64, f64)], // RFI range (MHz)
    sampling_speed: u32, // Sampling speed (Hz)
    fft_point: u32, // Number of FFT points
) -> Result<(Vec<Vec<f32>>, Vec<Vec<Complex<f32>>>), FrinZError> {

    if input_data.is_empty() || input_data[0].is_empty() {
        return Err(FrinZError::Logic("Input data is empty.".to_string()));
    }

    let n_rows_original = input_data.len();
    let n_cols_original = input_data[0].len();

    // --- Data preparation for FFT ---
    // 1. Create zero-padded 2D buffer
    let mut buffer_2d: Vec<Vec<Complex<f32>>> = vec![vec![Complex::new(0.0, 0.0); n_cols_padded]; n_rows_padded];

    // Calculate RFI index ranges
    let mut rfi_index_ranges: Vec<(usize, usize)> = Vec::new();
    if !rfi_ranges_mhz.is_empty() && sampling_speed > 0 && fft_point > 0 {
        let total_bandwidth_hz = sampling_speed as f64;
        let nyquist_bandwidth_hz = total_bandwidth_hz / 2.0;
        let num_channels_half_fft = fft_point / 2;
        let channel_resolution_hz = nyquist_bandwidth_hz / num_channels_half_fft as f64;

        if channel_resolution_hz > f64::EPSILON {
            for &(start_mhz, end_mhz) in rfi_ranges_mhz {
                let start_freq_hz = start_mhz * 1.0e6;
                let end_freq_hz = end_mhz * 1.0e6;

                let mut start_idx = (start_freq_hz / channel_resolution_hz).floor() as usize;
                let mut end_idx = (end_freq_hz / channel_resolution_hz).ceil() as usize - 1;

                start_idx = start_idx.max(0);
                end_idx = end_idx.min(num_channels_half_fft as usize - 1);

                if start_idx <= end_idx {
                    rfi_index_ranges.push((start_idx, end_idx));
                }
            }
        }
    }

    // 2. Copy original data to buffer, applying RFI mask and DC component cut
    for r in 0..n_rows_original {
        for c in 0..n_cols_original {
            let mut is_rfi_channel = false;
            for &(rfi_start, rfi_end) in &rfi_index_ranges {
                if c >= rfi_start && c <= rfi_end {
                    is_rfi_channel = true;
                    break;
                }
            }

            if is_rfi_channel || c == 0 { // c == 0 is DC component
                buffer_2d[r][c] = Complex::new(0.0, 0.0);
            } else {
                buffer_2d[r][c] = input_data[r][c];
            }
        }
    }

    let mut planner = FftPlanner::new();

    // --- Vertical (time axis) FFT ---
    // Apply FFT to each column
    let fft_vertical = planner.plan_fft_forward(n_rows_padded);
    for c in 0..n_cols_padded {
        // Collect data for column c into a temporary buffer
        let mut col_buffer: Vec<Complex<f32>> = (0..n_rows_padded).map(|r| buffer_2d[r][c]).collect();
        // Execute FFT
        fft_vertical.process(&mut col_buffer);
        // Write results back to original buffer
        for r in 0..n_rows_padded {
            buffer_2d[r][c] = col_buffer[r];
        }
    }

    // --- Horizontal (frequency axis) Inverse FFT ---
    let fft_horizontal_inv = planner.plan_fft_inverse(n_cols_padded);
    // Apply Inverse FFT to each row
    for r in 0..n_rows_padded {
        fft_horizontal_inv.process(&mut buffer_2d[r]);
    }

    // Additional scaling to match overall normalization of C++ version
    // C++ final scaling corresponds to 1 / N_rows_original_data
    // rustfft inverse FFT is already scaled by 1 / n_cols_padded
    // Required additional scaling is (1 / n_rows_original) / (1 / n_cols_padded) = n_cols_padded / n_rows_original
    let scaling_factor_to_match_cpp = n_cols_padded as f32 / n_rows_original as f32;
    for r in 0..n_rows_padded {
        for c in 0..n_cols_padded {
            buffer_2d[r][c] *= scaling_factor_to_match_cpp;
        }
    }

    // Clone post-FFT complex data for SNR calculation
    let fft_out_complex_data = buffer_2d.clone();

    // --- fftshift and amplitude calculation ---
    let mut shifted_amplitude = vec![vec![0.0f32; n_cols_padded]; n_rows_padded];
    let mid_rows = n_rows_padded / 2;
    let mid_cols = n_cols_padded / 2;

    for r in 0..n_rows_padded {
        for c in 0..n_cols_padded {
            let shifted_r = (r + mid_rows) % n_rows_padded;
            let shifted_c = (c + mid_cols) % n_cols_padded;

            // Inverse FFT scaling (rustfft does not automatically scale)
            let scaled_val = buffer_2d[r][c] / (n_cols_padded as f32);
            shifted_amplitude[shifted_r][shifted_c] = scaled_val.norm_sqr().sqrt(); // .norm_sqr().sqrt() は .norm() と同じ
        }
    }

    Ok((shifted_amplitude, fft_out_complex_data))
}


/// FFT後の振幅データからピークパラメータを計算する
pub fn calculate_fft_peak_parameters(
    fft_shifted_amplitude: &Vec<Vec<f32>>,
    fft_out_complex_data: &Vec<Vec<Complex<f32>>>,
    effective_integration_length: f32, // Required for rate calculation
    delay_search_min: f64,
    delay_search_max: f64,
    delay_search_range_specified: bool,
    rate_search_min: f64,
    rate_search_max: f64,
    rate_search_range_specified: bool,
) -> Result<FftPeakParameters, FrinZError> {

    if fft_shifted_amplitude.is_empty() || fft_shifted_amplitude[0].is_empty() {
        return Err(FrinZError::Logic("Amplitude data is empty.".to_string()));
    }

    let n_rows_padded = fft_shifted_amplitude.len();
    let n_cols_padded = fft_shifted_amplitude[0].len();

    let mut params = FftPeakParameters::default();
    let mut peak_found_in_range = false;

    // --- ピーク探索 (探索範囲フィルタリングあり) ---
    for r in 0..n_rows_padded {
        for c in 0..n_cols_padded {
            // 現在のピクセル (i, j) の物理レートとディレイを計算
            let current_physical_rate_hz = if effective_integration_length > 1e-9 && n_rows_padded > 0 {
                let max_fringe_freq = 1.0 / (2.0 * effective_integration_length);
                (r as f32 - n_rows_padded as f32 / 2.0)
                    * (2.0 * max_fringe_freq / n_rows_padded as f32)
            } else { 0.0 };

            let current_physical_delay_samples = if n_cols_padded > 0 {
                (n_cols_padded as f32 / 2.0) - c as f32
            } else { 0.0 };

            // 探索範囲内にあるかチェック
            let mut in_delay_range = true;
            if delay_search_range_specified {
                if current_physical_delay_samples < delay_search_min as f32 || current_physical_delay_samples > delay_search_max as f32 {
                    in_delay_range = false;
                }
            }

            let mut in_rate_range = true;
            if rate_search_range_specified {
                if current_physical_rate_hz < rate_search_min as f32 || current_physical_rate_hz > rate_search_max as f32 {
                    in_rate_range = false;
                }
            }

            if in_delay_range && in_rate_range {
                if fft_shifted_amplitude[r][c] > params.max_amplitude {
                    params.max_amplitude = fft_shifted_amplitude[r][c];
                    params.max_i = r;
                    params.max_j = c;
                    params.physical_rate_hz = current_physical_rate_hz;
                    params.physical_delay_samples = current_physical_delay_samples;
                    peak_found_in_range = true;
                }
            }
        }
    }

    if !peak_found_in_range {
        params.success = false;
        return Err(FrinZError::Logic("No peak found within the specified search ranges.".to_string()));
    }

    // --- SNR calculation (matching C++ noise_level_for_snr) ---
    let mut complex_sum_real = 0.0f64;
    let mut complex_sum_imag = 0.0f64;
    let mut complex_data_count = 0;

    for r in 0..n_rows_padded {
        for c in 0..n_cols_padded {
            // C++ version excludes DC component (0,0)
            //if r == 0 && c == 0 {
            //    continue;
            //}
            let val = fft_out_complex_data[r][c];
            complex_sum_real += val.re as f64;
            complex_sum_imag += val.im as f64;
            complex_data_count += 1;
        }
    }

    let mut noise_level_for_snr = 0.0f32;
    if complex_data_count > 0 {
        let mean_real = (complex_sum_real / complex_data_count as f64) as f32;
        let mean_imag = (complex_sum_imag / complex_data_count as f64) as f32;
        let mean_complex = Complex::new(mean_real, mean_imag);

        let mut sum_abs_diff = 0.0f64;
        for r in 0..n_rows_padded {
            for c in 0..n_cols_padded {
                let val = fft_out_complex_data[r][c];
                sum_abs_diff += (val - mean_complex).norm() as f64; // norm() is magnitude
            }
        }
        noise_level_for_snr = (sum_abs_diff / complex_data_count as f64 / n_cols_padded as f64) as f32;
    }

    // params.stddev_amplitude stores the calculated noise_level.
    params.stddev_amplitude = noise_level_for_snr;
    // Note: This variable originally stored standard deviation, but now stores noise level due to SNR calculation method change.
    // params.mean_amplitude does not seem to have a direct equivalent in the C++ output, so it is not set here or remains 0.0.
    // In the C++ version of `calculate_fft_peak_parameters` where `mean_amplitude` comes from `calculate_mean_stddev_from_float_vec`,
    // which is the mean of `fft_shifted_amplitude`, the current Rust code's `mean_amplitude` calculation is correct.
    // However, the user is requesting `StdDev` matching, so `stddev_amplitude` calculation is prioritized.
    // `mean_amplitude` is left as is.

    if params.stddev_amplitude > f32::EPSILON { // ゼロ除算を避ける
        params.snr = params.max_amplitude / params.stddev_amplitude;
    } else if params.max_amplitude.abs() < f32::EPSILON { // 最大振幅もほぼ0の場合
        params.snr = 0.0;
    } else { // ノイズがほぼ0で信号がある場合
        params.snr = f32::INFINITY;
    }

    params.success = true;
    Ok(params)
}

/// Applies phase correction to input data
pub fn apply_phase_correction(
    input_data: &[Vec<Complex<f32>>],
    rate_hz_for_correction: f32,
    delay_samples_for_correction: f32,
    effective_integration_length: f32,
    sampling_speed: u32,
    fft_point: u32,
) -> Vec<Vec<Complex<f32>>> {
    let mut corrected_data = input_data.to_vec();

    let n_rows_original = input_data.len();
    let n_cols_original = if n_rows_original > 0 { input_data[0].len() } else { 0 };

    let can_phase_correct = sampling_speed > 0 && fft_point >= 2 && effective_integration_length.abs() > 1e-9 && n_cols_original > 0;

    if can_phase_correct {
        let py_equiv_sampling_speed_mhz = sampling_speed as f32 / 1.0e6;
        let stop_val_for_linspace_mhz = (py_equiv_sampling_speed_mhz / 2.0).floor() - 1.0;

        for r_orig in 0..n_rows_original {
            let time_for_rate_corr_sec = r_orig as f32 * effective_integration_length;
            let rate_corr_factor = Complex::new(0.0, -2.0 * PI * rate_hz_for_correction * time_for_rate_corr_sec).exp();

            for c_orig in 0..n_cols_original {
                let mut original_val = corrected_data[r_orig][c_orig];

                // Delay correction factor
                let freq_k_hz_for_delay_corr = if n_cols_original > 1 {
                    (c_orig as f32 * stop_val_for_linspace_mhz as f32 / (n_cols_original - 1) as f32) * 1.0e6
                } else { 0.0 };
                let delay_seconds = delay_samples_for_correction / sampling_speed as f32;
                let delay_corr_factor = Complex::new(0.0, -2.0 * PI * delay_seconds * freq_k_hz_for_delay_corr).exp();
                
                original_val *= rate_corr_factor * delay_corr_factor;
                corrected_data[r_orig][c_orig] = original_val;
            }
        }
    }
    corrected_data
}
/*
#[cfg(test)]
mod tests {
    use super::*;
    use rustfft::num_complex::Complex;

    #[test]
    fn test_perform_2d_fft_basic() {
        // 2x2のシンプルなデータ
        let input_data = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)],
            vec![Complex::new(3.0, 0.0), Complex::new(4.0, 0.0)],
        ];
        let n_rows_padded = 2;
        let n_cols_padded = 2;
        let rfi_ranges_mhz = vec![];
        let sampling_speed = 1000;
        let fft_point = 2;

        let result = perform_2d_fft(
            &input_data,
            n_rows_padded,
            n_cols_padded,
            &rfi_ranges_mhz,
            sampling_speed,
            fft_point,
        );
        assert!(result.is_ok());
        let (amplitude, _) = result.unwrap();

        // Basic check for FFT result shape
        assert_eq!(amplitude.len(), n_rows_padded);
        assert_eq!(amplitude[0].len(), n_cols_padded);

        // Check that a peak exists (specific value depends on FFT properties)
        let mut max_amp = 0.0f32;
        for r in 0..n_rows_padded {
            for c in 0..n_cols_padded {
                if amplitude[r][c] > max_amp {
                    max_amp = amplitude[r][c];
                }
            }
        }
        assert!(max_amp > 0.0, "Max amplitude should be greater than 0");
    }

    #[test]
    fn test_perform_2d_fft_rfi_mask() {
        // 2x2のシンプルなデータ
        let input_data = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
        ];
        let n_rows_padded = 2;
        let n_cols_padded = 2;
        // 0-10MHzをRFIとして指定 (channel 0がRFIになるはず)
        let rfi_ranges_mhz = vec![(0.0, 10.0)]; 
        let sampling_speed = 20; // Nyquist = 10Hz, num_channels_half_fft = 1 (fft_point/2)
        let fft_point = 2; // channel 0, 1

        // perform_2d_fft creates buffer_2d internally and applies RFI mask
        let result = perform_2d_fft(
            &input_data,
            n_rows_padded,
            n_cols_padded,
            &rfi_ranges_mhz,
            sampling_speed,
            fft_point,
        );
        assert!(result.is_ok());
        let (_, fft_out_complex_data) = result.unwrap();

        // Verify that channels with RFI mask applied (channel 0) are zero
        // DC component should also be zero
        // Although it should check the corresponding part of buffer_2d before FFT,
        // perform_2d_fft does not return buffer_2d, so check indirectly from FFT results.
        // However, there is no guarantee that FFT results will be exactly 0,
        // so modify the test here.
        // Verify that FFT results for channels with RFI mask applied differ from those without RFI.
        // Or, compare FFT results between data before and after RFI mask application.

        // Here, verify that FFT results for channels with RFI mask applied
        // are at least different from values expected from original input data
        // (expect very small values, not strict 0)
        assert!(fft_out_complex_data[0][0].norm() < 1e-5);
        assert!(fft_out_complex_data[1][0].norm() < 1e-5);
    }

    #[test]
    fn test_calculate_fft_peak_parameters_basic() {
        let fft_shifted_amplitude = vec![
            vec![0.1, 0.2, 0.3],
            vec![0.4, 0.9, 0.5],
            vec![0.6, 0.7, 0.8],
        ];
        let fft_out_complex_data = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
        ];
        let effective_integration_length = 1.0;
        let delay_search_min = -100.0;
        let delay_search_max = 100.0;
        let delay_search_range_specified = false;
        let rate_search_min = -100.0;
        let rate_search_max = 100.0;
        let rate_search_range_specified = false;

        let result = calculate_fft_peak_parameters(
            &fft_shifted_amplitude,
            &fft_out_complex_data,
            effective_integration_length,
            delay_search_min,
            delay_search_max,
            delay_search_range_specified,
            rate_search_min,
            rate_search_max,
            rate_search_range_specified,
        );
        assert!(result.is_ok());
        let params = result.unwrap();

        assert_eq!(params.max_amplitude, 0.9);
        assert_eq!(params.max_i, 1);
        assert_eq!(params.max_j, 1);
        // Physical quantity calculation depends on padding size, so strict values are not checked here
        // assert!(params.physical_rate_hz != 0.0);
        // assert!(params.physical_delay_samples != 0.0);
        assert!(params.snr > 0.0);
        assert!(params.success);
    }

    #[test]
    fn test_apply_phase_correction_basic() {
        let input_data = vec![
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
            vec![Complex::new(1.0, 0.0), Complex::new(1.0, 0.0)],
        ];
        let rate_hz_for_correction = 1.0;
        let delay_samples_for_correction = 1.0;
        let effective_integration_length = 1.0;
        let sampling_speed = 10000; // Larger value
        let fft_point = 1024; // Larger value

        let corrected_data = apply_phase_correction(
            &input_data,
            rate_hz_for_correction,
            delay_samples_for_correction,
            effective_integration_length,
            sampling_speed,
            fft_point,
        );

        assert_eq!(corrected_data.len(), 2);
        assert_eq!(corrected_data[0].len(), 2);
        // Verify that phase correction has been applied (values have changed)
        assert_ne!(corrected_data[0][0], input_data[0][0]);
    }
}
*/
