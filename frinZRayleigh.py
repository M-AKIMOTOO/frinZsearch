#!/usr/bin/env python3 
# Made by Gemini
# 2025/06/05
#

from fitter import Fitter
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rayleigh
from scipy.optimize import curve_fit
import argparse
import os


#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["xtick.direction"]     = "in"       
plt.rcParams["ytick.direction"]     = "in"       
plt.rcParams["xtick.minor.visible"] = True       
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.top"]           = True
plt.rcParams["xtick.bottom"]        = True
plt.rcParams["ytick.left"]          = True
plt.rcParams["ytick.right"]         = True  
plt.rcParams["xtick.major.size"]    = 5          
plt.rcParams["ytick.major.size"]    = 5          
plt.rcParams["xtick.minor.size"]    = 3          
plt.rcParams["ytick.minor.size"]    = 3          
plt.rcParams["axes.grid"]           = False
plt.rcParams["grid.color"]          = "lightgray"
plt.rcParams["axes.labelsize"]      = 15
plt.rcParams["font.size"]           = 12



def parse_rayleigh_csv(filepath):
    """
    Parses the Rayleigh CSV file, extracting comments and data.
    Returns: df, snr, min_amp_data, actual_max_peak_amplitude, comments
    """
    comments = {}
    data_lines = []
    header_cols = None

    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    if ':' in line:
                        try:
                            key, value = line[1:].split(':', 1)
                            comments[key.strip()] = value.strip()
                        except ValueError:
                            # Handle comment lines without a colon gracefully
                            pass 
                elif header_cols is None:
                    header_cols = [col.strip() for col in line.split(',')]
                else:
                    data_lines.append([val.strip() for val in line.split(',')])
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found at: {filepath}")
    except Exception as e:
        raise IOError(f"Error reading CSV file: {e}")

    if not header_cols:
        raise ValueError("CSV header not found.")
    if not data_lines:
        # It's possible to have a header but no data lines
        print(f"Warning: No data rows found in CSV file: {filepath}")
        df = pd.DataFrame([], columns=header_cols) # Return empty DataFrame with headers
    else:
        df = pd.DataFrame(data_lines, columns=header_cols)
    
    # Convert data columns to numeric
    required_cols = ['AmplitudeBinCenter', 'FrequencyCount']
    if not all(col in df.columns for col in required_cols):
        missing = [col for col in required_cols if col not in df.columns]
        raise ValueError(f"Missing required columns in CSV: {', '.join(missing)}")

    try:
        for col in required_cols:
             df[col] = pd.to_numeric(df[col], errors='coerce') # Use coerce to turn non-numeric into NaN
        # Drop rows where required numeric columns failed to convert
        df.dropna(subset=required_cols, inplace=True)

        if 'CDF' in df.columns:
             df['CDF'] = pd.to_numeric(df['CDF'], errors='coerce')
             df.dropna(subset=['CDF'], inplace=True) # Drop if CDF conversion failed

    except Exception as e:
        # This catch is less likely with errors='coerce' but kept for safety
        raise ValueError(f"Could not convert data columns to numeric. Please check CSV format. Error: {e}")

    if df.empty and data_lines:
         raise ValueError("No valid numeric data rows found after parsing and conversion.")


    # Extract specific parameters from comments or data
    snr = float(comments.get('SNR (Max/NoiseLevel)', '0.0'))
    
    # Min Amplitude for Binning: Prefer comment, fallback to min bin center if data exists
    min_amp_data = 0.0 # Default
    parsed_min_amp = comments.get('Min Amplitude for Binning')
    if parsed_min_amp:
        try:
            min_amp_data = float(parsed_min_amp)
        except ValueError:
            print(f"Warning: Could not parse 'Min Amplitude for Binning' ('{parsed_min_amp}') from comments. Using default 0.0 or min bin center.")
    if df['AmplitudeBinCenter'].min() < min_amp_data: # Use actual min bin center if smaller
         min_amp_data = df['AmplitudeBinCenter'].min() if not df.empty else 0.0

    # Actual Max Peak Amplitude: Prefer comment.
    # If comment value is smaller than data max (e.g. signal removed), use comment.
    # If comment value is larger, or no comment, use data max.
    actual_max_peak_amplitude = 0.0 # Default
    parsed_max_amp = comments.get('Max Amplitude')
    if parsed_max_amp:
        try:
            actual_max_peak_amplitude = float(parsed_max_amp)
        except ValueError: # If parsing fails, use data max
            print(f"Warning: Could not parse 'Max Amplitude' ('{parsed_max_amp}') from comments. Using max amplitude from data.")
            actual_max_peak_amplitude = df['AmplitudeBinCenter'].max() if not df.empty else 0.0
    elif not df.empty : # If no comment, use data max
        actual_max_peak_amplitude = df['AmplitudeBinCenter'].max()


    return df, snr, min_amp_data, actual_max_peak_amplitude, comments

def calculate_weighted_stats(bin_centers, frequencies):
    """
    Calculates weighted mean and standard deviation.
    Assumes bin_centers and frequencies are pandas Series or numpy arrays.
    """
    if bin_centers.empty or frequencies.empty or frequencies.sum() <= 0:
        # print("Warning: Empty data or all frequencies are zero, cannot calculate weighted stats.")
        return 0.0, 0.0
    
    # Ensure inputs are numpy arrays for consistent behavior
    bin_centers_arr = bin_centers.to_numpy() if isinstance(bin_centers, pd.Series) else np.array(bin_centers)
    frequencies_arr = frequencies.to_numpy() if isinstance(frequencies, pd.Series) else np.array(frequencies)

    # Filter out zero frequencies to avoid issues with np.average weights
    positive_freq_indices = frequencies_arr > 0
    if not np.any(positive_freq_indices):
         # print("Warning: All frequencies are zero, cannot calculate weighted stats.")
         return 0.0, 0.0

    filtered_bin_centers = bin_centers_arr[positive_freq_indices]
    filtered_frequencies = frequencies_arr[positive_freq_indices]

    weighted_mean = np.average(filtered_bin_centers, weights=filtered_frequencies)
    weighted_variance = np.average((filtered_bin_centers - weighted_mean)**2, weights=filtered_frequencies)
    if weighted_variance < 0: # Handle potential floating point inaccuracies
        weighted_variance = 0
    weighted_std = np.sqrt(weighted_variance)
    return weighted_mean, weighted_std

def rayleigh_pdf_func(x, sigma):
    """
    Rayleigh probability density function.
    x: array-like (must be non-negative)
    sigma: scale parameter (sigma > 0)
    """
    x = np.asarray(x) # Ensure x is a numpy array
    # Handle x < 0 explicitly, PDF is 0 there
    pdf_values = np.zeros_like(x, dtype=float)
    positive_x_indices = x > 1e-9 # Use a small epsilon to handle floating point zeros
    
    if sigma <= 1e-9:
        # print(f"Warning: sigma is zero or negative ({sigma}), returning zero PDF.")
        return pdf_values # Return zeros if sigma is invalid

    # Calculate PDF only for positive x values
    x_positive = x[positive_x_indices]
    pdf_values[positive_x_indices] = (x_positive / sigma**2) * np.exp(-x_positive**2 / (2 * sigma**2))
    
    return pdf_values


def rebin_data(original_bin_centers, original_frequencies, original_bin_width, rebin_factor):
    """
    Rebins histogram data.
    original_bin_centers: pandas Series or numpy array of original bin centers.
    original_frequencies: pandas Series or numpy array of original frequencies.
    original_bin_width: float, the width of a single original bin.
    rebin_factor: int, the number of original bins to combine into one new bin.
    Returns: (new_bin_centers_series, new_frequencies_series, new_bin_width)
    """
    if not isinstance(rebin_factor, int) or rebin_factor <= 0:
        raise ValueError("rebin_factor must be a positive integer.")
    if rebin_factor == 1: # No rebinning needed
        return original_bin_centers, original_frequencies, original_bin_width

    n_original_bins = len(original_frequencies)
    if n_original_bins == 0:
        return pd.Series([], dtype=float), pd.Series([], dtype=float), original_bin_width * rebin_factor

    if rebin_factor > n_original_bins:
        print(f"Warning: rebin_factor ({rebin_factor}) is larger than the number of original bins ({n_original_bins}). Returning original data effectively.")
        return original_bin_centers, original_frequencies, original_bin_width


    n_new_bins = n_original_bins // rebin_factor
    # Handle case where original bins are not perfectly divisible by rebin_factor
    # We will drop the last few bins if they don't form a full new bin
    if n_original_bins % rebin_factor != 0:
         print(f"Warning: Number of original bins ({n_original_bins}) is not perfectly divisible by rebin_factor ({rebin_factor}). Dropping last {n_original_bins % rebin_factor} original bins.")

    new_frequencies_arr = np.zeros(n_new_bins)
    new_bin_centers_arr = np.zeros(n_new_bins)
    
    original_bin_centers_arr = original_bin_centers.to_numpy() if isinstance(original_bin_centers, pd.Series) else np.array(original_bin_centers)
    original_frequencies_arr = original_frequencies.to_numpy() if isinstance(original_frequencies, pd.Series) else np.array(original_frequencies)

    for i in range(n_new_bins):
        start_idx = i * rebin_factor
        end_idx = start_idx + rebin_factor

        new_frequencies_arr[i] = np.sum(original_frequencies_arr[start_idx:end_idx])
        
        # Calculate the center of the new, larger bin
        # The center is the average of the original bin centers within the new bin
        new_bin_centers_arr[i] = np.mean(original_bin_centers_arr[start_idx:end_idx])
        
        # Alternative center calculation based on edges (should be equivalent if original bins are uniform)
        # left_edge_of_first_original_bin = original_bin_centers_arr[start_idx] - (original_bin_width / 2.0)
        # new_large_bin_width = original_bin_width * rebin_factor
        # new_bin_centers_arr[i] = left_edge_of_first_original_bin + (new_large_bin_width / 2.0)


    return pd.Series(new_bin_centers_arr), pd.Series(new_frequencies_arr), original_bin_width * rebin_factor

def fit_rayleigh_to_data(bin_centers, frequencies, bin_width, overall_stats, fit_max_amp_arg=None):
    """
    Fits a Rayleigh distribution to the specified data range.
    Uses Fitter first, falls back to curve_fit.
    Returns the fitted sigma or None if fitting fails.
    """
    mean_val = overall_stats['mean']
    std_val = overall_stats['std']
    actual_max_peak_val = overall_stats['actual_max_peak']

    # Select data for fitting: positive bin centers, potentially limited by fit_max_amp
    # Define the upper limit for fitting based on 3-sigma of the *overall* data,
    # unless a specific fit_max_amp is provided and is lower.
    fit_upper_limit_from_stats = mean_val + 3 * std_val if std_val > 1e-9 else np.inf

    effective_fit_max_amp = fit_upper_limit_from_stats
    if fit_max_amp_arg is not None and fit_max_amp_arg > 0:
        effective_fit_max_amp = min(fit_upper_limit_from_stats, fit_max_amp_arg)
        print(f"Info: Rayleigh fitting will be limited to amplitudes <= {effective_fit_max_amp:.4g} (user-provided or 3-sigma).")
    elif fit_max_amp_arg is not None and fit_max_amp_arg <= 0:
         print(f"Warning: Invalid --fit-max-amp value ({fit_max_amp_arg}). Ignoring and using 3-sigma limit.")


    fit_indices = bin_centers > 1e-9 # Rayleigh is for x > 0
    if effective_fit_max_amp != np.inf :
        fit_indices = fit_indices & (bin_centers <= effective_fit_max_amp)

    fit_bin_centers = bin_centers[fit_indices]
    fit_frequencies = frequencies[fit_indices]

    fitted_sigma = None # Initialize

    if bin_width > 1e-9 and fit_frequencies.sum() > 1e-9:
        # Normalize frequencies of the *fitting data* to represent probability density
        normalized_fit_frequencies = fit_frequencies / (fit_frequencies.sum() * bin_width)

        if len(fit_bin_centers) > 1 and not fit_bin_centers.empty:
            # Attempt 1: Use Fitter library with reconstructed data
            reconstructed_data_for_fitter = []
            try:
                if not fit_bin_centers.empty and not fit_frequencies.empty and fit_frequencies.sum() > 0 and fit_bin_centers.nunique() > 0 :
                    total_points_in_fit_range = fit_frequencies.sum()
                    # Manageable threshold for reconstruction (e.g., 10 million points)
                    if total_points_in_fit_range > 0 and total_points_in_fit_range < 10_000_000: 
                        # Round frequencies before converting to int for np.repeat
                        int_fit_frequencies = fit_frequencies.round().astype(int) 
                        if (int_fit_frequencies < 0).any():
                            print("Warning: Negative frequencies found, cannot reconstruct data for Fitter.")
                        else:
                            reconstructed_data_for_fitter = np.repeat(
                                fit_bin_centers.to_numpy(),
                                int_fit_frequencies.to_numpy()
                            )
                
                if len(reconstructed_data_for_fitter) > 1:
                    print("Info: Attempting to fit Rayleigh distribution using Fitter library...")
                    # Fix location to 0 for Rayleigh. Use method='mle' explicitly.
                    fitter_params_config = {'rayleigh': {'floc': 0, 'method': 'mle'}} 
                    # Increase timeout if needed for large datasets
                    f = Fitter(reconstructed_data_for_fitter, distributions=['rayleigh'], timeout=120, fitter_params=fitter_params_config)
                    f.fit()
                    
                    if 'rayleigh' in f.fitted_param and f.fitted_param['rayleigh'] is not None:
                        params_fitter = f.fitted_param['rayleigh']
                        # For rayleigh.fit with floc=0, it returns (scale,).
                        if len(params_fitter) == 1:
                            temp_sigma = params_fitter[0]
                        # Handle case where fitter might return (loc, scale) even with floc=0 if loc is close to 0
                        elif len(params_fitter) == 2 and abs(params_fitter[0]) < 1e-6 : 
                            temp_sigma = params_fitter[1]
                        else:
                            print(f"Warning: Fitter returned unexpected parameters for Rayleigh: {params_fitter}. Will try curve_fit.")
                            temp_sigma = None

                        # Sanity check on the fitted sigma
                        # It should be positive and not excessively large compared to the data range
                        max_fit_amp = fit_bin_centers.max() if not fit_bin_centers.empty else actual_max_peak_val
                        upper_bound_check = max(max_fit_amp * 2, actual_max_peak_val * 1.5) if max_fit_amp > 1e-9 else actual_max_peak_val * 1.5
                        
                        if temp_sigma is not None and temp_sigma > 1e-9 and temp_sigma < upper_bound_check:
                            fitted_sigma = temp_sigma
                            print(f"Info: Rayleigh sigma from Fitter: {fitted_sigma:.4g}")
                        elif temp_sigma is not None:
                            print(f"Warning: Fitter resulted in potentially invalid sigma ({temp_sigma:.4g}). Will try curve_fit.")
                    else:
                        print("Warning: Fitter did not successfully fit Rayleigh or returned None params. Will try curve_fit.")
                elif not fit_bin_centers.empty : # Only print if there was data to begin with
                    print("Info: Not enough reconstructed data points for Fitter. Will use curve_fit.")
            
            except ImportError:
                print("Warning: Fitter library not installed (pip install fitter). Falling back to curve_fit.")
            except Exception as e_fitter:
                print(f"Warning: Fitter library fitting failed: {e_fitter}. Falling back to curve_fit.")

            # Attempt 2: Fallback to curve_fit if Fitter failed or was not used
            if fitted_sigma is None:
                print("Info: Using curve_fit for Rayleigh distribution.")
                # --- Initial guess for sigma (for curve_fit) ---
                initial_sigma_guess = 1.0  # Default fallback
                
                # Calculate stats *only* for the data used for fitting
                mean_fit, std_fit = calculate_weighted_stats(fit_bin_centers, fit_frequencies)

                # 1. Try sigma from mode of fitting data (Rayleigh mode = sigma)
                if not fit_frequencies.empty and fit_frequencies.sum() > 0:
                    # Find the bin center corresponding to the highest frequency in the fitting range
                    mode_index = fit_frequencies.idxmax() 
                    sigma_from_mode = fit_bin_centers.loc[mode_index] if isinstance(fit_bin_centers, pd.Series) else fit_bin_centers[mode_index]
                    if sigma_from_mode > 1e-9: initial_sigma_guess = sigma_from_mode
                
                # 2. Fallback to sigma from mean of fitting data if mode-based was not good
                if initial_sigma_guess <= 1e-9 and mean_fit > 1e-9: # Check if mode-based was default or invalid
                    sigma_from_mean_fit_data = mean_fit / np.sqrt(np.pi / 2)
                    if sigma_from_mean_fit_data > 1e-9: initial_sigma_guess = sigma_from_mean_fit_data
                
                # 3. Fallback to sigma from std of fitting data if others were not good
                if initial_sigma_guess <= 1e-9 and std_fit > 1e-9: # Check if previous were default or invalid
                    denominator_std_fit_data = np.sqrt((4 - np.pi) / 2)
                    if denominator_std_fit_data > 1e-9:
                        sigma_from_std_fit_data = std_fit / denominator_std_fit_data
                        if sigma_from_std_fit_data > 1e-9: initial_sigma_guess = sigma_from_std_fit_data
                
                # 4. Try MLE initial guess for curve_fit's p0 using reconstructed data (if manageable)
                try:
                    if not fit_bin_centers.empty and not fit_frequencies.empty and fit_frequencies.sum() > 0 and fit_bin_centers.nunique() > 0 :
                        total_points_in_fit_range_mle = fit_frequencies.sum()
                        if total_points_in_fit_range_mle > 0 and total_points_in_fit_range_mle < 10_000_000: # Manageable threshold
                            int_fit_frequencies_mle = fit_frequencies.round().astype(int) # Round frequencies before converting to int
                            if not (int_fit_frequencies_mle < 0).any():
                                reconstructed_data_for_mle_init = np.repeat(
                                    fit_bin_centers.to_numpy(),
                                    int_fit_frequencies_mle.to_numpy()
                                )
                                if len(reconstructed_data_for_mle_init) > 1:
                                    # rayleigh.fit with floc=0 returns (scale,)
                                    mle_scale_param_init = rayleigh.fit(reconstructed_data_for_mle_init, floc=0)
                                    mle_sigma_init = mle_scale_param_init[0] if isinstance(mle_scale_param_init, tuple) else mle_scale_param_init
                                    
                                    # Sanity check for MLE initial guess
                                    max_fit_amp = fit_bin_centers.max() if not fit_bin_centers.empty else actual_max_peak_val
                                    upper_bound_check_init = max(max_fit_amp * 2, actual_max_peak_val * 1.5) if max_fit_amp > 1e-9 else actual_max_peak_val * 1.5

                                    if mle_sigma_init > 1e-9 and mle_sigma_init < upper_bound_check_init:
                                        # Only use MLE guess if current guess is default or MLE is significantly different/better
                                        if initial_sigma_guess <= 1e-9 or abs(initial_sigma_guess - mle_sigma_init) / initial_sigma_guess > 0.2: # Heuristic: Use MLE if it's >20% different or current is default
                                             initial_sigma_guess = mle_sigma_init
                                             # print(f"Debug: Using MLE-based initial sigma for curve_fit: {initial_sigma_guess:.4g}")
                                    # else:
                                        # print(f"Debug: MLE initial sigma for curve_fit ({mle_sigma_init:.4g}) out of range or moment-based preferred.")

                except Exception as e_mle_init:
                    # print(f"Debug: MLE for initial guess for curve_fit failed: {e_mle_init}")
                    pass 

                # Final check on initial guess before curve_fit
                if initial_sigma_guess <= 1e-9: initial_sigma_guess = 1.0

                try:
                    # Fit using the normalized frequencies of the *fitting data*
                    params_cv, covariance = curve_fit(rayleigh_pdf_func, fit_bin_centers, normalized_fit_frequencies, p0=[initial_sigma_guess], bounds=(1e-9, np.inf))
                    fitted_sigma = params_cv[0]
                    print(f"Info: Rayleigh sigma from curve_fit: {fitted_sigma:.4g}")
                except (RuntimeError, ValueError) as e_cv:
                    print(f"Warning: curve_fit for Rayleigh fitting also failed. {e_cv}")
                    fitted_sigma = None # Explicitly set

        else:
            print("Warning: Not enough data points in the fitting range to attempt Rayleigh fit.")
    else:
        print("Warning: Cannot perform Rayleigh fit due to zero bin width or zero total frequency in fitting range.")

    return fitted_sigma


def plot_amplitude_data(bin_centers, frequencies, bin_width, overall_stats, fitted_sigma, snr_value, plot_type, plot_title, input_filename, x_limit=None, output_filename_base="plot"):
    """
    Plots the amplitude data (histogram or scatter) and the fitted Rayleigh curve.
    """
    mean_val = overall_stats['mean']
    std_val = overall_stats['std']
    min_data_val = overall_stats['min_data'] 
    actual_max_peak_val = overall_stats['actual_max_peak']

    plt.figure(figsize=(14, 8))
    
    # --- Plotting Data ---
    if not bin_centers.empty and not frequencies.empty:
        if plot_type == 'histogram':
            # Calculate bin edges for plt.hist
            # Assumes bin_centers are sorted and bins are contiguous with width bin_width
            np_bin_centers = bin_centers.to_numpy() if isinstance(bin_centers, pd.Series) else np.array(bin_centers)
            np_frequencies = frequencies.to_numpy() if isinstance(frequencies, pd.Series) else np.array(frequencies)

            bin_edges = np.append(np_bin_centers - bin_width / 2.0,
                                  np_bin_centers[-1] + bin_width / 2.0)

            plt.hist(np_bin_centers, bins=bin_edges, weights=np_frequencies,
                     alpha=0.7, label='Histogram', edgecolor='black', linewidth=0.5, rwidth=0.9)
        elif plot_type == 'scatter':
             # Plot bin centers vs frequencies as scatter points
             plt.scatter(bin_centers, frequencies, 
                         alpha=0.6, s=10, label='Data Points', color='blue') # s is marker size
        else:
             print(f"Error: Unknown plot_type '{plot_type}'. Skipping data plot.")
             # Still proceed to plot fit and stats lines if available
    else:
        print("Warning: No data to plot.")

    # Note: Statistical lines (Mean, Std Dev, Min Data, Max Peak) are based on the *original* data
    # as calculated in main(), which might include signal peaks.
    # This is intentional to show the overall data characteristics.
    plt.axvline(mean_val, color='red', linestyle='dashed', linewidth=1.5, label=f'Mean: {mean_val:.4g}')
    
    num_sigma_lines = 5 if snr_value > 100 else 1
    if std_val > 1e-9:
        sigma_colors = ['darkorange', 'gold', 'lightgreen', 'deepskyblue', 'violet']
        for i in range(1, num_sigma_lines + 1):
            color = sigma_colors[i-1] if i-1 < len(sigma_colors) else 'gray'
            
            current_label = None
            if i == 1: # Label for +/-1 sigma
                current_label = f'Mean ± {i}σ'
            elif i == num_sigma_lines and num_sigma_lines > 1: # Label for outermost sigma if not 1st
                current_label = f'Mean ± {i}σ'
            
            plt.axvline(mean_val + i * std_val, color=color, linestyle='dotted', linewidth=1.2, label=current_label)
            if i * std_val > 1e-9: # Avoid plotting on top of itself if std_val is effectively zero
                 plt.axvline(mean_val - i * std_val, color=color, linestyle='dotted', linewidth=1.2) # No label for the second line of the pair

    plt.axvline(min_data_val, color='blue', linestyle='dashdot', linewidth=1.5, label=f'Min Data Amp: {min_data_val:.4g}')
    plt.axvline(actual_max_peak_val, color='purple', linestyle='dashdot', linewidth=1.5, label=f'Max Peak Amp: {actual_max_peak_val:.4g}')

    # --- Plotting Fitted Rayleigh Distribution ---
    if fitted_sigma is not None and fitted_sigma > 1e-9 and not bin_centers.empty:
        # Generate points for the fitted PDF. Use a denser linspace for a smoother curve.
        # Ensure the range covers the data range, starting from 0.
        x_plot_fit = np.linspace(0, bin_centers.max() * 1.1, 500) # Extend slightly beyond max data for visualization
        
        y_fit_pdf = rayleigh_pdf_func(x_plot_fit, fitted_sigma) # PDF values
        
        # Scale the PDF to match the data's frequency count scale
        # Use the total frequency count and bin width of the *plotted* data (after rebinning)
        total_plotted_frequency = frequencies.sum()
        if total_plotted_frequency > 0 and bin_width > 1e-9:
             y_fit_counts = y_fit_pdf * total_plotted_frequency * bin_width
             plt.plot(x_plot_fit, y_fit_counts, color='cyan', linestyle='-', linewidth=2, label=f'Rayleigh Fit (σ={fitted_sigma:.4g})')
        else:
             print("Warning: Cannot scale fitted PDF to frequency counts (zero total frequency or bin width).")


    plt.title(plot_title, fontsize=16)
    plt.xlabel('Amplitude', fontsize=14)
    plt.ylabel('Frequency Count', fontsize=14) # Y-label is Frequency Count for both histogram and scatter of binned data
    
    # Create a comprehensive legend, place it carefully
    handles, labels = plt.gca().get_legend_handles_labels()
    # Remove duplicate labels if any (e.g. from multiple sigma lines)
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=20) # Increased fontsize
    
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.tick_params(axis='both', which='major', labelsize=12)

    # --- Set Plot Limits and Save ---
    plot_filename = output_filename_base # Default filename base

    if x_limit:
        plt.xlim(x_limit)
        plot_filename = f"{output_filename_base}_zoomed.png"
        # For zoomed plot, adjust y-limit to focus on visible data
        if x_limit is not None and not bin_centers[(bin_centers >= x_limit[0]) & (bin_centers <= x_limit[1])].empty:
            relevant_frequencies = frequencies[(bin_centers >= x_limit[0]) & (bin_centers <= x_limit[1])]
            if not relevant_frequencies.empty:
                 plt.ylim(0, relevant_frequencies.max() * 1.15 if relevant_frequencies.max() > 0 else 10)
            else: # No data in the zoomed x-range
                 plt.ylim(0, 10) # Set a small default y-limit
        else: # No bin centers in the zoomed x-range
             plt.ylim(0, 10) # Set a small default y-limit

    else: # Full plot
        plot_filename = f"{output_filename_base}.png" # Use the base filename directly for the full plot
        if not frequencies.empty:
            # Set y-limit based on the full (possibly rebinned) data
            plt.ylim(0, frequencies.max() * 1.15 if frequencies.max() > 0 else 10)
        else:
             plt.ylim(0, 10) # Set a small default y-limit


    plt.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make space for suptitle
    plt.suptitle(f"File: {os.path.basename(input_filename)} (SNR: {snr_value:.2f})", fontsize=10, y=0.99)
    
    try:
        plt.savefig(plot_filename)
        print(f"Saved plot to {plot_filename}")
    except Exception as e:
        print(f"Error saving plot to {plot_filename}: {e}")
    plt.show()
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Analyze Rayleigh CSV file from frinZsearch and plot amplitude distribution.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("csv_filepath", help="Path to the .rayleigh.csv file")
    parser.add_argument(
        "--rebin",
        type=int,
        default=1,
        help="Factor by which to rebin the histogram (number of original bins to combine into one new bin). Default: 1 (no rebinning)."
        )
    parser.add_argument(
        "--fit-max-amp",
        type=float,
        default=None,
        help="Maximum amplitude to include in the Rayleigh fitting. Useful to exclude signal peaks. Default: None (fit all positive amplitudes)."
    )
    parser.add_argument(
        "--plot-type",
        type=str,
        choices=['histogram', 'scatter'],
        default='histogram',
        help="Type of plot for the data: 'histogram' or 'scatter'. Default: 'histogram'."
    )
    args = parser.parse_args()

    try:
        df_original, snr, min_amp_data, actual_max_peak_amp, comments = parse_rayleigh_csv(args.csv_filepath)
    except Exception as e:
        print(f"Error parsing CSV file '{args.csv_filepath}': {e}")
        return

    if df_original.empty:
        print("No data to process from CSV.")
        return

    # These stats are calculated from the original, un-rebinned data
    original_bin_centers = df_original['AmplitudeBinCenter']
    original_frequencies = df_original['FrequencyCount']
    
    # Determine original_bin_width from comments or data
    original_bin_width = 0.0
    parsed_bw_comment = comments.get('Bin Width')
    if parsed_bw_comment:
        try:
            original_bin_width = float(parsed_bw_comment)
        except ValueError:
            print(f"Warning: Could not parse 'Bin Width' ('{parsed_bw_comment}') from comments.")
    
    if original_bin_width <= 1e-9: 
        print(f"Warning: 'Bin Width' from comments is invalid ({original_bin_width}). Attempting calculation from original_bin_centers.")
        if len(original_bin_centers) > 1:
            calculated_width = original_bin_centers.iloc[1] - original_bin_centers.iloc[0]
            if calculated_width > 1e-9:
                original_bin_width = calculated_width
                print(f"Info: Calculated original_bin_width as {original_bin_width:.4g} from data.")
            else:
                original_bin_width = 1.0 
                print(f"Critical Warning: Could not determine valid original_bin_width from data. Using fallback: {original_bin_width}.")
        elif len(original_bin_centers) == 1:
             original_bin_width = 1.0
             print(f"Warning: Single bin in data and 'Bin Width' comment invalid/unusable. Using fallback width for this bin: {original_bin_width}.")
        else: 
            original_bin_width = 1.0
            print(f"Warning: No bin data to process. Using fallback width: {original_bin_width}.")


    overall_stats_dict = {
        'mean': calculate_weighted_stats(original_bin_centers, original_frequencies)[0],
        'std': calculate_weighted_stats(original_bin_centers, original_frequencies)[1],
        'min_data': min_amp_data,          # Min amplitude of the binned data
        'actual_max_peak': actual_max_peak_amp # Actual max amplitude from FFT peak params
    }
    
    # --- Rebinning ---
    current_bin_centers, current_frequencies, current_bin_width = rebin_data(
        original_bin_centers, original_frequencies, original_bin_width, args.rebin
    )

    # --- Rayleigh Fitting ---
    fitted_sigma = fit_rayleigh_to_data(
        current_bin_centers, current_frequencies, current_bin_width, 
        overall_stats_dict, args.fit_max_amp
    )

    # Generate a base name for output files
    output_plot_filename_base = os.path.splitext(args.csv_filepath)[0] # Use filename without extension

    # --- Plotting ---
    # Plot full range
    plot_amplitude_data(
        current_bin_centers, current_frequencies, current_bin_width, 
        overall_stats_dict, fitted_sigma, snr, args.plot_type,
        plot_title=f"Amplitude Distribution (Full Range)",
        input_filename=args.csv_filepath,
        output_filename_base=output_plot_filename_base
    )

    # Plot zoomed range if SNR is high
    if snr > 100 and overall_stats_dict['std'] > 1e-9:
        lower_limit_5sigma = overall_stats_dict['mean'] - 5 * overall_stats_dict['std']
        upper_limit_5sigma = overall_stats_dict['mean'] + 5 * overall_stats_dict['std']
        
        # Ensure the 5-sigma range is sensible (e.g., lower_limit < upper_limit and lower_limit >= 0)
        lower_limit_5sigma = max(0, lower_limit_5sigma) # Amplitude cannot be negative
        
        if lower_limit_5sigma >= upper_limit_5sigma:
            print(f"Warning: Calculated 5-sigma range ({lower_limit_5sigma:.4g} to {upper_limit_5sigma:.4g}) is invalid or starts at/below zero. Skipping zoomed plot.")
        else:
            plot_amplitude_data(
                current_bin_centers, current_frequencies, current_bin_width, 
                overall_stats_dict, fitted_sigma, snr, args.plot_type,
                plot_title=f"Amplitude Distribution (Mean ± 5σ Range)",
                input_filename=args.csv_filepath,
                x_limit=(lower_limit_5sigma, upper_limit_5sigma),
                output_filename_base=output_plot_filename_base
            ) 

if __name__ == "__main__":
    main()
