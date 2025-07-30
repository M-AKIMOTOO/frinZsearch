use clap::{Parser, CommandFactory};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use num_complex::Complex;
use plotters::prelude::*;
use plotters::style::colors::colormaps::ViridisRGB;

// frinZsearchの他のモジュールを読み込む
mod frinz_read;
mod frinz_error;

use frinz_error::FrinZError;
use frinz_read::{read_binary_file_data, HeaderRegion};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Specify input binary file (required)
    #[arg(short, long, aliases = ["i", "in", "inp", "inpu"])]
    input: String,

    /// Output the read spectrum to a text file
    #[arg(short, long, aliases = ["o", "ou", "out", "outp", "outpu"])]
    output: bool,

    /// Generate and save spectrum heatmaps (amplitude and phase)
    #[arg(short, long, aliases = ["p", "pl", "plo"])]
    plot: bool,
}

fn main() {
    let args = match Args::try_parse() {
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

    let output_dir = if args.output || args.plot {
        let input_path = Path::new(&args.input);
        let parent_dir = input_path.parent().unwrap_or_else(|| Path::new("."));
        let frinz_dir = parent_dir.join("frinZ").join("frinZrawvis");
        if !frinz_dir.exists() {
            fs::create_dir_all(&frinz_dir).expect("Failed to create output directory");
        }
        Some(frinz_dir)
    } else {
        None
    };

    match read_binary_file_data(&args.input) {
        Ok(file_data) => {
            println!("File reading successful!");

            if let Some(dir) = output_dir {
                let input_basename = Path::new(&args.input).file_stem().unwrap().to_str().unwrap();

                if args.output {
                    let output_filename = format!("{}.spectrum.txt", input_basename);
                    let output_filepath = dir.join(output_filename);
                    // For text output, we still write the first spectrum as a sample
                    if let Some(first_spectrum) = file_data.spectrum_data.get(0) {
                        match write_spectrum_to_text(&output_filepath, &file_data.header, first_spectrum) {
                            Ok(_) => println!("Spectrum data written to {}", output_filepath.display()),
                            Err(e) => eprintln!("Failed to write spectrum to text file: {}", e),
                        }
                    } else {
                        eprintln!("No spectrum data available to write to text file.");
                    }
                }

                if args.plot {
                    let amp_filename = format!("{}_heatmap_amp.png", input_basename);
                    let phase_filename = format!("{}_heatmap_phs.png", input_basename);
                    let amp_filepath = dir.join(amp_filename);
                    let phase_filepath = dir.join(phase_filename);

                    match plot_spectrum_heatmaps(&amp_filepath, &phase_filepath, &file_data.spectrum_data) {
                        Ok(_) => println!(
                            "Spectrum heatmaps saved to {} and {}",
                            amp_filepath.display(),
                            phase_filepath.display()
                        ),
                        Err(e) => eprintln!("Failed to plot spectrum heatmaps: {}", e),
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("File reading error: {}", e);
        }
    }
}

pub fn write_spectrum_to_text<P: AsRef<Path>>(
    output_path: P,
    header: &HeaderRegion,
    spectrum_data: &Vec<Complex<f32>>,
) -> Result<(), FrinZError> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);

    // Copy fields from packed struct to local variables to avoid alignment issues.
    let magic_word = header.magic_word;
    let header_version = header.header_version;
    let software_version = header.software_version;
    let sampling_speed = header.sampling_speed;
    let observed_sky_freq = header.observed_sky_freq;
    let fft_point = header.fft_point;
    let num_sector = header.num_sector;
    let source_ra = header.source_ra;
    let source_dec = header.source_dec;

    writeln!(writer, "# Header Information")?;
    writeln!(writer, "# Magic Word: {:#X}", magic_word)?;
    writeln!(writer, "# Header Version: {}", header_version)?;
    writeln!(writer, "# Software Version: {}", software_version)?;
    writeln!(writer, "# Sampling Speed: {} Hz", sampling_speed)?;
    writeln!(writer, "# Observed Sky Freq: {} Hz", observed_sky_freq)?;
    writeln!(writer, "# FFT Point: {}", fft_point)?;
    writeln!(writer, "# Number of Sectors in file: {}", num_sector)?;
    let station1_name = String::from_utf8_lossy(&header.station1_name).trim_end_matches('\0').to_string();
    writeln!(writer, "# Station1 Name: {}", station1_name)?;
    let station2_name = String::from_utf8_lossy(&header.station2_name).trim_end_matches('\0').to_string();
    writeln!(writer, "# Station2 Name: {}", station2_name)?;
    let source_name = String::from_utf8_lossy(&header.source_name).trim_end_matches('\0').to_string();
    writeln!(writer, "# Source Name: {}", source_name)?;
    writeln!(writer, "# Source RA: {} rad", source_ra)?;
    writeln!(writer, "# Source Dec: {} rad", source_dec)?;
    writeln!(writer, "# --- Spectrum Data ---")?;
    writeln!(writer, "# FFT_Point, Real, Imag, Amplitude, Phase(deg)")?;

    for (i, data) in spectrum_data.iter().enumerate() {
        writeln!(
            writer,
            "{}, {}, {}, {}, {}",
            i, data.re, data.im, data.norm(), data.arg().to_degrees()
        )?;
    }

    Ok(())
}

pub fn plot_spectrum_heatmaps<P: AsRef<Path>>(
    output_path_amplitude: P,
    output_path_phase: P,
    spectrum_data: &Vec<Vec<Complex<f32>>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let output_path_amplitude = output_path_amplitude.as_ref().to_path_buf();
    let output_path_phase = output_path_phase.as_ref().to_path_buf();
    if spectrum_data.is_empty() || spectrum_data[0].is_empty() {
        return Err("Spectrum data for heatmap is empty".into());
    }

    let num_sectors = spectrum_data.len();
    let fft_points = spectrum_data[0].len();

    // --- Amplitude Heatmap ---
    let color_bar_width = 50;
    let color_bar_padding = 20; // Padding between chart and color bar
    let main_chart_width = 1000;
    let total_width_amp = main_chart_width + color_bar_width + color_bar_padding +70;
    let total_height_amp = 768;

    let root_amp = BitMapBackend::new(&output_path_amplitude, (total_width_amp, total_height_amp));
    let root_amp_drawing_area = root_amp.into_drawing_area();
    root_amp_drawing_area.fill(&WHITE)?;

    let (chart_area_amp, color_bar_area_for_amp) = root_amp_drawing_area.split_horizontally(main_chart_width);

    let amplitudes: Vec<f32> = spectrum_data.iter().flatten().map(|c| c.norm()).collect();
    let max_amp = amplitudes.iter().cloned().fold(0.0, f32::max);
    let min_amp = amplitudes.iter().cloned().fold(f32::MAX, f32::min);

    let mut chart_amp = ChartBuilder::on(&chart_area_amp)
        .caption("Complex Spectrum - Amplitude Heatmap", ("sans-serif", 20).into_font())
        .margin(10)
        .x_label_area_size(70)
        .y_label_area_size(70)
        .build_cartesian_2d(0..fft_points, 0..num_sectors)?;

    chart_amp.configure_mesh()
        .x_desc("Channels")
        .y_desc("PP")
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()?;

    chart_amp.draw_series(
        (0..fft_points).flat_map(|x| (0..num_sectors).map(move |y| (x, y)))
        .map(|(x, y)| {
            let amp = spectrum_data[y][x].norm();
            let color_value = if max_amp > min_amp { (amp - min_amp) / (max_amp - min_amp) } else { 0.0 };
            let color = ViridisRGB.get_color(color_value as f64);
            Rectangle::new([(x, y), (x + 1, y + 1)], color.filled())
        })
    )?;

    // Draw color bar for amplitude
    let color_bar_height_for_drawing = total_height_amp as i32 - (10 + 10 + 70 + 25); // total_height - (top_margin + bottom_margin + x_label_area_size)
    let color_bar_y_offset_for_drawing = 35; // top_margin

    for i in 0..color_bar_height_for_drawing {
        let color_value = i as f64 / color_bar_height_for_drawing as f64;
        let color = ViridisRGB.get_color(color_value);
        color_bar_area_for_amp.draw(&Rectangle::new(
            [(0, color_bar_y_offset_for_drawing + i), (color_bar_width as i32, color_bar_y_offset_for_drawing + i + 1)],
            color.filled(),
        ))?;
    }

    // Add labels to the color bar
    let num_labels = 5;
    for i in 0..num_labels {
        let value = min_amp + (max_amp - min_amp) * (i as f32 / (num_labels - 1) as f32);
        let y_pos_for_drawing = color_bar_y_offset_for_drawing + color_bar_height_for_drawing - (i as i32 * color_bar_height_for_drawing / (num_labels - 1) as i32);
        color_bar_area_for_amp.draw_text(
            &format!("{:.1e}", value),
            &TextStyle::from(("sans-serif", 25).into_font()).color(&BLACK),
            ((color_bar_width + 5) as i32, y_pos_for_drawing - 7),
        )?;
    }
    
    root_amp_drawing_area.present()?;
    


    // --- Phase Heatmap ---
    let color_bar_width = 50;
    let color_bar_padding = 30; // Padding between chart and color bar
    let main_chart_width = 1000;
    let total_width_phase = main_chart_width + color_bar_width + color_bar_padding +50;
    let total_height_phase = 768;

    let root_phase = BitMapBackend::new(&output_path_phase, (total_width_phase, total_height_phase));
    let root_phase_drawing_area = root_phase.into_drawing_area();
    root_phase_drawing_area.fill(&WHITE)?;

    let (chart_area_phase, color_bar_area_for_phase) = root_phase_drawing_area.split_horizontally(main_chart_width);

    let mut chart_phase = ChartBuilder::on(&chart_area_phase)
        .caption("Complex Spectrum - Phase Heatmap", ("sans-serif", 20).into_font())
        .margin(10)
        .x_label_area_size(70)
        .y_label_area_size(70)
        .build_cartesian_2d(0..fft_points, 0..num_sectors)?;

    chart_phase.configure_mesh()
        .x_desc("Channels")
        .y_desc("PP")
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()?;

    chart_phase.draw_series(
        (0..fft_points).flat_map(|x| (0..num_sectors).map(move |y| (x, y)))
        .map(|(x, y)| {
            let phase_deg = spectrum_data[y][x].arg().to_degrees();
            // Normalize phase from -180..180 to 0..1 for the color map
            let color_value = (phase_deg + 180.0) / 360.0;
            let color = ViridisRGB.get_color(color_value as f64);
            Rectangle::new([(x, y), (x + 1, y + 1)], color.filled())
        })
    )?;

    // Draw color bar for phase
    let color_bar_height_i32 = total_height_phase as i32 - (10 + 10 + 70 + 25); // total_height - (top_margin + bottom_margin + x_label_area_size)
    let color_bar_y_offset_i32 = 35; // top_margin

    for i in 0..color_bar_height_i32 {
        let color_value = i as f64 / color_bar_height_i32 as f64;
        let color = ViridisRGB.get_color(color_value);
        color_bar_area_for_phase.draw(&Rectangle::new(
            [(0, color_bar_y_offset_i32 + i), (color_bar_width as i32, color_bar_y_offset_i32 + i + 1)],
            color.filled(),
        ))?;
    }

    // Add labels to the color bar
    let num_labels = 9;
    for i in 0..num_labels {
        let value = -180.0 + (360.0) * (i as f32 / (num_labels - 1) as f32);
        let y_pos_i32 = color_bar_y_offset_i32 + color_bar_height_i32 - (i as i32 * color_bar_height_i32 / (num_labels - 1) as i32);
        color_bar_area_for_phase.draw_text(
            &format!("{:.0}", value),
            &TextStyle::from(("sans-serif", 30).into_font()).color(&BLACK),
            ((color_bar_width + 5) as i32, y_pos_i32 - 7),
        )?;
    }

    root_phase_drawing_area.present()?;

    Ok(())
}
