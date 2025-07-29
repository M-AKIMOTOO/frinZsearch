use clap::Parser;
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
    #[arg(short, long)]
    input: String,

    /// Output the read spectrum to a text file
    #[arg(long)]
    output: bool,

    /// Generate and save spectrum heatmaps (amplitude and phase)
    #[arg(long)]
    plot: bool,
}

fn main() {
    let args = Args::parse();

    let output_dir = if args.output || args.plot {
        let input_path = Path::new(&args.input);
        let parent_dir = input_path.parent().unwrap_or_else(|| Path::new("."));
        let frinz_dir = parent_dir.join("frinZ").join("frinZread");
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
    if spectrum_data.is_empty() || spectrum_data[0].is_empty() {
        return Err("Spectrum data for heatmap is empty".into());
    }

    let num_sectors = spectrum_data.len();
    let fft_points = spectrum_data[0].len();

    // --- Amplitude Heatmap ---
    let root_amp = BitMapBackend::new(&output_path_amplitude, (1024, 768)).into_drawing_area();
    root_amp.fill(&WHITE)?;

    let amplitudes: Vec<f32> = spectrum_data.iter().flatten().map(|c| c.norm()).collect();
    let max_amp = amplitudes.iter().cloned().fold(0.0, f32::max);

    let mut chart_amp = ChartBuilder::on(&root_amp)
        .caption("Complex Spectrum - Amplitude Heatmap", ("sans-serif", 20).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0..fft_points, 0..num_sectors)?;

    chart_amp.configure_mesh()
        .x_desc("FFT Point")
        .y_desc("Sector")
        .draw()?;

    chart_amp.draw_series(
        (0..fft_points).flat_map(|x| (0..num_sectors).map(move |y| (x, y)))
        .map(|(x, y)| {
            let amp = spectrum_data[y][x].norm();
            let color_value = if max_amp > 0.0 { amp / max_amp } else { 0.0 };
            let color = ViridisRGB.get_color(color_value as f64);
            Rectangle::new([(x, y), (x + 1, y + 1)], color.filled())
        })
    )?;
    
    root_amp.present()?;

    // --- Phase Heatmap ---
    let root_phase = BitMapBackend::new(&output_path_phase, (1024, 768)).into_drawing_area();
    root_phase.fill(&WHITE)?;

    let mut chart_phase = ChartBuilder::on(&root_phase)
        .caption("Complex Spectrum - Phase Heatmap", ("sans-serif", 20).into_font())
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(0..fft_points, 0..num_sectors)?;

    chart_phase.configure_mesh()
        .x_desc("FFT Point")
        .y_desc("Sector")
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

    root_phase.present()?;

    Ok(())
}

