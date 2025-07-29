use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use byteorder::{LittleEndian, ByteOrder};
use num_complex::Complex;
use chrono::{TimeZone, Utc};

use crate::frinz_error::FrinZError;

/// Reads and returns the content of a specified file from the local filesystem.
#[repr(C, packed)]
#[derive(Debug, Clone, Copy)]
pub struct HeaderRegion {
    pub magic_word: u32,
    pub header_version: u32,
    pub software_version: u32,
    pub sampling_speed: u32,
    pub observed_sky_freq: f64,
    pub fft_point: u32,
    pub num_sector: u32,
    pub station1_name: [u8; 16],
    pub station1_pos_x: f64,
    pub station1_pos_y: f64,
    pub station1_pos_z: f64,
    pub station1_key: [u8; 8],
    pub station2_name: [u8; 16],
    pub station2_pos_x: f64,
    pub station2_pos_y: f64,
    pub station2_pos_z: f64,
    pub station2_key: [u8; 8],
    pub source_name: [u8; 16],
    pub source_ra: f64,
    pub source_dec: f64,
    // Fill reserved area
    _reserved: [u8; 100],
}

/// Reads and returns the content of a specified file from the local filesystem.
#[derive(Debug)]
pub struct FileData {
    pub header: HeaderRegion,
    pub spectrum_data: Vec<Vec<Complex<f32>>>,
    pub sector_start_times_utc: Vec<String>, // Corresponds to C++ version
    pub first_effective_integration_length: f32, // Corresponds to C++ version
    // pub first_sector_epoch_time: u32, // C++版の対応
}

/// Reads and returns the content of a specified file from the local filesystem.
fn read_header(file: &mut File) -> Result<HeaderRegion, FrinZError> {
    // HeaderRegionと完全に同じサイズのバッファを確保
    let mut buffer = [0u8; std::mem::size_of::<HeaderRegion>()];
    file.read_exact(&mut buffer)?;

    // Convert byte array to HeaderRegion struct
    // This operation is similar to memory mapping a C/C++ struct directly.
    // It's an unsafe operation, so it's enclosed in an unsafe block.
    let header = unsafe {
        std::ptr::read(buffer.as_ptr() as *const _)
    };
    Ok(header)
}

pub fn read_binary_file_data<P: AsRef<Path>>(filepath: P) -> Result<FileData, FrinZError> {
    let mut file = File::open(filepath)?;

    // Function to read header
    let header = read_header(&mut file)?;

    let mut spectrum_data = Vec::with_capacity(header.num_sector as usize);
    let mut sector_start_times_utc = Vec::with_capacity(header.num_sector as usize);
    let mut first_effective_integration_length = 0.0f32;
    // let mut first_sector_epoch_time = 0u32;
    let mut first_sector_processed = false;

    // Calculate the size of one sector data block (excluding header and initial 256 bytes)
    let sector_data_size = 128 + 4 * header.fft_point; // 128 bytes for metadata + 4 bytes/complex * fft_point

    // Read sector data
    for s in 0..header.num_sector {
        // Seek to the start of each sector
        let sector_offset = 256 + s * sector_data_size;
        file.seek(SeekFrom::Start(sector_offset as u64))?;

        // Read the entire sector data block into a buffer
        let mut sector_buffer = vec![0u8; sector_data_size as usize];
        file.read_exact(&mut sector_buffer)?;

        // Directly read from slices of sector_buffer
        let mut offset_in_buffer = 0;

        let corr_start_sec_raw = LittleEndian::read_u32(&sector_buffer[offset_in_buffer..]);
        offset_in_buffer += 4; // sizeof(u32)

        offset_in_buffer += 108; // Skip 108 bytes

        let effective_integration_length_current = LittleEndian::read_f32(&sector_buffer[offset_in_buffer..]);
        offset_in_buffer += 4; // sizeof(f32)

        if !first_sector_processed {
            first_effective_integration_length = effective_integration_length_current;
            first_sector_processed = true;
        }

        offset_in_buffer += 12; // Skip 12 bytes

        let mut current_sector_spectrum = Vec::with_capacity((header.fft_point / 2) as usize);
        for _ in 0..(header.fft_point / 2) {
            let real = LittleEndian::read_f32(&sector_buffer[offset_in_buffer..]);
            offset_in_buffer += 4;
            let imag = LittleEndian::read_f32(&sector_buffer[offset_in_buffer..]);
            offset_in_buffer += 4;
            current_sector_spectrum.push(Complex::new(real, imag));
        }
        spectrum_data.push(current_sector_spectrum);

        // Generate UTC time string (matching C++ version format)
        let datetime = Utc.timestamp_opt(corr_start_sec_raw as i64, 0).unwrap();
        let time_buffer = datetime.format("%Y-%j %H:%M:%S UTC").to_string();
        let final_buffer = format!("{} (Epoch: {}) (Int.Length: {:.6} s)", 
                                   time_buffer, corr_start_sec_raw, effective_integration_length_current);
        sector_start_times_utc.push(final_buffer);
    }

    Ok(FileData {
        header,
        spectrum_data,
        sector_start_times_utc,
        first_effective_integration_length,
        // first_sector_epoch_time,
    })
}
