use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use byteorder::{LittleEndian, ReadBytesExt};
use num_complex::Complex;
use chrono::{TimeZone, Utc};

use crate::frinZerror::FrinZError;

// C++のHeaderRegion構造体に対応
// #[repr(C, packed)] を使い、C言語と同じメモリレイアウトにすることを保証する
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
    // 予約領域を埋める
    _reserved: [u8; 100],
}

// C++のFileData構造体に対応
#[derive(Debug)]
pub struct FileData {
    pub header: HeaderRegion,
    pub spectrum_data: Vec<Vec<Complex<f32>>>,
    pub sector_start_times_utc: Vec<String>, // C++版の対応
    pub first_effective_integration_length: f32, // C++版の対応
    pub first_sector_epoch_time: u32, // C++版の対応
}

// ヘッダーを読み込む関数
fn read_header(file: &mut File) -> Result<HeaderRegion, FrinZError> {
    // HeaderRegionと完全に同じサイズのバッファを確保
    let mut buffer = [0u8; std::mem::size_of::<HeaderRegion>()];
    file.read_exact(&mut buffer)?;

    // バイト列をHeaderRegion構造体に変換
    // C/C++の構造体をそのままメモリにマップするのと同様の操作
    // 安全でない操作なので unsafe ブロックで囲む
    let header = unsafe {
        std::ptr::read(buffer.as_ptr() as *const _)
    };
    Ok(header)
}

pub fn read_binary_file_data<P: AsRef<Path>>(filepath: P) -> Result<FileData, FrinZError> {
    let mut file = File::open(filepath)?;

    // ヘッダー読み込み
    let header = read_header(&mut file)?;

    let mut spectrum_data = Vec::with_capacity(header.num_sector as usize);
    let mut sector_start_times_utc = Vec::with_capacity(header.num_sector as usize);
    let mut first_effective_integration_length = 0.0f32;
    let mut first_sector_epoch_time = 0u32;
    let mut first_sector_processed = false;

    // セクターデータの読み込み
    for s in 0..header.num_sector {
        // 各セクターの開始位置にシーク
        let sector_offset = 256 + s * (128 + 4 * header.fft_point);
        file.seek(SeekFrom::Start(sector_offset as u64))?;

        let corr_start_sec_raw = file.read_u32::<LittleEndian>()?;
        if s == 0 {
            first_sector_epoch_time = corr_start_sec_raw;
        }

        // C++版のコードでは、corr_start_sec_rawの後に108バイトスキップしている
        file.seek(SeekFrom::Current(108))?;

        let effective_integration_length_current = file.read_f32::<LittleEndian>()?;
        if !first_sector_processed {
            first_effective_integration_length = effective_integration_length_current;
            first_sector_processed = true;
        }

        // C++版のコードでは、effective_integration_lengthの後に12バイトスキップしている
        file.seek(SeekFrom::Current(12))?;

        let mut current_sector_spectrum = Vec::with_capacity((header.fft_point / 2) as usize);
        for _ in 0..(header.fft_point / 2) {
            let real = file.read_f32::<LittleEndian>()?;
            let imag = file.read_f32::<LittleEndian>()?;
            current_sector_spectrum.push(Complex::new(real, imag));
        }
        spectrum_data.push(current_sector_spectrum);

        // UTC時刻文字列の生成 (C++版の形式に合わせる)
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
        first_sector_epoch_time,
    })
}
