#ifndef FRINZARGS_HPP
#define FRINZARGS_HPP

#include <string>
#include <vector>
#include <limits>
#include <filesystem>

// ProgramOptions 構造体の定義
struct ProgramOptions {
    std::string input_filename;
    double specified_length_sec = -1.0;
    double skip_seconds = 0.0;
    int loop_count = 1;
    int iterations = 3;
    bool iter_explicitly_set = false;
    std::string output_dir_final;
    std::string output_calibrated_filename_base = "calibrated_fft_amplitude";
    bool frequency_mode = false;
    std::vector<std::pair<double, double>> rfi_ranges_mhz;
    double initial_delay_samples = 0.0;
    bool initial_delay_specified = false;
    double initial_rate_hz = 0.0;
    bool initial_rate_specified = false;
    double delay_search_min = -std::numeric_limits<double>::infinity();
    double delay_search_max = std::numeric_limits<double>::infinity();
    bool delay_search_range_specified = false;
    double rate_search_min = -std::numeric_limits<double>::infinity();
    double rate_search_max = std::numeric_limits<double>::infinity();
    bool rate_search_range_specified = false;
    bool output_header_info = false;
    bool quick_mode = false;
    bool noconsole = false;
    bool enable_text_log_output = false;
    int delay_fit_points = 3; // Number of points for delay fit (center +/- (N-1)/2 )
    int rate_fit_points = 3;  // Number of points for rate fit (center +/- (N-1)/2 )
    bool output_rayleigh_csv = false; // For --Rayleigh option
};

// 引数を解析し、ProgramOptions 構造体に設定する関数
// 成功した場合は true、ヘルプ/バージョン表示やエラーの場合は false を返す
bool parse_arguments(int argc, char* argv[], ProgramOptions& params);

// パース後にオプションを最終処理する関数
void post_process_options(ProgramOptions& params);

#endif // FRINZARGS_HPP
