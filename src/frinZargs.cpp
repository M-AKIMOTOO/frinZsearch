#include "frinZargs.hpp"
#include <iostream> // For std::cerr
#include <filesystem> // For directory creation
#include <string>
#include <vector>
#include <algorithm> // For std::find_if if needed
#include <stdexcept>  // For std::invalid_argument, std::out_of_range

// ヘルパー関数: 文字列をdoubleに変換 (エラーチェック付き)
bool string_to_double(const std::string& s, double& value) {
    try {
        size_t processed_chars;
        value = std::stod(s, &processed_chars);
        // Ensure the entire string was processed, to catch cases like "1.2.3" or "12abc"
        return processed_chars == s.length();
    } catch (const std::invalid_argument& ia) {
        std::cerr << "エラー: 無効な数値形式です: " << s << std::endl;
        return false;
    } catch (const std::out_of_range& oor) {
        std::cerr << "エラー: 数値が範囲外です: " << s << std::endl;
        return false;
    }
}

// ヘルパー関数: 文字列をintに変換 (エラーチェック付き)
bool string_to_int(const std::string& s, int& value) {
    try {
        size_t processed_chars;
        value = std::stoi(s, &processed_chars);
        return processed_chars == s.length();
    } catch (const std::invalid_argument& ia) {
        std::cerr << "エラー: 無効な整数形式です: " << s << std::endl;
        return false;
    } catch (const std::out_of_range& oor) {
        std::cerr << "エラー: 整数が範囲外です: " << s << std::endl;
        return false;
    }
}

void print_help(const char* prog_name) {
    std::cout << "frinZsearch:" << std::endl;
    std::cout << ">> frinZ.py の拡張・補助ツール．frinZ.py は精密にフリンジの位置を推定できないので，それを補うのが frinZsearch である．" << std::endl;
    std::cout << ">> frinZsearch が推定した delay と rate を用いて frinZ.py を実行することで，遅延較正・位相較正が可能となる．" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "使用法: frinZsearch --input <file>" << " [オプション]" << std::endl;
    std::cout << "オプション:" << std::endl;
    std::cout << "    --input <file>        入力バイナリファイ，ルを指定 (必須)" << std::endl;
    std::cout << "    --length <seconds>    FFTセグメントの積分時間を指定 (秒，デフォルト: 全データ)" << std::endl;
    std::cout << "    --skip <seconds>      データの先頭からスキップする時間を指定 (秒，デフォルト: 0.0)" << std::endl;
    std::cout << "    --loop <count>        FFT処理をセグメントごとにループする回数 (デフォルト: 1)" << std::endl;
    std::cout << "    --iter <count>        位相較正のイテレーション回数 (デフォルト: 3)" << std::endl;
    std::cout << "    --output              コンソールメッセージとデータファイルを入力ファイル基準の 'frinZ/frinZsearch' サブディレクトリに出力" << std::endl;
    std::cout << "    --frequency           最初のFFTのみ実行し，較正/イテレーションなし" << std::endl;
    std::cout << "    --rfi <start,end> ... RFI周波数範囲をMHzで指定 (例: --rfi 10,50 100,120)" << std::endl;
    std::cout << "    --delay <samples>     初期遅延補正 (サンプル単位，デフォルト: 0.0)" << std::endl;
    std::cout << "    --rate <Hz>           初期レート補正 (Hz単位，デフォルト: 0.0)" << std::endl;
    std::cout << "    --drange <min max>    較正時の遅延探索範囲を制限 [min_samp max_samp]" << std::endl;
    std::cout << "    --rrange <min max>    較正時のレート探索範囲を制限 [min_Hz max_Hz]" << std::endl;
    std::cout << "    --dstep <points>      遅延二次フィッティングの点数 (3, 5 or 7，デフォルト: 3)" << std::endl;
    std::cout << "    --rstep <points>      レート二次フィッティングの点数 (3, 5 or 7，デフォルト: 3)" << std::endl;
    std::cout << "    --header              入力ファイルからヘッダー情報を出力" << std::endl;
    std::cout << "    --quick               位相較正をスキップし，初期FFT結果のみ出力" << std::endl;
    std::cout << "    --noconsole           コンソール出力を抑制 (--output 指定時はファイル出力)" << std::endl;
    std::cout << "    --help                このヘルプメッセージを表示して終了" << std::endl;
    std::cout << "    --version             バージョン情報を表示して終了" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "(c) M.AKIMOTO with Gemini in 2025/06/03" << std::endl;
    std::cout << "This program is licensed under the MIT License (see https://opensource.org/license/mit)" << std::endl;
}

void print_version() {
    std::cout << "frinZsearch v0.1 in 2025/05/31 (Made by M.AKIMOTO with Gemini)    \n"
                 "frinZsearch v1.0 in 2025/06/01   \n"
                 "frinZsearch v2.0 in 2025/06/03 " << std::endl;
}

bool parse_arguments(int argc, char* argv[], ProgramOptions& params) {
    if (argc == 1) { // 引数なしの場合
        print_help(argv[0]);
        return false; 
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help") {
            print_help(argv[0]);
            return false; // ヘルプ表示後，処理を続行しない
        } else if (arg == "--version") {
            print_version();
            return false; // バージョン表示後，処理を続行しない
        } else if (arg == "--input") {
            if (i + 1 < argc) {
                params.input_filename = argv[++i];
                if (!std::filesystem::exists(params.input_filename) || !std::filesystem::is_regular_file(params.input_filename)) {
                    std::cerr << "エラー: 入力ファイルが見つからないか，通常ファイルではありません: " << params.input_filename << std::endl;
                    return false;
                }
            } else {
                std::cerr << "エラー: --input にはファイル名が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--length") {
            if (i + 1 < argc) {
                if (!string_to_double(argv[++i], params.specified_length_sec)) return false;
            } else {
                std::cerr << "エラー: --length には秒数が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--skip") {
            if (i + 1 < argc) {
                if (!string_to_double(argv[++i], params.skip_seconds)) return false;
                if (params.skip_seconds < 0) {
                    std::cerr << "エラー: --skip の値は非負である必要があります．" << std::endl;
                    return false;
                }
            } else {
                std::cerr << "エラー: --skip には秒数が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--loop") {
            if (i + 1 < argc) {
                if (!string_to_int(argv[++i], params.loop_count)) return false;
                if (params.loop_count <= 0) {
                    std::cerr << "エラー: --loop の値は正の整数である必要があります．" << std::endl;
                    return false;
                }
            } else {
                std::cerr << "エラー: --loop には回数が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--iter") {
            if (i + 1 < argc) {
                if (!string_to_int(argv[++i], params.iterations)) return false;
                if (params.iterations < 0) {
                    std::cerr << "エラー: --iter の値は非負の整数である必要があります．" << std::endl;
                    return false;
                }
                params.iter_explicitly_set = true;
            } else {
                std::cerr << "エラー: --iter には回数が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--output") {
            params.enable_text_log_output = true;
        } else if (arg == "--frequency") {
            params.frequency_mode = true;
        } else if (arg == "--rfi") {
            while (i + 1 < argc && (argv[i+1][0] != '-' || (argv[i+1][0] == '-' && argv[i+1][1] != '-'))) {
                 // 次の引数が "--" で始まらない限り，RFI範囲として解釈
                std::string range_str = argv[++i];
                size_t comma_pos = range_str.find(',');
                if (comma_pos != std::string::npos) {
                    try {
                        double start_mhz = std::stod(range_str.substr(0, comma_pos));
                        double end_mhz = std::stod(range_str.substr(comma_pos + 1));
                        if (start_mhz < 0 || end_mhz < 0 || start_mhz > end_mhz) {
                            std::cerr << "エラー: --rfi の範囲が無効です: " << range_str << "．周波数は非負で start <= end である必要があります．" << std::endl;
                            return false;
                        }
                        params.rfi_ranges_mhz.emplace_back(start_mhz, end_mhz);
                    } catch (const std::exception& e) {
                        std::cerr << "エラー: --rfi の範囲形式が無効です: " << range_str << "．期待される形式: start_mhz,end_mhz．エラー: " << e.what() << std::endl;
                        return false;
                    }
                } else {
                    std::cerr << "エラー: --rfi の範囲形式が無効です: " << range_str << "．カンマ区切りのペアを期待します: start_mhz,end_mhz" << std::endl;
                    return false;
                }
            }
        } else if (arg == "--delay") {
            if (i + 1 < argc) {
                if (!string_to_double(argv[++i], params.initial_delay_samples)) return false;
                params.initial_delay_specified = true;
            } else {
                std::cerr << "エラー: --delay にはサンプル数が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--rate") {
            if (i + 1 < argc) {
                if (!string_to_double(argv[++i], params.initial_rate_hz)) return false;
                params.initial_rate_specified = true;
            } else {
                std::cerr << "エラー: --rate にはHz値が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--drange") {
            if (i + 2 < argc) {
                if (!string_to_double(argv[++i], params.delay_search_min)) return false;
                if (!string_to_double(argv[++i], params.delay_search_max)) return false;
                if (params.delay_search_min > params.delay_search_max) {
                    std::cerr << "エラー: --drange の min_samp (" << params.delay_search_min << ") > max_samp (" << params.delay_search_max << ") です．" << std::endl;
                    return false;
                }
                params.delay_search_range_specified = true;
            } else {
                std::cerr << "エラー: --drange には 2 つの数値 (min_samp max_samp) が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--rrange") {
            if (i + 2 < argc) {
                if (!string_to_double(argv[++i], params.rate_search_min)) return false;
                if (!string_to_double(argv[++i], params.rate_search_max)) return false;
                 if (params.rate_search_min > params.rate_search_max) {
                    std::cerr << "エラー: --rrange の min_Hz (" << params.rate_search_min << ") > max_Hz (" << params.rate_search_max << ") です．" << std::endl;
                    return false;
                }
                params.rate_search_range_specified = true;
            } else {
                std::cerr << "エラー: --rrange には2つの数値 (min_Hz max_Hz) が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--dstep") {
            if (i + 1 < argc) {
                if (!string_to_int(argv[++i], params.delay_fit_points)) return false;
                if (params.delay_fit_points < 3 || params.delay_fit_points % 2 == 0) {
                    std::cerr << "エラー: --dstep の値は（3, 5, 7）のいずれである必要があります．" << std::endl;
                    return false;
                }
            } else {
                std::cerr << "エラー: --dstep には点数（3, 5, 7）が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--rstep") {
            if (i + 1 < argc) {
                if (!string_to_int(argv[++i], params.rate_fit_points)) return false;
                 if (params.rate_fit_points < 3 || params.rate_fit_points % 2 == 0) {
                    std::cerr << "エラー: --rstep の値は（3, 5, 7）のいずれである必要があります．" << std::endl;
                    return false;
                }
            } else {
                std::cerr << "エラー: --rstep には点数（3, 5, 7）が必要です．" << std::endl;
                return false;
            }
        } else if (arg == "--header") {
            params.output_header_info = true;
        } else if (arg == "--quick") {
            params.quick_mode = true;
        } else if (arg == "--noconsole") {
            params.noconsole = true;
        }
         else {
            std::cerr << "エラー: 不明なオプションです: " << arg << std::endl;
            std::cout << "frinZsearch --help を実行してオプションを確認してください．" << std::endl;
            return false;
        }
    }

    // 必須オプションのチェック
    if (params.input_filename.empty()) {
        std::cerr << "エラー: --input オプションは必須です．" << std::endl;
        return false;
    }

    return true; // パース成功
}

void post_process_options(ProgramOptions& params) {
    if (params.frequency_mode || params.quick_mode) {
        params.iterations = 0;
    }
    // 元の finalize_options にあった --delay や --rate 指定時の --iter の挙動に関するコメントは，
    // 実際のコードでは frequency_mode または quick_mode の場合にのみ iterations を 0 に設定していました．
    // その挙動を維持します．もし --delay や --rate が指定され，かつ --iter が明示的に指定されていない場合に
    // iterations を 0 にしたい場合は，以下のようなロジックを追加します．
    // if ((params.initial_delay_specified || params.initial_rate_specified) && !params.iter_explicitly_set) {
    //     params.iterations = 0;
    // }

    if (params.enable_text_log_output) {
        if (params.input_filename.empty()) { 
            // このケースは parse_arguments でチェックされるはずですが，念のため
            std::cerr << "致命的なエラー: --output が有効ですが，入力ファイル名が空です．これは予期しない状態です．" << std::endl;
            params.enable_text_log_output = false; 
        } else {
            namespace fs = std::filesystem;
            fs::path input_file_path_fs(params.input_filename);
            fs::path absolute_input_path = fs::absolute(input_file_path_fs);
            fs::path base_output_dir = absolute_input_path.parent_path();
            
            params.output_dir_final = (base_output_dir / "frinZ" / "frinZsearch").string();
            try {
                fs::create_directories(params.output_dir_final);
            } catch (const fs::filesystem_error& e) {
                std::cerr << "エラー: 出力ディレクトリ '" << params.output_dir_final << "' の作成に失敗しました: " << e.what() << std::endl;
                params.enable_text_log_output = false; // Disable if dir creation fails
                params.output_dir_final.clear();
            }
        }
    }
}