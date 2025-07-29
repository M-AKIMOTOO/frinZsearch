#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <filesystem>
#include "frinZread.hpp"
#include "frinZoutput.hpp"
#include "frinZlogger.hpp"
#include "frinZargs.hpp"

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    ProgramOptions params;
    Logger logger; // Logger instance

    std::string input_filepath;
    bool output_spectrum = false;
    bool plot_spectrum = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--input" || arg == "-i") {
            if (i + 1 < argc) {
                input_filepath = argv[++i];
            }
        } else if (arg == "--output") {
            output_spectrum = true;
        } else if (arg == "--plot") {
            plot_spectrum = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0] << " --input <filepath> [--output-spectrum] [--plot-spectrum]" << std::endl;
            return 0;
        }
    }

    if (input_filepath.empty()) {
        std::cerr << "Error: Input file not specified. Use --input <filepath>." << std::endl;
        return 1;
    }

    logger.setup(true, false, ""); // Console output enabled, file output disabled

    FileData loaded_data = read_binary_file_data(input_filepath, logger, false);
    if (!loaded_data.success) {
        std::cerr << "Error reading binary file: " << loaded_data.error_message << std::endl;
        return 1;
    }

    fs::path input_path_fs(input_filepath);
    fs::path parent_dir = input_path_fs.parent_path();
    fs::path frinz_dir = parent_dir / "frinZ" / "frinZread";

    if (!fs::exists(frinz_dir)) {
        fs::create_directories(frinz_dir);
    }

    std::string input_basename = input_path_fs.stem().string();

    if (output_spectrum) {
        std::string output_filename = input_basename + ".spectrum.txt";
        std::string output_filepath = (frinz_dir / output_filename).string();
        if (!loaded_data.spectrum_data.empty()) {
            write_spectrum_to_text(output_filepath, loaded_data.header, loaded_data.spectrum_data);
            std::cout << "Spectrum data written to: " << output_filepath << std::endl;
        } else {
            std::cerr << "No spectrum data available to write to text file." << std::endl;
        }
    }

    if (plot_spectrum) {
        std::string data_filename = input_basename + ".spectrum.txt";
        std::string data_filepath = (frinz_dir / data_filename).string();

        if (!fs::exists(data_filepath)) {
            if (!loaded_data.spectrum_data.empty()) {
                write_spectrum_to_text(data_filepath, loaded_data.header, loaded_data.spectrum_data);
                std::cout << "Generated data file for plotting: " << data_filepath << std::endl;
            } else {
                std::cerr << "No spectrum data available to generate for plotting." << std::endl;
                return 1;
            }
        }

        std::string amp_plot_filename = input_basename + "_heatmap_amp.png";
        std::string phase_plot_filename = input_basename + "_heatmap_phs.png";
        std::string amp_plot_filepath = (frinz_dir / amp_plot_filename).string();
        std::string phase_plot_filepath = (frinz_dir / phase_plot_filename).string();

        plot_spectrum_heatmaps_with_gnuplot(
            amp_plot_filepath,
            phase_plot_filepath,
            loaded_data.spectrum_data,
            loaded_data.header.fft_point,
            loaded_data.header.num_sector
        );
        std::cout << "Spectrum plots generated: " << amp_plot_filepath << " and " << phase_plot_filepath << std::endl;
    }

    return 0;
}
