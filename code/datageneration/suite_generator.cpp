#include "generators.h"
#include <exception>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>

namespace fs = std::filesystem;

namespace DataGeneration {

bool make_test_suite(std::string &path, size_t start, size_t step, size_t end, size_t cases, std::mt19937 &rd,
										 float min_length, float max_length) {
	if (fs::exists(path)) {
		return false;
	}

	try {
		fs::create_directory(path);
		fs::path well_behaved = fs::path(path) / "well-behaved";
		fs::path non_well_behaved = fs::path(path) / "non-well-behaved";

		fs::create_directory(well_behaved);
		fs::create_directory(non_well_behaved);

		for (unsigned int i = start; i <= end; i += step) {
			fs::create_directory(well_behaved  / std::to_string(i));
			fs::create_directory(non_well_behaved  / std::to_string(i));

			for (unsigned int j = 0; j < cases; j++) {
				auto const polyline_w = DataGeneration::make_polyline(i, 2, min_length, max_length, 60 / 180.0 * std::numbers::pi, rd);
				DataGeneration::write_to_file(*polyline_w, well_behaved / std::to_string(i) / std::to_string(j));
				auto const polyline_n = DataGeneration::make_polyline(i, 2, min_length, max_length, 180 / 180.0 * std::numbers::pi, rd);
				DataGeneration::write_to_file(*polyline_n, non_well_behaved / std::to_string(i) / std::to_string(j));
			}
		}

		return true;
	} catch (std::exception const &e) {
		std::cerr << "Error creating directories: " << e.what() << std::endl;
		return false;
	}
}

}
