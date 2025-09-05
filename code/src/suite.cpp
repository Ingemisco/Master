#include "datastructures.h"
#include "log.h"
#include "simplification.h"
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using DataStructures::Polyline;
using std::vector;
using std::unique_ptr;

constexpr float _epsilon = 2.0;
constexpr float _epsilon2 = _epsilon * _epsilon;
AlgorithmConfiguration config = {};

namespace Log {

template <Log::Algorithm algorithm, Simplification::Simplification _simplify(Polyline const &, float, AlgorithmConfiguration &)>
void use_algorithm(
	vector<vector<unique_ptr<Polyline>>> const &w_polylines,
	vector<vector<unique_ptr<Polyline>>> const &n_polylines,
	Log::PerformanceLogger &log,
	float epsilon
) {
	for (auto const &vec : w_polylines) {
		auto const point_count = vec[0]->point_count;
		log.begin_data_set(algorithm, vec[0]->dimension, point_count, 
										Log::algorithm_name(algorithm) + ":w:" + std::to_string(point_count));

		for (auto const &poly : vec) {
			_simplify(*poly, epsilon, config);
		}
	}

	for (auto const &vec : n_polylines) {
		auto const point_count = vec[0]->point_count;
		log.begin_data_set(algorithm, vec[0]->dimension, point_count, 
										Log::algorithm_name(algorithm) + ":n:" + std::to_string(point_count));

		for (auto const &poly : vec) {
			_simplify(*poly, epsilon, config);
		}
	}


	log.emit();
}

void process_category(fs::path const &category_path, vector<vector<unique_ptr<Polyline>>> &data) {
	// Iterate over each numbered subdirectory
	for (const auto& entry : fs::directory_iterator(category_path)) {
		if (!entry.is_directory()) continue;

		vector<unique_ptr<Polyline>> polyline_subdir;

		// Collect all files in the numbered subdirectory
		for (const auto& file_entry : fs::directory_iterator(entry.path())) {
			if (file_entry.is_regular_file()) {
				polyline_subdir.push_back(Polyline::from_file(file_entry.path()));
			}
		}

		data.push_back(std::move(polyline_subdir));
	}
}

void measure_suite(std::filesystem::path path) {
	config.logger = std::make_optional<Log::PerformanceLogger>(Log::PerformanceLogger());
	auto &log = config.logger.value();
	config.output_visualization = false;
	vector<vector<unique_ptr<Polyline>>> w_polylines;
	vector<vector<unique_ptr<Polyline>>> n_polylines;
	fs::path well_behaved_path = fs::path(path) / "well-behaved";
	fs::path non_well_behaved_path = fs::path(path) / "non-well-behaved";
	process_category(well_behaved_path, w_polylines);
	process_category(non_well_behaved_path, n_polylines);

	use_algorithm<Log::Algorithm::SIMPLIFICATION_SIMPLE_MANHATTAN, Simplification::simplification_naive_manhattan>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_SIMPLE_CHEBYSHEV, Simplification::simplification_naive_chebyshev>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_SIMPLE_EUCLIDEAN, Simplification::simplification_naive_euclidean>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_SIMPLE_IMPLICIT_EUCLIDEAN, Simplification::simplification_naive_euclidean_implicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_SIMPLE_SEMIEXPLICIT_EUCLIDEAN, Simplification::simplification_naive_euclidean_semiexplicit>(w_polylines, n_polylines, log, _epsilon2);


	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_MANHATTAN, Simplification::simplification_advanced_manhattan_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_CHEBYSHEV, Simplification::simplification_advanced_chebyshev_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_EUCLIDEAN, Simplification::simplification_advanced_euclidean_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_SEMIEXPLICIT_EUCLIDEAN, Simplification::simplification_advanced_euclidean_semiexplicit>(w_polylines, n_polylines, log, _epsilon2);
}

void measure_suite_advanced(std::filesystem::path path) {
	config.logger = std::make_optional<Log::PerformanceLogger>(Log::PerformanceLogger());
	auto &log = config.logger.value();
	config.output_visualization = false;
	vector<vector<unique_ptr<Polyline>>> w_polylines;
	vector<vector<unique_ptr<Polyline>>> n_polylines;
	fs::path well_behaved_path = fs::path(path) / "well-behaved";
	fs::path non_well_behaved_path = fs::path(path) / "non-well-behaved";
	process_category(well_behaved_path, w_polylines);
	process_category(non_well_behaved_path, n_polylines);

	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_MANHATTAN, Simplification::simplification_advanced_manhattan_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_CHEBYSHEV, Simplification::simplification_advanced_chebyshev_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_EUCLIDEAN, Simplification::simplification_advanced_euclidean_explicit>(w_polylines, n_polylines, log, _epsilon);
	use_algorithm<Log::Algorithm::SIMPLIFICATION_ADVANCED_SEMIEXPLICIT_EUCLIDEAN, Simplification::simplification_advanced_euclidean_semiexplicit>(w_polylines, n_polylines, log, _epsilon2);
}

}
