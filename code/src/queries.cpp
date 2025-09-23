#include "queries.h"
#include "datastructures.h"
#include "distance.h"
#include "simplification.h"
#include "log.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <string>

using DataStructures::Polyline, DataStructures::unnormalized_euclidean_distance;

Simplification::Simplification &simplify_query(SimplificationQuerier &, float);
AlgorithmConfiguration empty_config = {
		.output_visualization = false,
		.logger = std::nullopt,
};

void static _binary_search_build_ds(Polyline &polyline, float *events, SimplificationQuerier &ds,
																		size_t event_count, size_t min_to_set, size_t max_to_set) {
	if (event_count == 0 || min_to_set > max_to_set) {
		return;
	} 

	size_t mid_index = event_count / 2;
	// float const epsilon = std::sqrt(events[mid_index]); // events are all squared
	float const epsilon = events[mid_index]; // events are all squared
	
	Simplification::Simplification simplification;
	if (epsilon > 1e-5) {
		// auto simplification = Simplification::simplification_advanced_euclidean_explicit(polyline, epsilon + 1e-7, empty_config);
		simplification = Simplification::simplification_advanced_euclidean_semiexplicit(polyline, epsilon + 1e-7, empty_config);
	} else {
		// if epsilon is zero we only need to remove colinear points (somehow this does not work with the regular simplifications algorithms??)
		simplification = Simplification::simplification_filter_colinear(polyline, empty_config);
	}
	// this "-1" is to get the size of the simplification, i.e., the amount of line segments instead of amount of points 
	size_t const simplification_size = simplification->size() - 1;

	if (simplification_size <= max_to_set) {
		// this "-1" is to get 0 indexing, as the simplification sizes are from 1 to n - 1, not 0 to n - 2
		ds.simplifications[simplification_size - 1] = std::move(simplification);
		ds.epsilons[simplification_size - 1] = epsilon;
	}

	_binary_search_build_ds(polyline, events, ds, mid_index, simplification_size, max_to_set);
	_binary_search_build_ds(polyline, events + mid_index + 1, ds, event_count - mid_index - 1, min_to_set, simplification_size - 1);
}


std::unique_ptr<SimplificationQuerier> build_querier_simple(Polyline &polyline, AlgorithmConfiguration &config) {
	auto start = std::chrono::high_resolution_clock::now();
	size_t const n = polyline.point_count;

	// Create array of events 
	float *events = new float[n * n * (n-1) * (n + 1) / 4];
	size_t event_count = 0;
	for (size_t i = 0; i < n - 1; i++) {
		for (size_t j = i + 1; j < n; j++) {
			float const dij = unnormalized_euclidean_distance(polyline, i, j);
			for (size_t k=0; k < n; k++) {
				float const dik = unnormalized_euclidean_distance(polyline, i, k);
				float event = dik;
				event = std::min(unnormalized_euclidean_distance(polyline, j, k), event);
				float const scalar_prod = scalar_product(polyline, i, j, k); 
				float const t = scalar_prod / dij;
				if (0 <= t && t <= 1) {
					// float inaccuracies can make this number negative
					event = std::min(event, std::max(dik - scalar_prod * t, 0.0f) );
				}
				events[event_count++] = event;
			}

			for (size_t u = 0; u < n - 1; u++) {
				for (size_t v = i + 1; v < n; v++) {   
					float const scalar_prod = scalar_product(polyline, u, v, j, i);
					if (scalar_prod == 0) {
					  continue;
					}
					float const diu = unnormalized_euclidean_distance(polyline, i, u);
					float const div = unnormalized_euclidean_distance(polyline, i, v);
					float const t = (diu - div) / (2 * scalar_prod);
					if (t > 1 || t < 0) {
					  continue;
					}
					events[event_count++] = std::max(diu + t * (2 * scalar_product(polyline, i, j, u) + t * dij), 0.0f);
				}
			}
		}
	}

	std::sort(events, events + event_count);
	size_t k = 0;
	for (size_t i = 1; i < event_count; i++) {
		if (1e-06 < events[i] - events[k]) {
			k++;
			events[k] = events[i];
		}
	}
	event_count = k+1;

  auto datastructure = std::make_unique<SimplificationQuerier>(n);
	auto temp_datastructure = SimplificationQuerier(n);
	_binary_search_build_ds(polyline, events, temp_datastructure, event_count, 1, n - 1);
	size_t index = 0;
	for (size_t i = 0; i < n - 1; i++) {
		if (temp_datastructure.epsilons[i] != -1.0) {
			datastructure->epsilons[index] = temp_datastructure.epsilons[i];
			datastructure->simplifications[index++] = std::move(temp_datastructure.simplifications[i]);
		}
	}

	for (size_t i = 0; i < index; i++) {
		std::cout << "simplification size: " << datastructure->simplifications[i]->size() << ", epsilon: " << datastructure->epsilons[i] << ", Simplification: ";
		for (auto v: *datastructure->simplifications[i]) {
			std::cout << v << ", ";
		}
		std::cout << std::endl;
	}
	delete[] events;


	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			config.logger.value().add_data(0, end - start, "");
		}
	}


	return datastructure;
}

void SimplificationQuerier::save_datastructure_to_file(std::filesystem::path path) {
  std::ofstream output;
  output.open(path);
  if (!output.is_open()) {
    std::cerr << "Could not write to file " << path << std::endl;
    exit(1);
  }


	for (size_t i = 0; i < this->epsilons.size(); i++) {
		if (!simplifications[i]) {
			// the simplification size of i is always skipped in favor for a smaller simplification no matter the epsilon so we skip it here
			continue; 
		}
		output << epsilons[i] << " " << simplifications[i]->size();
		for (auto v: *simplifications[i]) {
			output << " " << v;
		}
		output << '\n';
	}

	output.close();
}


std::unique_ptr<SimplificationQuerier> SimplificationQuerier::from_file(std::filesystem::path path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("Failed to open file: " + path.string());
  }

  auto ds = std::make_unique<SimplificationQuerier>();
	std::string line;
	float epsilon;
	size_t simplification_size;
	size_t point;
	while (std::getline(input, line)) {
		std::istringstream line_stream(line);
		line_stream >> epsilon >> simplification_size;
		auto polyline = std::make_unique<std::vector<size_t>>();

		polyline->reserve(simplification_size);
		while ((line_stream >> point)) {
			polyline->push_back(point);
		}
		
		ds->simplifications.push_back(std::move(polyline));
		ds->epsilons.push_back(epsilon);
	}
	return ds;
}



std::unique_ptr<SimplificationQuerier> build_querier(Polyline &);











void SimplificationQuerier::print() {
	for (unsigned int i = 0; i < epsilons.size(); i++) {
		std::cout << epsilons[i] << ": ";
		for (auto &v : *simplifications[i]) {
			std::cout << v << " ";
		}
		std::cout << std::endl;
	}
}
