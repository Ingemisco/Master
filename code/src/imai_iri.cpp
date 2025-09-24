#include "distance.h"
#include "log.h"
#include "simplification.h"
#include <cassert>

namespace Simplification {

// tests if shortcut from start to end is valid local shortcut. start and end are inclusive
static inline bool is_local_shortcut_euclidean(Polyline const &polyline, size_t start, size_t end, float epsilon) {
  float first_reachable = 0;
  // first points already matched so start with 1 instead of 0
  // do not need to go to end of polyline so last point exluded
  for (unsigned int index = start + 1; index < end; index++) {
    auto data = DataStructures::solve_euclidean(polyline, start, end, index, epsilon);
    if (data == DataStructures::EMPTY_INTERVAL_EXPLICIT || !(first_reachable <= data.second)) {
      return false;
    }

    first_reachable = std::max(first_reachable, data.first);
  }
  return true;
}

Simplification simplification_imai_iri_euclidean(Polyline const &polyline, float epsilon, AlgorithmConfiguration &config) {
	auto start = std::chrono::high_resolution_clock::now();
	unsigned int *sizes = new unsigned int[polyline.point_count];
	unsigned int *path  = new unsigned int[polyline.point_count];
	sizes[0] = 0;
	for (unsigned int i = 1; i < polyline.point_count; i++) {
		sizes[i] = sizes[i-1] + 1;
		path[i] = i - 1;
		for (unsigned int j = 0; j < i; j++) {
			if (sizes[j] + 1 < sizes[i] && is_local_shortcut_euclidean(polyline, j, i, epsilon)) {
				sizes[i] = sizes[j] + 1;
				path[i] = j;
			}
		}
	}

	unsigned int k = polyline.point_count - 1;
	unsigned int index = sizes[k];
	Simplification result = std::make_unique<std::vector<size_t>>(index + 1);
	while (index > 0) {
		assert(k > 0 && k < polyline.point_count);
		(*result)[index] = k;
		k = path[k];
		index--;
	}
	
	(*result)[0] = 0;


	if (config.logger.has_value()) {
		auto end = std::chrono::high_resolution_clock::now();
		if (config.logger.has_value()) {
			// last entry is number how many shortcuts are needed, + 1 because point is one more
			config.logger.value().add_data(sizes[polyline.point_count-1] + 1, end - start, "");
		}
	}

	delete[] sizes;
	delete[] path;
	return result;
}
}
