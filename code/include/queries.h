#ifndef INCLUDE_QUERIES_H
#define INCLUDE_QUERIES_H

#include "datastructures.h"
#include "simplification.h"
#include <cmath>
#include <memory>
#include <vector>

struct SimplificationQuerier final {
	std::vector<float> epsilons;
	std::vector<Simplification::Simplification> simplifications;
	size_t simplification_count;

	SimplificationQuerier(size_t n) : epsilons(n - 1, -1.0f), simplifications(n - 1),
		simplification_count(n-1) {
		for (size_t i = 0; i < n - 1; i++) {
			simplifications[i] = nullptr;
		}
	}

	SimplificationQuerier(SimplificationQuerier &) = delete;
	SimplificationQuerier(SimplificationQuerier &&) = delete;
	SimplificationQuerier &operator=(SimplificationQuerier &) = delete;
	SimplificationQuerier &operator=(SimplificationQuerier &&) = delete;
};

// all only for Euclidean distance 
Simplification::Simplification &simplify_query(SimplificationQuerier &, float);
std::unique_ptr<SimplificationQuerier> build_querier_simple(DataStructures::Polyline &);
std::unique_ptr<SimplificationQuerier> build_querier(DataStructures::Polyline &);

#endif // INCLUDE_QUERIES_H
