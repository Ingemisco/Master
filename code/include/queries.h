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

	SimplificationQuerier() {}

	SimplificationQuerier(size_t n) : epsilons(n - 1, -1.0f), simplifications(n - 1) {
		for (size_t i = 0; i < n - 1; i++) {
			simplifications[i] = nullptr;
		}
	}

	SimplificationQuerier(SimplificationQuerier &) = delete;
	SimplificationQuerier(SimplificationQuerier &&) = delete;
	SimplificationQuerier &operator=(SimplificationQuerier &) = delete;
	SimplificationQuerier &operator=(SimplificationQuerier &&) = delete;

	void save_datastructure_to_file(std::filesystem::path);
	static std::unique_ptr<SimplificationQuerier> from_file(std::filesystem::path);

	void print();
};

// all only for Euclidean distance 
Simplification::Simplification &simplify_query(SimplificationQuerier &, float);
std::unique_ptr<SimplificationQuerier> build_querier_simple(DataStructures::Polyline &, AlgorithmConfiguration &);
std::unique_ptr<SimplificationQuerier> build_querier(DataStructures::Polyline &);


#endif // INCLUDE_QUERIES_H
