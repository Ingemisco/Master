#define BOOST_TEST_MODULE SimplificationTests

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/unit_test.hpp>
#include <optional>
#include <random>

#include "datastructures.h"
#include "log.h"
#include "simplification.h"
#include "generators.h"
#include "global.h"

using DataStructures::Polyline;
using Simplification::simplification_advanced_euclidean_explicit;
using Simplification::simplification_advanced_chebyshev_explicit;
using Simplification::simplification_advanced_manhattan_explicit;
using Simplification::simplification_advanced_euclidean_semiexplicit;
using Simplification::simplification_naive_euclidean;
using Simplification::simplification_naive_chebyshev;
using Simplification::simplification_naive_manhattan;
using Simplification::simplification_naive_euclidean_semiexplicit;
using Simplification::simplification_naive_euclidean_implicit;

constexpr unsigned long RANDOM_SEED = 0; // needs to be deterministic
constexpr unsigned int RANDOM_TEST_CASES_COUNT = 20;
constexpr Dimension dimension = 2;


template <void test_function(Polyline const &, float)>
static inline void perform_on_polylines() {
	// auto const p0 = Polyline::from_file("test_data/p0");
	auto const p0 = Polyline::from_file("test_data/fail.poly");
	auto const &polyline0 = *p0;
	test_function(polyline0, 1);
	test_function(polyline0, 3);
	test_function(polyline0, 3);
	test_function(polyline0, 2);
	test_function(polyline0, 2.5);
	test_function(polyline0, 1.72);

	std::mt19937 gen(RANDOM_SEED);
	for (unsigned int i = 0; i < RANDOM_TEST_CASES_COUNT; i++) {
		auto p = DataGeneration::make_polyline(10 + 5 * i, dimension, 2, 10, std::numbers::pi, gen);
		// auto p = DataGeneration::make_polyline(10 + 5 * i, dimension, 2, 10, 175 * std::numbers::pi / 180.0, gen);
		test_function(*p, 1);
		test_function(*p, 2);
	}
}

Log::AlgorithmConfiguration config {
	.output_visualization = false,
	.logger = std::nullopt,
};

static inline void compare_semiexplicit(Polyline const &polyline, float epsilon) {
	size_t const euclidean_size_advanced = simplification_advanced_euclidean_explicit(polyline, epsilon, config)->size();

	float const epsilon2 = epsilon * epsilon;
	size_t const euclidean_size_naive_semi = simplification_naive_euclidean_semiexplicit(polyline, epsilon2, config)->size();
	size_t const euclidean_size_advanced_semi = simplification_advanced_euclidean_semiexplicit(polyline, epsilon2, config)->size();


	BOOST_CHECK_EQUAL(euclidean_size_naive_semi, euclidean_size_advanced_semi);
	BOOST_CHECK_EQUAL(euclidean_size_advanced_semi, euclidean_size_advanced);
}

// Tests that the semiexplicit versions with epsilon squared and the regular ones with epsilon yield the same size 
BOOST_AUTO_TEST_CASE(EpsilonVSEpsilonSquared) {
	perform_on_polylines<compare_semiexplicit>();
}

static inline void compare_sizes(Polyline const &polyline, float epsilon) {
	size_t const manhattan_size_naive    = simplification_naive_manhattan(polyline, epsilon, config)->size();
	size_t const manhattan_size_advanced = simplification_advanced_manhattan_explicit(polyline, epsilon, config)->size();
	size_t const chebyshev_size_naive    = simplification_naive_chebyshev(polyline, epsilon, config)->size();
	size_t const chebyshev_size_advanced = simplification_advanced_chebyshev_explicit(polyline, epsilon, config)->size();
	size_t const euclidean_size_naive    = simplification_naive_euclidean(polyline, epsilon, config)->size();
	size_t const euclidean_size_advanced = simplification_advanced_euclidean_explicit(polyline, epsilon, config)->size();
	size_t const euclidean_size_implicit = simplification_naive_euclidean_implicit(polyline, epsilon, config)->size();

	size_t const euclidean_size_local    = Simplification::simplification_imai_iri_euclidean(polyline, epsilon, config)->size();
	size_t const euclidean_sizeg_heuristic = Simplification::simplification_global_imai_iri_euclidean(polyline, epsilon, config)->size();
  
	BOOST_CHECK_EQUAL(euclidean_size_advanced, euclidean_size_naive);
	BOOST_CHECK_EQUAL(euclidean_size_advanced, euclidean_size_implicit);
	BOOST_CHECK_EQUAL(manhattan_size_advanced, manhattan_size_naive);
	BOOST_CHECK_EQUAL(chebyshev_size_advanced, chebyshev_size_naive);
  
	BOOST_CHECK_LE(euclidean_size_advanced, manhattan_size_advanced);
	BOOST_CHECK_LE(chebyshev_size_advanced, euclidean_size_advanced);

	BOOST_CHECK_LE(euclidean_size_advanced, euclidean_sizeg_heuristic);
	BOOST_CHECK_LE(euclidean_sizeg_heuristic, euclidean_size_local);
}

// Compares the sizes of the simplifications. It must always hold for a fixed epsilon that the Manhattan is the largest and Chebyshev is the smallest one. 
BOOST_AUTO_TEST_CASE(RelativeSizes) {
	perform_on_polylines<compare_sizes>();
}

