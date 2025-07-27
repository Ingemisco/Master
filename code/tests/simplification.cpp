#include "simplification.h"
#include <optional>
#define BOOST_TEST_MODULE SimplificationTests

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/unit_test.hpp>

#include "datastructures.h"
#include "log.h"

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

Log::AlgorithmConfiguration config {
	.output_visualization = false,
	.logger = std::nullopt,
};

// Tests some hand-computed polylines for their simplification sizes 
BOOST_AUTO_TEST_CASE(SimplificationSizes) {
	auto const p = Polyline::from_file("test_data/p0");
	auto const &polyline = *p;

}



static inline void compare_semiexplicit(std::filesystem::path polyline_file_path, float epsilon) {
	auto const p = Polyline::from_file(polyline_file_path);
	auto const &polyline = *p;

	size_t const euclidean_size_advanced = simplification_advanced_euclidean_explicit(polyline, epsilon, config)->size();

	float const epsilon2 = epsilon * epsilon;
	size_t const euclidean_size_naive_semi = simplification_naive_euclidean_semiexplicit(polyline, epsilon2, config)->size();
	size_t const euclidean_size_advanced_semi = simplification_advanced_euclidean_semiexplicit(polyline, epsilon2, config)->size();

	BOOST_CHECK_EQUAL(euclidean_size_naive_semi, euclidean_size_advanced_semi);
	BOOST_CHECK_EQUAL(euclidean_size_advanced_semi, euclidean_size_advanced);
}

// Tests that the semiexplicit versions with epsilon squared and the regular ones with epsilon yield the same size 
BOOST_AUTO_TEST_CASE(EpsilonVSEpsilonSquared) {
	compare_semiexplicit("test_data/p0", 1);
	compare_semiexplicit("test_data/p0", 3);
	compare_semiexplicit("test_data/p0", 3);
	compare_semiexplicit("test_data/p0", 2.5);
	compare_semiexplicit("test_data/p0", 1.72);
}











static inline void compare_sizes(std::filesystem::path polyline_file_path, float epsilon) {
	auto const p = Polyline::from_file(polyline_file_path);
	auto const &polyline = *p;
  
	size_t const manhattan_size_naive    = simplification_naive_manhattan(polyline, epsilon, config)->size();
	size_t const manhattan_size_advanced = simplification_advanced_manhattan_explicit(polyline, epsilon, config)->size();
	size_t const chebyshev_size_naive    = simplification_naive_chebyshev(polyline, epsilon, config)->size();
	size_t const chebyshev_size_advanced = simplification_advanced_chebyshev_explicit(polyline, epsilon, config)->size();
	size_t const euclidean_size_naive    = simplification_naive_euclidean(polyline, epsilon, config)->size();
	size_t const euclidean_size_advanced = simplification_advanced_euclidean_explicit(polyline, epsilon, config)->size();
  
	BOOST_CHECK_EQUAL(euclidean_size_advanced, euclidean_size_naive);
	BOOST_CHECK_EQUAL(manhattan_size_advanced, manhattan_size_naive);
	BOOST_CHECK_EQUAL(chebyshev_size_advanced, chebyshev_size_naive);
  
	BOOST_CHECK_LE(euclidean_size_advanced, manhattan_size_advanced);
	BOOST_CHECK_LE(chebyshev_size_advanced, euclidean_size_advanced);
}


// Compares the sizes of the simplifications. It must always hold for a fixed epsilon that the Manhattan is the largest and Chebyshev is the smallest one. 
BOOST_AUTO_TEST_CASE(RelativeSizes) {
	compare_sizes("test_data/p0", 1);
	compare_sizes("test_data/p0", 2);
	compare_sizes("test_data/p0", 3);
	compare_sizes("test_data/p0", 4);
	compare_sizes("test_data/p0", 2.5);
}

