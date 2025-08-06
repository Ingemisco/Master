#define BOOST_TEST_MODULE AltGodauTests

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/unit_test.hpp>
#include "simplification.h"
#include <random>
#include <cstdlib>
#include "generators.h"

#include "distance.h"
#include "datastructures.h"

constexpr unsigned long RANDOM_SEED = 0; // needs to be deterministic
constexpr unsigned int RANDOM_TEST_CASES_COUNT = 1000000;
constexpr Dimension dimension = 2;

using DataStructures::Polyline;

static inline void test_implicit_explicit(Polyline const &polyline, size_t j_, size_t j, size_t i_, size_t i, size_t r_, float epsilon)  {
	auto const [t_, _unused] = DataStructures::solve_euclidean(polyline, j_, j_ + 1, r_, epsilon);
	auto const [i_reachable, _unused2] = DataStructures::solve_euclidean(polyline, j_, j_ + 1, i_, epsilon);
	if (t_ == DataStructures::EXPLICIT_UNREACHABLE || i_reachable == DataStructures::EXPLICIT_UNREACHABLE) {
		// comparison would not make any sense if there is no solution
		return;
	}
	auto const r = DataStructures::alt_godau_euclidean_implicit(polyline, j_, j, i_, i, r_, epsilon * epsilon);
	
	auto const t  = DataStructures::alt_godau_euclidean(polyline, j_, j, t_, i_, i, epsilon);

	float t_compare;
	if (r != DataStructures::IMPLICIT_UNREACHABLE) {
		t_compare = DataStructures::solve_euclidean(polyline, j, j + 1, r, epsilon).first;
	} else {
		t_compare = DataStructures::EXPLICIT_UNREACHABLE;
	}

	BOOST_CHECK_CLOSE(t, t_compare, 1e-3);
}

BOOST_AUTO_TEST_CASE(ValuesTest) {

}

BOOST_AUTO_TEST_CASE(ExplicitImplicit) {
	std::mt19937 gen(RANDOM_SEED);
	for (unsigned int i = 0; i < RANDOM_TEST_CASES_COUNT; i++) {
		auto p = DataGeneration::make_polyline(10, dimension, 2, 10, std::numbers::pi, gen);
  
		test_implicit_explicit(*p, 4, 4, 2, 4, 0, 2);
		test_implicit_explicit(*p, 4, 6, 2, 4, 0, 2);
		test_implicit_explicit(*p, 4, 6, 2, 4, 1, 2);
		test_implicit_explicit(*p, 4, 6, 2, 4, 2, 2);
	}

	auto p = Polyline::from_file("test_data/fail.poly");
	test_implicit_explicit(*p, 7, 11, 8, 12, 12, 2);
}

