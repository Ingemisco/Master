#define BOOST_TEST_MODULE GeneralTests

#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>

#include "distance.h"
#include "datastructures.h"

using DataStructures::Polyline;

using DataStructures::integer_exponentiation;
using DataStructures::unnormalized_euclidean_distance;
using DataStructures::euclidean_distance;
using DataStructures::manhattan_distance;
using DataStructures::chebyshev_distance;
using DataStructures::unnormalized_minkowski_distance;
using DataStructures::minkowski_distance;

BOOST_AUTO_TEST_CASE(IntegerExponentiation) {
	BOOST_CHECK_EQUAL(integer_exponentiation(2,5), 32);
	BOOST_CHECK_EQUAL(integer_exponentiation(3,4), 81);
	BOOST_CHECK_EQUAL(integer_exponentiation(1.5f, 2), 2.25f); // both numbers can be represented with no inaccuracies so equal check should work
}

BOOST_AUTO_TEST_CASE(PointDistances) {
	auto p = Polyline::from_file("test_data/p0");
	auto &polyline = *p;

	BOOST_CHECK_CLOSE(unnormalized_euclidean_distance(polyline, 0, 1), 5, 0);
	BOOST_CHECK_CLOSE(manhattan_distance(polyline, 0, 1), 3, 0);
	BOOST_CHECK_CLOSE(chebyshev_distance(polyline, 0, 1), 2, 0);

	BOOST_CHECK_CLOSE(euclidean_distance(polyline, 0, 1), std::sqrt(5), 1e-2);
	BOOST_CHECK_CLOSE(unnormalized_minkowski_distance(polyline, 0, 1, 2), 5, 0);
	BOOST_CHECK_CLOSE(unnormalized_minkowski_distance(polyline, 0, 1, 4), 17, 0);

	BOOST_CHECK_CLOSE(minkowski_distance(polyline, 0, 1, 2), sqrt(5), 1e-2);
	BOOST_CHECK_CLOSE(minkowski_distance(polyline, 0, 1, 4), std::pow(17, 0.25), 1e-2);

	auto [f, l] = DataStructures::solve_chebyshev(polyline, 0, 1, 2, 3);
	BOOST_CHECK_CLOSE(f, 0.5f, 1e-2);
	BOOST_CHECK_CLOSE(l, 1.0f, 1e-2);
}

