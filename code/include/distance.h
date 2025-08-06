#ifndef INCLUDE_INCLUDE_DISTANCE_H_
#define INCLUDE_INCLUDE_DISTANCE_H_

#include "datastructures.h"
#include <cmath>
#include <ostream>
#include <utility>

namespace DataStructures {
#define NO_POINT_REACHABLE -1

typedef float Distance(Polyline const &, size_t, size_t);

struct LRValue final {
	float a; 
	float d;

	bool operator==(LRValue const &) const;
};

struct FRValue final {
	float a; 
	float d;

	bool operator<=(FRValue const &) const;
	bool operator==(FRValue const &) const;
	bool operator<(FRValue const &) const;

	bool operator<=(LRValue const &) const;
};

inline std::ostream &operator<<(std::ostream &os, FRValue const &val) {
	os << "(" << val.a << ", " << val.d << ")";
	return os;
}

inline std::ostream &operator<<(std::ostream &os, LRValue const &val) {
	os << "(" << val.a << ", " << val.d << ")";
	return os;
}

typedef std::pair<float, float> ReachabilityData;
typedef std::pair<FRValue, LRValue> SEReachabilityData; // for semiexplicit euclidean

// first two points form a line segment and must not be the same, third point is
// the point from which we solve the equation
typedef ReachabilityData Solver(Polyline const &, size_t, size_t, size_t, float);

typedef float AltGodau(Polyline const &, size_t, size_t, float, size_t, size_t, float);

// used to designate that no point on a line segment is reachable.
// in ReachabilityData first and last will be set to this in this case
float constexpr EXPLICIT_UNREACHABLE = -1;

size_t constexpr IMPLICIT_UNREACHABLE = (size_t)-1;
size_t constexpr INDEX_UNREACHABLE = (size_t)-1;
size_t constexpr IMPLICIT_NEVER_REACHABLE = (size_t)-2;


std::pair<float const, float const> constexpr EMPTY_INTERVAL_EXPLICIT(EXPLICIT_UNREACHABLE, EXPLICIT_UNREACHABLE);
std::pair<FRValue const, LRValue const> constexpr EMPTY_INTERVAL_SEMIEXPLICIT(FRValue(-1,-1), LRValue(-1,-1));
std::pair<FRValue const, LRValue const> constexpr NONEMPTY_INTERVAL_SEMIEXPLICIT(FRValue(0.0f, 0.0f), LRValue(0.0f, -1.0f));

FRValue constexpr SEMIEXPLICIT_UNREACHABLE = EMPTY_INTERVAL_SEMIEXPLICIT.first;

float integer_exponentiation(float, int);
int integer_exponentiation(int, int);

float unnormalized_minkowski_distance(Polyline const &, size_t, size_t, int);
float minkowski_distance(Polyline const &, size_t, size_t, int);

float unnormalized_euclidean_distance(Polyline const &, size_t, size_t);
float euclidean_distance(Polyline const &, size_t, size_t);

float chebyshev_distance(Polyline const &, size_t, size_t);
float manhattan_distance(Polyline const &, size_t, size_t);

float alt_godau_manhattan(Polyline const &, size_t, size_t, float, size_t, size_t, float);
float alt_godau_euclidean(Polyline const &, size_t, size_t, float, size_t, size_t, float);

size_t alt_godau_euclidean_implicit(Polyline const &, size_t _, size_t, size_t, size_t, size_t, float);
FRValue alt_godau_euclidean_semiexplicit(Polyline const &, size_t, size_t, FRValue, size_t, size_t, float);

float alt_godau_chebyshev(Polyline const &, size_t, size_t, float, size_t, size_t, float);
size_t alt_godau_minkowski_implicit(Polyline const &, size_t, size_t, float, unsigned int);

ReachabilityData solve_manhattan(Polyline const &, size_t, size_t, size_t, float);
ReachabilityData solve_chebyshev(Polyline const &, size_t, size_t, size_t, float);
ReachabilityData solve_euclidean(Polyline const &, size_t, size_t, size_t, float);

bool solve_implicit_euclidean(Polyline const &, size_t, size_t, size_t, size_t, float);

SEReachabilityData solve_euclidean_se(Polyline const &, size_t, size_t, size_t, float);

bool is_line_reachable_euclidean(Polyline const &, size_t, size_t, size_t, float);

} // namespace DataStructures
#endif // INCLUDE_INCLUDE_DISTANCE_H_
