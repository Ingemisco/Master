#ifndef INCLUDE_INCLUDE_DISTANCE_H_
#define INCLUDE_INCLUDE_DISTANCE_H_

#include "datastructures.h"

namespace DataStructures {
#define NO_POINT_REACHABLE -1

typedef float Distance(Point const &, Point const &);
typedef struct {
  float first;
  float last;
} ReachabilityData;

typedef float Distance(Point const &, Point const &);

// first two points form a line segment and must not be the same, third point is
// the point from which we solve the equation
typedef ReachabilityData Solver(Point const &, Point const &, Point const &,
                                float);

typedef float AltGodau(PolylineRange, LineSegment, float);

// used to designate that no point on a line segment is reachable.
// in ReachabilityData first and last will be set to this in this case
float const EXPLICIT_UNREACHABLE = -1;

size_t const IMPLICIT_UNREACHABLE = (size_t)-1;
size_t const EXPLICIT_INDEX_UNREACHABLE = (size_t)-1;
size_t const IMPLICIT_NEVER_REACHABLE = (size_t)-2;

float integer_exponentiation(float, int);
int integer_exponentiation(int, int);

float unnormalized_minkowski_distance(Point const &, Point const &, int);
float minkowski_distance(Point const &, Point const &, int);

float unnormalized_euclidean_distance(Point const &, Point const &);
float euclidean_distance(Point const &, Point const &);

float chebyshev_distance(Point const &, Point const &);
float manhattan_distance(Point const &, Point const &);

float alt_godau_manhattan(PolylineRange, LineSegment, float);
float alt_godau_euclidean(PolylineRange, LineSegment, float);

size_t alt_godau_euclidean_implicit(Polyline const &, size_t _, size_t, size_t,
                                    size_t, size_t, float);

float alt_godau_chebyshev(PolylineRange, LineSegment, float);
size_t alt_godau_minkowski_implicit(Polyline const &, Point const &,
                                    Point const &, float, unsigned int);

ReachabilityData solve_manhattan(Point const &, Point const &, Point const &,
                                 float);
ReachabilityData solve_chebyshev(Point const &, Point const &, Point const &,
                                 float);
ReachabilityData solve_euclidean(Point const &, Point const &, Point const &,
                                 float);

bool solve_implicit_euclidean(LineSegment, Point const &, Point const &, float);

size_t solve_implicit_euclidean_in(LineSegment, Point const &, Point const &,
                                   float);

bool is_line_reachable_euclidean(LineSegment, Point const &, float);

} // namespace DataStructures
#endif // INCLUDE_INCLUDE_DISTANCE_H_
