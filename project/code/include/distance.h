#ifndef INCLUDE_INCLUDE_DISTANCE_H_
#define INCLUDE_INCLUDE_DISTANCE_H_

#include "datastructures.h"

namespace DataStructures {

typedef float Distance(Point const &, Point const &);
typedef struct {
  float first;
  float last;
} ReachabilityData;

// first two points form a line segment and must not be the same, third point is
// the point from which we solve the equation
typedef ReachabilityData Solver(Point const &, Point const &, Point const &,
                                float);

// used to designate that no point on a line segment is reachable.
// in ReachabilityData first and last will be set to this in this case
float const UNREACHABLE = -1;

float integer_exponentiation(float, int);
int integer_exponentiation(int, int);

float unnormalized_minkowski_distance(Point const &, Point const &, int);
float minkowski_distance(Point const &, Point const &, int);

float unnormalized_euclidean_distance(Point const &, Point const &);
float euclidean_distance(Point const &, Point const &);

float maximum_norm_distance(Point const &, Point const &);
float manhattan_distance(Point const &, Point const &);

bool frechet_distance_decision(DataStructures::Polyline const &,
                               DataStructures::Polyline const &, float,
                               Distance = euclidean_distance);

bool frechet_distance_decision_eucliean_sqrt(DataStructures::Polyline const &,
                                             DataStructures::Polyline const &,
                                             float);
bool frechet_distance_decision_eucliean_nosqrt(DataStructures::Polyline const &,
                                               DataStructures::Polyline const &,
                                               float);

float frechet_distance(DataStructures::Polyline const &,
                       DataStructures::Polyline const &,
                       Distance = euclidean_distance);

ReachabilityData solve_manhattan(Point const &, Point const &, Point const &,
                                 float);
ReachabilityData solve_maximum(Point const &, Point const &, Point const &,
                               float);
ReachabilityData solve_euclidean(Point const &, Point const &, Point const &,
                                 float);

} // namespace DataStructures
#endif // INCLUDE_INCLUDE_DISTANCE_H_
