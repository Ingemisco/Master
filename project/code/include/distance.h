#ifndef INCLUDE_INCLUDE_DISTANCE_H_
#define INCLUDE_INCLUDE_DISTANCE_H_

#include "datastructures.h"

namespace DataStructures {

typedef float Distance(Point &, Point &);

float euclidean_distance(Point &, Point &);

bool frechet_distance_at_most(DataStructures::Polyline const &,
                              DataStructures::Polyline const &, float,
                              Distance = euclidean_distance);

float frechet_distance(DataStructures::Polyline const &,
                       DataStructures::Polyline const &,
                       Distance = euclidean_distance);

} // namespace DataStructures
#endif // INCLUDE_INCLUDE_DISTANCE_H_
