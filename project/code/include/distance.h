#ifndef INCLUDE_INCLUDE_DISTANCE_H_
#define INCLUDE_INCLUDE_DISTANCE_H_

#include "datastructures.h"

namespace DataStructures {

float frechet_distance(DataStructures::Polyline const &,
                       DataStructures::Polyline const &);

bool frechet_distance_at_most(DataStructures::Polyline const &,
                       DataStructures::Polyline const &, float);

}
#endif // INCLUDE_INCLUDE_DISTANCE_H_
