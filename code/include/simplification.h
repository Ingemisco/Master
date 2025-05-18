#ifndef INCLUDE_INCLUDE_SIMPLIFICATION_H_
#define INCLUDE_INCLUDE_SIMPLIFICATION_H_

#include "datastructures.h"
#include <iostream>
#include <memory>
#include <vector>

namespace Simplification {

// vector of indices of the points included in the simplification in ascending order
typedef std::unique_ptr<std::vector<size_t>> Simplification;

Simplification simplification_naive_euclidean(DataStructures::Polyline &, float);
Simplification simplification_naive_manhattan(DataStructures::Polyline &, float);
Simplification simplification_naive_chebyshev(DataStructures::Polyline &, float);

Simplification simplification_naive_euclidean_implicit(DataStructures::Polyline &, float);
Simplification simplification_naive_euclidean_semiexplicit(DataStructures::Polyline &, float);

Simplification simplification_advanced_manhattan_explicit(DataStructures::Polyline &, float);
Simplification simplification_advanced_euclidean_explicit(DataStructures::Polyline &, float);
Simplification simplification_advanced_chebyshev_explicit(DataStructures::Polyline &, float);

Simplification simplification_advanced_euclidean_implicit(DataStructures::Polyline &, float);
Simplification simplification_advanced_euclidean_semiexplicit(DataStructures::Polyline &, float);


} // namespace Simplification

#endif // INCLUDE_INCLUDE_SIMPLIFICATION_H_
