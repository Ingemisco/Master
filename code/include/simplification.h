#ifndef INCLUDE_INCLUDE_SIMPLIFICATION_H_
#define INCLUDE_INCLUDE_SIMPLIFICATION_H_

#include "datastructures.h"
#include <memory>
#include <vector>

namespace Log {
	struct AlgorithmConfiguration;
}
using DataStructures::Polyline;
using Log::AlgorithmConfiguration;

namespace Simplification {


// vector of indices of the points included in the simplification in ascending order
typedef std::unique_ptr<std::vector<size_t>> Simplification;

Simplification simplification_naive_euclidean(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_naive_manhattan(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_naive_chebyshev(Polyline const &, float, AlgorithmConfiguration &);

Simplification simplification_naive_euclidean_implicit(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_naive_euclidean_semiexplicit(Polyline const &, float, AlgorithmConfiguration &);

Simplification simplification_advanced_manhattan_explicit(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_advanced_euclidean_explicit(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_advanced_chebyshev_explicit(Polyline const &, float, AlgorithmConfiguration &);

Simplification simplification_advanced_euclidean_implicit(Polyline const &, float, AlgorithmConfiguration &);
Simplification simplification_advanced_euclidean_semiexplicit(Polyline const &, float, AlgorithmConfiguration &);


} // namespace Simplification

#endif // INCLUDE_INCLUDE_SIMPLIFICATION_H_
