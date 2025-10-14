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

// logically equivalent to any simpification algorithm with epsilon = 0, but more stable and linear runtime
Simplification simplification_filter_colinear(Polyline const &, AlgorithmConfiguration &);

Simplification simplification_imai_iri_euclidean(Polyline const &, float, AlgorithmConfiguration &);

#if __has_include(<generator>)
Simplification simplification_global_imai_iri_euclidean(Polyline const &, float, AlgorithmConfiguration &);
#endif

} // namespace Simplification

#endif // INCLUDE_INCLUDE_SIMPLIFICATION_H_
