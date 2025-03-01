#ifndef INCLUDE_INCLUDE_SIMPLIFICATION_H_
#define INCLUDE_INCLUDE_SIMPLIFICATION_H_

#include "datastructures.h"
#include <memory>
#include <vector>

namespace Simplification {

// vector of indices of the points included in the simplification in ascending
// order
typedef std::unique_ptr<std::vector<size_t>> Simplification;

Simplification simplification_naive_euclidean(DataStructures::Polyline &,
                                              float);

} // namespace Simplification

#endif // INCLUDE_INCLUDE_SIMPLIFICATION_H_
