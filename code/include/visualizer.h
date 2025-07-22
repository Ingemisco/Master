#ifndef INCLUDE_INCLUDE_VISUALIZER_H_
#define INCLUDE_INCLUDE_VISUALIZER_H_

#include "datastructures.h"
#include "distance.h"
#include <vector>
namespace VisualizationLog {

enum class Distance { EUCLIDEAN, MANHATTAN, CHEBYSHEV, EUCLIDEAN_IMPLICIT, EUCLIDEAN_SEMIEXPLICIT };

struct VisualizationData final {
  size_t const k;
  size_t const i;
  size_t const j;
  size_t const i_;
  size_t const j_;
	union {
		float const t;                   // for explicit denotes value in [0, 1]
		size_t const r;                  // restriction for implicit 
		DataStructures::FRValue const f; // same as t but for semiexplicit 
	};

  VisualizationData(size_t, size_t, size_t, size_t, size_t, float);
  VisualizationData(size_t, size_t, size_t, size_t, size_t, size_t);
  VisualizationData(size_t, size_t, size_t, size_t, size_t, DataStructures::FRValue);
};

class VisualizationLogger final {
private:
  DataStructures::Polyline const &polyline;
  float const epsilon;
  std::vector<VisualizationData> data;
  Distance const distance;

public:
  VisualizationLogger(DataStructures::Polyline const &, float, Distance);

  void leave_range(size_t);
  void add_use(VisualizationData);
  void emit();
};
} // namespace VisualizationLog

#endif // INCLUDE_INCLUDE_VISUALIZER_H_
