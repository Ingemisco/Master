#ifndef INCLUDE_INCLUDE_VISUALIZER_H_
#define INCLUDE_INCLUDE_VISUALIZER_H_

#include "datastructures.h"
#include <vector>
namespace VisualizationLog {

enum class Distance { EUCLIDEAN, MANHATTAN, CHEBYSHEV, EUCLIDEAN_IMPLICIT };

struct VisualizationData final {
  size_t const k;
  size_t const i;
  size_t const j;
  size_t const i_;
  size_t const j_;
  float const t;

  VisualizationData(size_t, size_t, size_t, size_t, size_t, float);
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
