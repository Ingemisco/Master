#include "distance.h"
#include "datastructures.h"
#include <cmath>

namespace DataStructures {

typedef struct {
  float min_line1;
  float min_line2;
} ReachabilityTable;

float euclidean_distance(Point &point1, Point &point2) {
  assert_compatible_points(point1, point2);
  float result = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    result += (point1[i] - point2[i]) * (point1[i] - point2[i]);
  }
  return std::sqrt(result);
}

// TODO: Use algorithm from Alt & Godau for the decision variant
bool frechet_distance_at_most(DataStructures::Polyline const &polyline1,
                              DataStructures::Polyline const &polyline2,
                              [[maybe_unused]] float epsilon,
                              [[maybe_unused]] Distance distance) {
  assert_compatible_polylines(polyline1, polyline2);

  Polyline const *smaller_line = &polyline2;
  Polyline const *larger_line = &polyline1;
  if (polyline1.point_count < polyline2.point_count) {
    smaller_line = &polyline1;
    larger_line = &polyline1;
  }
  size_t const table_size = smaller_line->point_count;

  ReachabilityTable *reachable_table = new ReachabilityTable[table_size];

  // TODO: complete this
  // MAYBE TODO: optimize space to O(min(m, n)) ? / think if it works
  for (unsigned int i = 0; i < polyline1.point_count; i++) {
    for (unsigned int j = 0; j < polyline2.point_count; j++) {
    }
  }

  delete[] reachable_table;
  return false;
}

// TODO: Use algorithm from Alt & Godau
float frechet_distance(Polyline const &polyline1, Polyline const &polyline2,
                       [[maybe_unused]] Distance distance) {
  assert_compatible_polylines(polyline1, polyline2);
  return 0;
}
} // namespace DataStructures
