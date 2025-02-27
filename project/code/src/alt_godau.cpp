// This file contains implementations for the algorithm from Alt and Godau to
// solve the decisional Fréchet distance problem where the second polyline is a
// single line segment. Different implementations are for different distance
// measures (Manhattan, Euclidean, Chebyshev) as well other differences in
// implementation

// explicit versions return the earliest reachable point on the last line
// segment of the polyline if it exists and otherwise UNREACHABLE.
// Implicit versions only return the index of the point whose distance to the
// last line segments bounds the first reachable point and UNREACHABLE if no
// such point exists.

#include "config.h"
#include "datastructures.h"
#include "distance.h"
#include <algorithm>
#include <cmath>

namespace DataStructures {

float alt_godau_manhattan(Polyline const &polyline, Point const &line_start,
                          Point const &line_end, float epsilon);

float alt_godau_euclidean(Polyline const &polyline, Point const &line_start,
                          Point const &line_end, float epsilon) {
#if DEBUG
  assert_compatible_points(line_start, line_end);
  assert_compatible_points(polyline.get_point(0), line_end);
#endif

  if (unnormalized_euclidean_distance(polyline.get_point(0), line_start) >
      epsilon * epsilon) {
    return UNREACHABLE;
  }

  float first_reachable = 0;

  // first points already matched so start with 1 instead of 0
  for (unsigned int i = 1; i < polyline.point_count; i++) {
    auto data =
        solve_euclidean(line_start, line_end, polyline.get_point(i), epsilon);
    if (data.first == UNREACHABLE || data.last < first_reachable) {
      return UNREACHABLE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }

  return first_reachable;
}

size_t alt_godau_euclidean_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

float alt_godau_chebyshev(Polyline const &polyline, Point const &line_start,
                          Point const &line_end, float epsilon);

size_t alt_godau_minkowski_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

} // namespace DataStructures
