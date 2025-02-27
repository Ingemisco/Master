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

template <ReachabilityData _solver(Point const &, Point const &, Point const &,
                                   float)>
static inline float _alt_godau_main(PolylineRange polyline, LineSegment line,
                                    float epsilon) {
  float first_reachable = 0;

  // first points already matched so start with 1 instead of 0
  for (unsigned int i = polyline.start_point_index + 1;
       i <= polyline.end_point_index; i++) {
    auto data =
        _solver(line.start, line.end, polyline.polyline.get_point(i), epsilon);
    if (data.first == UNREACHABLE || data.last < first_reachable) {
      return UNREACHABLE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }
  return first_reachable;
}

float alt_godau_manhattan(PolylineRange polyline, LineSegment line,
                          float epsilon);

float alt_godau_euclidean(PolylineRange polyline, LineSegment line,
                          float epsilon) {
#if DEBUG
  assert_compatible_points(polyline.polyline.get_point(0), line.start);
#endif

  // compute distance of first point on polyline (p[start] + t(p[start + 1] -
  // p[start])) to first point of line segment and compare if distance greater
  // then epsilon, if so, the line segment is too far away
  float initial_dist = 0;
  for (unsigned int i = 0; i < line.start.dimension; i++) {
    float const coord =
        polyline.polyline[polyline.start_point_index, i] +
        polyline.start_point_offset *
            (polyline.polyline[polyline.start_point_index + 1, i] -
             polyline.polyline[polyline.start_point_index, i]) -
        line.start[i];
    initial_dist += coord * coord;
  }
  if (initial_dist > epsilon * epsilon) {
    return UNREACHABLE;
  }

  return _alt_godau_main<solve_euclidean>(polyline, line, epsilon);
}

size_t alt_godau_euclidean_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

float alt_godau_chebyshev(PolylineRange polyline, LineSegment line,
                          float epsilon);

size_t alt_godau_minkowski_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

} // namespace DataStructures
