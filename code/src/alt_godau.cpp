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
#include <cstdlib>

namespace DataStructures {

template <ReachabilityData _solver(Point const &, Point const &, Point const &,
                                   float)>
static inline float _alt_godau_main(PolylineRange &polyline, LineSegment &line,
                                    float epsilon) {
  if (polyline.start_point_index + 1 == polyline.end_point_index) {
    auto data = _solver(polyline.polyline.get_point(polyline.start_point_index),
                        polyline.polyline.get_point(polyline.end_point_index),
                        line.end, epsilon);
    if (data.first == UNREACHABLE || polyline.start_point_offset > data.last) {
      return UNREACHABLE;
    }
    return std::max(data.first, polyline.start_point_offset);
  }

  float first_reachable = 0;
  // first points already matched so start with 1 instead of 0
  // do not need to go to end of polyline so last point exluded
  for (unsigned int i = polyline.start_point_index + 1;
       i < polyline.end_point_index; i++) {
    auto data =
        _solver(line.start, line.end, polyline.polyline.get_point(i), epsilon);
    if (data.first == UNREACHABLE || data.last < first_reachable) {
      return UNREACHABLE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }

  // find on last line segment first reachable point
  auto data = _solver(polyline.polyline.get_point(polyline.end_point_index - 1),
                      polyline.polyline.get_point(polyline.end_point_index),
                      line.end, epsilon);

  return data.first;
}

float alt_godau_manhattan(PolylineRange polyline, LineSegment line,
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
    initial_dist += std::abs(coord);
  }
  if (initial_dist > epsilon) {
    return UNREACHABLE;
  }

  return _alt_godau_main<solve_manhattan>(polyline, line, epsilon);
}

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

// input epsilon squared
size_t alt_godau_euclidean_implicit(Polyline const &polyline, size_t j_,
                                    size_t j, size_t i_, size_t i,
                                    size_t restriction, float epsilon2) {
  if (1 != solve_implicit_euclidean_in(
               LineSegment(polyline.get_point(j_), polyline.get_point(j_ + 1)),
               polyline.get_point(restriction), polyline.get_point(i_),
               epsilon2)) {
    return (size_t)-1;
  } else if (j_ == j) {
    auto const res = solve_implicit_euclidean_in(
        LineSegment(polyline.get_point(j), polyline.get_point(j + 1)),
        polyline.get_point(restriction), polyline.get_point(i), epsilon2);
    size_t const results[3] = {(size_t)-1, restriction, i};
    return results[res];
  }

  size_t new_restriction = i_;
  for (unsigned int x = j_ + 1; x <= j; x++) {
    auto res = solve_implicit_euclidean_in(
        LineSegment(polyline.get_point(i_), polyline.get_point(i)),
        polyline.get_point(new_restriction), polyline.get_point(x), epsilon2);
    if (res == 0) {
      return (size_t)-1;
    }
    size_t const results[3] = {0, restriction, i};
    new_restriction = results[res];
  }

  if (is_line_reachable_euclidean(
          LineSegment(polyline.get_point(j), polyline.get_point(j + 1)),
          polyline.get_point(i), epsilon2)) {
    return i;
  }
  return (size_t)-1;
}

float alt_godau_chebyshev(PolylineRange polyline, LineSegment line,
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
    initial_dist = std::max(initial_dist, std::abs(coord));
  }
  if (initial_dist > epsilon) {
    return UNREACHABLE;
  }

  return _alt_godau_main<solve_chebyshev>(polyline, line, epsilon);
}

size_t alt_godau_minkowski_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

} // namespace DataStructures
