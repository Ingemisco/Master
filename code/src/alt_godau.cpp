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
static inline float _alt_godau_main(SubPolyline &polyline, LineSegment &line,
                                    float epsilon) {
  if (polyline.start_point_index + 1 == polyline.end_point_index) {
    auto data = _solver(polyline.polyline.get_point(polyline.start_point_index),
                        polyline.polyline.get_point(polyline.end_point_index),
                        line.end, epsilon);
    if (data.first == EXPLICIT_UNREACHABLE || polyline.start_point_offset > data.second) {
      return EXPLICIT_UNREACHABLE;
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
    if (data.first == EXPLICIT_UNREACHABLE || data.second < first_reachable) {
      return EXPLICIT_UNREACHABLE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }

  // find on last line segment first reachable point
  auto data = _solver(polyline.polyline.get_point(polyline.end_point_index - 1),
                      polyline.polyline.get_point(polyline.end_point_index),
                      line.end, epsilon);

  return data.first;
}

float alt_godau_manhattan(SubPolyline polyline, LineSegment line,
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
    return EXPLICIT_UNREACHABLE;
  }

  return _alt_godau_main<solve_manhattan>(polyline, line, epsilon);
}

float alt_godau_euclidean(SubPolyline polyline, LineSegment line, float epsilon) {
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
    return EXPLICIT_UNREACHABLE;
  }

  return _alt_godau_main<solve_euclidean>(polyline, line, epsilon);
}

// computes the product <v - u | u - w>
static inline float scalar_product(Point const &u, Point const &v,
                                   Point const &w) {
  float dot_product = 0;
  for (unsigned int i = 0; i < u.dimension; i++) {
    dot_product += (v[i] - u[i]) * (u[i] - w[i]);
  }
  return dot_product;
}

static inline bool _alt_godau_euclidean_implicit_init(Point j_0, Point j_1,
                                                      Point r, Point i_,
                                                      float epsilon2) {
  float const r0dist = unnormalized_euclidean_distance(j_0, r);
  float const i_0dist = unnormalized_euclidean_distance(j_0, i_);
  if (r0dist <= epsilon2) {
    return i_0dist <= epsilon2;
  }

  float const ar0 = r0dist - epsilon2;
  float const ar1 = 2 * scalar_product(j_0, j_1, r);

  float const ai0 = i_0dist - epsilon2;
  float const ai1 = 2 * scalar_product(j_0, j_1, i_);

  float const a2 = unnormalized_euclidean_distance(j_0, j_1);

  float const dr = ar1 * ar1 - 4 * ar0 * a2;
  float const di = ai1 * ai1 - 4 * ai0 * a2;

  float const x = ai1 - ar1;
  float const y = dr + di - x * x;
  float const y2 = y * y;
  float const dprod = 4 * di * dr;

  float const i_1dist = unnormalized_euclidean_distance(j_1, i_);
  bool const first_condition =
      i_0dist <= epsilon2 ||
      (di >= 0 &&
       ((dr <= di && x >= 0) || (dr <= di && x < 0 && y >= 0 && y2 >= dprod) ||
        (dr > di && x >= 0 && (y <= 0 || y2 <= dprod))));

  bool const second_condition =
      i_1dist <= epsilon2 || (x <= 0 || y >= 0 || y2 <= dprod);

  return first_condition && second_condition;
}

static inline bool _is_in_01(float a1, float a2, float discriminant) {
  float const z = 2 * a2 + a1;
  return a1 <= 0 && discriminant <= a1 * a1 &&
         (z >= 0 || z * z <= discriminant);
}

// similar to the other implicit euclidean function. Returns 0 (unreachable), 1
// or
// 2. for the solutions [t01, t11], [t02, t12] on the line segment returns the
// index of the bigger of t01, t02 with same additional checks as above but also
// checks that t01 <= t12 and returns 0 if this is not the case
size_t solve_implicit_euclidean_in(LineSegment line, Point const &restriction,
                                   Point const &point, float epsilon2) {
#if DEBUG
  assert_compatible_points(line.start, restriction);
  assert_compatible_points(line.start, point);
#endif
  float const restriction_dist = unnormalized_euclidean_distance(restriction, line.start);

  float const point_dist = unnormalized_euclidean_distance(point, line.start);
  float const a0_1 = restriction_dist - epsilon2;
  float const a1_1 = 2 * scalar_product(line.start, line.end, restriction);

  float const a0_2 = point_dist - epsilon2;
  float const a1_2 = 2 * scalar_product(line.start, line.end, point);

  float const a2 = unnormalized_euclidean_distance(line.start, line.end);

  float const discriminant_1 = a1_1 * a1_1 - 4 * a0_1 * a2;
  float const discriminant_2 = a1_2 * a1_2 - 4 * a0_2 * a2;

  float const x = a1_1 - a1_2;
  float const y = discriminant_1 + discriminant_2 - x * x;
  float const y2 = y * y;
  float const discr_prod = 4 * discriminant_1 * discriminant_2;

  if (discriminant_2 < 0 ||
      (point_dist > epsilon2 && !_is_in_01(a1_2, a2, discriminant_2)) ||
      (x < 0 && y < 0 && discr_prod < y2)) {
    return 0;
  } else if (point_dist <= epsilon2 ||
             !((x >= 0 && discriminant_2 >= discriminant_1 &&
                (y <= 0 || y2 <= discr_prod)) ||
               (x < 0 && discriminant_2 < discriminant_1 &&
                (y >= 0 && y2 >= discr_prod)))) {
    return 1;
  }
  return 2;
}

// input epsilon squared
size_t alt_godau_euclidean_implicit(Polyline const &polyline, size_t j_,
                                    size_t j, size_t i_, size_t i,
                                    size_t restriction, float epsilon2) {

  if (!_alt_godau_euclidean_implicit_init(
          polyline.get_point(j_), polyline.get_point(j_ + 1),
          polyline.get_point(restriction), polyline.get_point(i_), epsilon2)) {
    return IMPLICIT_UNREACHABLE;
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
      return IMPLICIT_UNREACHABLE;
    }
    size_t const results[3] = {0, new_restriction, x};
    new_restriction = results[res];
  }

  if (is_line_reachable_euclidean(
          LineSegment(polyline.get_point(j), polyline.get_point(j + 1)),
          polyline.get_point(i), epsilon2)) {
    return i;
  }
  return IMPLICIT_UNREACHABLE;
}

float alt_godau_chebyshev(SubPolyline polyline, LineSegment line, float epsilon) {
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
    return EXPLICIT_UNREACHABLE;
  }

  return _alt_godau_main<solve_chebyshev>(polyline, line, epsilon);
}

size_t alt_godau_minkowski_implicit(Polyline const &polyline,
                                    Point const &line_start,
                                    Point const &line_end, float epsilon);

} // namespace DataStructures
