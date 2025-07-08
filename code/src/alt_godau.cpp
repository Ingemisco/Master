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
#include <utility>

namespace DataStructures {

template <typename F, typename L, std::pair<F, L> _solver(Polyline const &, size_t, size_t, size_t, float), std::pair<F const, L const> const UNREACHABLE_INTERVAL, F const UNREACHABLE_VALUE, F zero>
static inline F _alt_godau_main(Polyline const &polyline, size_t i_, size_t i, F t, size_t line_start, size_t line_end, float epsilon) {
  if (i_ + 1 == i) {
		std::pair<F, L> data = _solver(polyline, i_, i, line_end, epsilon);
    if (data == UNREACHABLE_INTERVAL || !(t <= data.second)) {
      return UNREACHABLE_VALUE;
    }
    return std::max(data.first, t);
  }

  F first_reachable = zero;
  // first points already matched so start with 1 instead of 0
  // do not need to go to end of polyline so last point exluded
  for (unsigned int j = i_ + 1; j < i; j++) {
    auto data = _solver(polyline, line_start, line_end, j, epsilon);
    if (data == UNREACHABLE_INTERVAL || !(first_reachable <= data.second)) {
      return UNREACHABLE_VALUE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }

  // find on last line segment first reachable point
  auto data = _solver(polyline, i - 1, i, line_end, epsilon);

  return data.first;
}

float alt_godau_manhattan(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
  // compute distance of first point on polyline (p[start] + t(p[start + 1] -
  // p[start])) to first point of line segment and compare if distance greater
  // then epsilon, if so, the line segment is too far away
  float initial_dist = 0;
  for (unsigned int j = 0; j < polyline.dimension; j++) {
    float const coord = polyline[i_, j] + t * (polyline[i_ + 1, j] - polyline[i_, j]) - polyline[line_start, j];
    initial_dist += std::abs(coord);
  }
  if (initial_dist > epsilon) {
    return EXPLICIT_UNREACHABLE;
  }

  return _alt_godau_main<float, float, solve_manhattan, EMPTY_INTERVAL_EXPLICIT, EXPLICIT_UNREACHABLE, 0.0f>(polyline, i_, i, t, line_start, line_end, epsilon);
}

float alt_godau_euclidean(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
  // compute distance of first point on polyline (p[start] + t(p[start + 1] -
  // p[start])) to first point of line segment and compare if distance greater
  // then epsilon, if so, the line segment is too far away
	// ALL of this is not required???
  // float initial_dist = 0;
  // for (unsigned int j = 0; j < polyline.dimension; j++) {
  //   float const coord =
  //       polyline[i_, j] + t * (polyline[i_ + 1, j] - polyline[i_, j]) - polyline[line_start, j];
  //   initial_dist += coord * coord;
  // }
  // if (initial_dist > epsilon * epsilon) {
  //   return EXPLICIT_UNREACHABLE;
  // }

  return _alt_godau_main<float, float, solve_euclidean, EMPTY_INTERVAL_EXPLICIT, EXPLICIT_UNREACHABLE, 0.0f>(polyline, i_, i, t, line_start, line_end, epsilon);
}

// computes the product <v - u | u - w>
static inline float scalar_product(Polyline const &polyline, size_t u, size_t v, size_t w) {
  float dot_product = 0;
  for (unsigned int i = 0; i < polyline.dimension; i++) {
    dot_product += (polyline[v, i] - polyline[u, i]) * (polyline[u, i] - polyline[w, i]);
  }
  return dot_product;
}

static inline bool _alt_godau_euclidean_implicit_init(Polyline const &polyline, size_t j_0, size_t j_1, size_t r, size_t i_, float epsilon2) {
  float const r0dist = unnormalized_euclidean_distance(polyline, j_0, r);
  float const i_0dist = unnormalized_euclidean_distance(polyline, j_0, i_);
  if (r0dist <= epsilon2) {
    return i_0dist <= epsilon2;
  }

  float const ar0 = r0dist - epsilon2;
  float const ar1 = scalar_product(polyline, j_0, j_1, r);

  float const ai0 = i_0dist - epsilon2;
  float const ai1 = scalar_product(polyline, j_0, j_1, i_);

  float const a2 = unnormalized_euclidean_distance(polyline, j_0, j_1);

  float const dr = ar1 * ar1 - ar0 * a2;
  float const di = ai1 * ai1 - ai0 * a2;

  float const x = ai1 - ar1;
  float const y = dr + di - x * x;
  float const y2 = y * y;
  float const dprod = 4 * di * dr;

  float const i_1dist = unnormalized_euclidean_distance(polyline, j_1, i_);
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
size_t solve_implicit_euclidean_in(Polyline const &polyline, size_t line_start, size_t line_end, size_t restriction, size_t point, float epsilon2) {
  float const restriction_dist = unnormalized_euclidean_distance(polyline, restriction, line_start);

  float const point_dist = unnormalized_euclidean_distance(polyline, point, line_start);
  float const a0_1 = restriction_dist - epsilon2;
  float const a1_1 = 2 * scalar_product(polyline, line_start, line_end, restriction);

  float const a0_2 = point_dist - epsilon2;
  float const a1_2 = 2 * scalar_product(polyline, line_start, line_end, point);

  float const a2 = unnormalized_euclidean_distance(polyline, line_start, line_end);

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
size_t alt_godau_euclidean_implicit(Polyline const &polyline, size_t j_, size_t j, size_t i_, size_t i, size_t restriction, float epsilon2) {
  if (!_alt_godau_euclidean_implicit_init(polyline, j_, j_ + 1, restriction, i_, epsilon2)) {
    return IMPLICIT_UNREACHABLE;
  } else if (j_ == j) {
    auto const res = solve_implicit_euclidean_in(polyline, j, j + 1, restriction, i, epsilon2);
    size_t const results[3] = {(size_t)-1, restriction, i};
    return results[res];
  }

  size_t new_restriction = i_;
  for (unsigned int x = j_ + 1; x <= j; x++) {
    auto res = solve_implicit_euclidean_in(polyline, i_, i, new_restriction, x, epsilon2);
    if (res == 0) {
      return IMPLICIT_UNREACHABLE;
    }
    size_t const results[3] = {0, new_restriction, x};
    new_restriction = results[res];
  }

  if (is_line_reachable_euclidean(polyline, j, j + 1, i, epsilon2)) {
    return i;
  }
  return IMPLICIT_UNREACHABLE;
}

float alt_godau_chebyshev(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
  // compute distance of first point on polyline (p[start] + t(p[start + 1] -
  // p[start])) to first point of line segment and compare if distance greater
  // then epsilon, if so, the line segment is too far away
  float initial_dist = 0;
  for (unsigned int j = 0; j < polyline.dimension; j++) {
    float const coord =
        polyline[i_, j] + t * (polyline[i_ + 1, j] - polyline[i_, j]) - polyline[line_start, j];
    initial_dist = std::max(initial_dist, std::abs(coord));
  }
  if (initial_dist > epsilon) {
    return EXPLICIT_UNREACHABLE;
  }

  return _alt_godau_main<float, float, solve_chebyshev, EMPTY_INTERVAL_EXPLICIT, EXPLICIT_UNREACHABLE, 0.0f>(polyline, i_, i, t, line_start, line_end, epsilon);
}

size_t alt_godau_minkowski_implicit(Polyline const &polyline, size_t line_start, size_t line_end, float epsilon);

FRValue alt_godau_euclidean_semiexplicit(Polyline const &polyline, size_t i_, size_t i, FRValue t, size_t line_start, size_t line_end, float epsilon) {
  return _alt_godau_main<FRValue, LRValue, solve_euclidean_se, EMPTY_INTERVAL_SEMIEXPLICIT, SEMIEXPLICIT_UNREACHABLE, FRValue(0,0)>(polyline, i_, i, t, line_start, line_end, epsilon);
}

} // namespace DataStructures
