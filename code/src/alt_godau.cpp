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
#include "simplification.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <utility>

namespace DataStructures {

template <typename F, typename L, std::pair<F, L> _solver(Polyline const &, size_t, size_t, size_t, float), std::pair<F const, L const> const UNREACHABLE_INTERVAL, F const UNREACHABLE_VALUE, F zero>
static inline F _alt_godau_main(Polyline const &polyline, size_t j_, size_t j, F t, size_t i_, size_t i, float epsilon) {
  if (j_ == j) {
		std::pair<F, L> data = _solver(polyline, j, j + 1, i, epsilon);
    if (data == UNREACHABLE_INTERVAL || !(t <= data.second)) {
      return UNREACHABLE_VALUE;
    }
    return std::max(data.first, t);
  }

  F first_reachable = zero;
  // first points already matched so start with 1 instead of 0
  // do not need to go to end of polyline so last point exluded
  for (unsigned int index = j_ + 1; index < j + 1; index++) {
    auto data = _solver(polyline, i_, i, index, epsilon);
    if (data == UNREACHABLE_INTERVAL || !(first_reachable <= data.second)) {
      return UNREACHABLE_VALUE;
    }

    first_reachable = std::max(first_reachable, data.first);
  }

  // find on last line segment first reachable point
  auto data = _solver(polyline, j, j + 1, i, epsilon);

  return data.first;
}

float alt_godau_manhattan(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
  return _alt_godau_main<float, float, solve_manhattan, EMPTY_INTERVAL_EXPLICIT, EXPLICIT_UNREACHABLE, 0.0f>(polyline, i_, i, t, line_start, line_end, epsilon);
}

float alt_godau_euclidean(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
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

static inline size_t _ag_implicit_single_line(Polyline const &polyline, size_t line_start, size_t line_end, size_t restriction, size_t u, float epsilon2) {
  float const r_dist = unnormalized_euclidean_distance(polyline, restriction, line_start);
  float const u_dist = unnormalized_euclidean_distance(polyline, u, line_start);
  float const u_end_dist = unnormalized_euclidean_distance(polyline, u, line_end);

	// parameters to the equations alpha t^2 + 2 beta t + gammma = 0, alpha for both the same
  float const gamma_r = r_dist - epsilon2;
  float const beta_r = scalar_product(polyline, line_start, line_end, restriction);

  float const gamma_u = u_dist - epsilon2;
  float const beta_u = scalar_product(polyline, line_start, line_end, u);

  float const alpha = unnormalized_euclidean_distance(polyline, line_start, line_end);

  float const D_r = beta_r * beta_r - gamma_r * alpha;
  float const D_u = beta_u * beta_u - gamma_u * alpha;

  float const x = beta_r - beta_u;
  float const y = D_r + D_u - x * x;
  float const y2 = y * y;
  float const discr_prod = 4 * D_r * D_u;

	// u does not reach the line segment 
	if (!(u_dist <= epsilon2 || u_end_dist <= epsilon2 || 
		(0 >= beta_u && beta_u >= -alpha && u_dist * alpha - beta_u * beta_u <= epsilon2 * alpha))) {
		return IMPLICIT_UNREACHABLE;
	}

	// is first sol of r before first sol of u 
	if ((x >= 0 && D_u <= D_r) || 
		(x >= 0 && D_u >= D_r && (y <= 0 || y2 <= discr_prod)) || 
		(x <= 0 && D_u <= D_r && (y >= 0 && y2 >= discr_prod))) {
		return u;
	}

	// first sol of r is after that of u, so check if second sol of u is after first of r 
	if (x >= 0 || y >= 0 || y2 <= discr_prod) {
		return restriction;
	}
	
	return IMPLICIT_UNREACHABLE;
}

static inline size_t _ag_implicit_multi_line(Polyline const &polyline, size_t i_, size_t i, size_t j_, size_t j, float epsilon2) {
	size_t restriction = i_;
	for (size_t k = j_ + 1; k <= j; k++) {
		restriction = _ag_implicit_single_line(polyline, i_, i, restriction, k, epsilon2);
		if (restriction == IMPLICIT_UNREACHABLE) {
			return IMPLICIT_UNREACHABLE;
		}
	}

	return _ag_implicit_single_line(polyline, j, j + 1, j, i, epsilon2);
}

// input epsilon squared
size_t alt_godau_euclidean_implicit(Polyline const &polyline, size_t j_, size_t j, size_t i_, size_t i, size_t restriction, float epsilon2) {
	if (j_ == j) {
		return _ag_implicit_single_line(polyline, j, j + 1, restriction, i, epsilon2);
  }

	return _ag_implicit_multi_line(polyline, i_, i, j_, j, epsilon2);
}

float alt_godau_chebyshev(Polyline const &polyline, size_t i_, size_t i, float t, size_t line_start, size_t line_end, float epsilon) {
  return _alt_godau_main<float, float, solve_chebyshev, EMPTY_INTERVAL_EXPLICIT, EXPLICIT_UNREACHABLE, 0.0f>(polyline, i_, i, t, line_start, line_end, epsilon);
}

size_t alt_godau_minkowski_implicit(Polyline const &polyline, size_t line_start, size_t line_end, float epsilon);

// epsilon must be squared
FRValue alt_godau_euclidean_semiexplicit(Polyline const &polyline, size_t i_, size_t i, FRValue t, size_t line_start, size_t line_end, float epsilon2) {
  return _alt_godau_main<FRValue, LRValue, solve_euclidean_se, EMPTY_INTERVAL_SEMIEXPLICIT, SEMIEXPLICIT_UNREACHABLE, FRValue(0,0)>(polyline, i_, i, t, line_start, line_end, epsilon2);
}

} // namespace DataStructures
