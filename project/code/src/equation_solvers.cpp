#include "config.h"
#include "datastructures.h"
#include "distance.h"
#include <cmath>
#include <stdexcept>

namespace DataStructures {

static inline void _solver_sanity_check(Point const &point1,
                                        Point const &point2,
                                        Point const &point3) {
  assert_compatible_points(point1, point2);
  assert_compatible_points(point1, point3);
  if (manhattan_distance(point1, point2) == 0) {
    throw std::runtime_error(
        "The two points given as a line segment are the same point.");
  }
}

ReachabilityData solve_manhattan(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif
  // TODO: Implement
  return {UNREACHABLE, UNREACHABLE};
}

ReachabilityData solve_maximum(Point const &point1, Point const &point2,
                               Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif
  // TODO: Implement
  return {UNREACHABLE, UNREACHABLE};
}

ReachabilityData solve_euclidean(Point const &point1, Point const &point2,
                                 Point const &point3, float epsilon) {
#if DEBUG
  _solver_sanity_check(point1, point2, point3);
#endif
  float dot_product = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    dot_product += (point2[i] - point1[i]) * (point1[i] - point3[i]);
  }

  float const a2 =
      DataStructures::unnormalized_euclidean_distance(point1, point2);
  float const a1 = 2 * dot_product;
  float const a0 =
      DataStructures::unnormalized_euclidean_distance(point1, point3) -
      epsilon * epsilon;

  float const discriminant = a1 * a1 - 4 * a2 * a0;
  if (discriminant < 0) {
    return {UNREACHABLE, UNREACHABLE};
  }

  float const root = std::sqrt(discriminant);
  float const t0 = (-a1 - root) / (2 * a2);
  float const t1 = (-a1 + root) / (2 * a2);
  if (t0 > 1 || t1 < 0) {
    return {UNREACHABLE, UNREACHABLE};
  }
  return {t0 < 0 ? 0 : t0, t1 > 1 ? 1 : t1};
}
} // namespace DataStructures
