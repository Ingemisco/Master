#include "distance.h"
#include "config.h"
#include "datastructures.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace DataStructures {

#define NO_POINT_REACHABLE -1

typedef struct {
  float upper_first;
  float right_first;
} ReachabilityTable;

// equivalent to unnormalized_minkowski_distance(point1, point2, 2)
float unnormalized_euclidean_distance(Point const &point1,
                                      Point const &point2) {
#if DEBUG
  assert_compatible_points(point1, point2);
#endif

  float result = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float const diff = point1[i] - point2[i];
    result += diff * diff;
  }
  return result;
}

// equivalent to minkowski_distance(point1, point2, 2)
float euclidean_distance(Point const &point1, Point const &point2) {
  return std::sqrt(unnormalized_euclidean_distance(point1, point2));
}

float manhattan_distance(Point const &point1, Point const &point2) {
#if DEBUG
  assert_compatible_points(point1, point2);
#endif

  float result = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float const diff = point1[i] - point2[i];
    result += std::abs(diff);
  }
  return result;
}

float maximum_norm_distance(Point const &point1, Point const &point2) {
#if DEBUG
  assert_compatible_points(point1, point2);
#endif

  float result = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float const diff = point1[i] - point2[i];
    result = std::max(std::abs(diff), result);
  }
  return result;
}

// assumes n >= 0
template <typename T> static inline T _integer_exponentiation(T a, int n) {
  T result = 1;
  while (n) {
    if (n & 1) {
      result *= a;
    }
    a *= a;
    n >>= 1;
  }
  return result;
}

// assumes n >= 0
float integer_exponentiation(float a, int n) {
  return _integer_exponentiation(a, n);
}

// assumes n >= 0
int integer_exponentiation(int a, int n) {
  return _integer_exponentiation(a, n);
}

float unnormalized_minkowski_distance(Point const &point1, Point const &point2,
                                      int p) {
#if DEBUG
  assert_compatible_points(point1, point2);
  if (p > 0) {
    throw std::runtime_error("Order of the norm must be greater than zero!");
  }
#endif

  float result = 0;
  for (unsigned int i = 0; i < point1.dimension; i++) {
    float const diff = point1[i] - point2[i];
    result += std::abs(integer_exponentiation(diff, p));
  }
  return result;
}

float minkowski_distance(Point const &point1, Point const &point2, int p) {
  return std::pow(unnormalized_minkowski_distance(point1, point2, p), 1.0 / p);
}

float first_reachable_euclidean(Polyline const &polyline1,
                                Polyline const &polyline2, float epsilon,
                                size_t line_start, size_t line_end,
                                size_t point) {
  auto const v1 = polyline1.get_point(line_start);
  auto const v2 = polyline1.get_point(line_end);
  auto const u = polyline2.get_point(point);

  float epsilon_squared = epsilon * epsilon;
  float const dist_u_v1 = unnormalized_euclidean_distance(v1, u);
  if (dist_u_v1 <= epsilon_squared) {
    return 0;
  }

  // get coefficients of the polynomial of degree 2 whose roots are the
  // intersections with the line segment if they exist
  float const a0 = unnormalized_euclidean_distance(v2, v1) - epsilon_squared;
  float a1 = 0;
  float const a2 = dist_u_v1;
  for (unsigned int i = 0; i < u.dimension; i++) {
    a1 += (u[i] - v1[i]) * (v2[i] - v1[i]);
  }
  a1 *= 2;

  float const val1 = 1 / a2;
  float const val2 = a1 * val1 / 2;
  float const val3 = a0 * val1;
  float const discriminant = val2 * val2 - val3;
  if (discriminant < 0) {
    return NO_POINT_REACHABLE;
  }

  return -val2 - std::sqrt(discriminant);
}

bool frechet_distance_at_most_euclidean_exact(
    DataStructures::Polyline const &polyline1,
    DataStructures::Polyline const &polyline2, float epsilon) {
#if DEBUG
  assert_compatible_polylines(polyline1, polyline2);
#endif

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

// TODO: Use algorithm from Alt & Godau for the decision variant
bool frechet_distance_decision(DataStructures::Polyline const &polyline1,
                               DataStructures::Polyline const &polyline2,
                               [[maybe_unused]] float epsilon,
                               [[maybe_unused]] Distance distance) {
#if DEBUG
  assert_compatible_polylines(polyline1, polyline2);
#endif

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

// uses Alt and Godau algorithm with square roots to explicitly compute first
// reachable points and track them
bool frechet_distance_decision_eucliean_sqrt(
    DataStructures::Polyline const &polyline1,
    DataStructures::Polyline const &polyline2, float epsilon) {
#if DEBUG
  assert_compatible_polylines(polyline1, polyline2);
#endif

  Polyline const *smaller_line = &polyline2;
  Polyline const *larger_line = &polyline1;
  if (polyline1.point_count < polyline2.point_count) {
    smaller_line = &polyline1;
    larger_line = &polyline1;
  }
  size_t const table_size = smaller_line->point_count;

  ReachabilityTable *reachable_table = new ReachabilityTable[table_size];
  for (unsigned int i = 0; i < polyline1.point_count; i++) {
    for (unsigned int j = 0; j < polyline2.point_count; j++) {
      // TODO: complete
    }
  }

  delete[] reachable_table;
  return false;
}

// TODO: Use algorithm from Alt & Godau
float frechet_distance(Polyline const &polyline1, Polyline const &polyline2,
                       [[maybe_unused]] Distance distance) {
#if DEBUG
  assert_compatible_polylines(polyline1, polyline2);
#endif
  return 0;
}
} // namespace DataStructures
