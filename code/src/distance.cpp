#include "distance.h"
#include "config.h"
#include "datastructures.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace DataStructures {

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

float chebyshev_distance(Point const &point1, Point const &point2) {
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

} // namespace DataStructures
