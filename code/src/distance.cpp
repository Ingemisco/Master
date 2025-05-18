#include "distance.h"
#include "config.h"
#include "datastructures.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

using DataStructures::Polyline;

namespace DataStructures {

// equivalent to unnormalized_minkowski_distance(point1, point2, 2)
float unnormalized_euclidean_distance(Polyline const &polyline, size_t point1, size_t point2) {
  float result = 0;
  for (unsigned int i = 0; i < polyline.dimension; i++) {
    float const diff = polyline[point1, i] - polyline[point2,i];
    result += diff * diff;
  }
  return result;
}

// equivalent to minkowski_distance(point1, point2, 2)
float euclidean_distance(Polyline const &polyline, size_t point1, size_t point2) {
  return std::sqrt(unnormalized_euclidean_distance(polyline, point1, point2));
}

float manhattan_distance(Polyline const &polyline, size_t point1, size_t point2) {
  float result = 0;
  for (unsigned int i = 0; i < polyline.dimension; i++) {
    float const diff = polyline[point1, i] - polyline[point2, i];
    result += std::abs(diff);
  }
  return result;
}

float chebyshev_distance(Polyline const &polyline, size_t point1, size_t point2) {
  float result = 0;
  for (unsigned int i = 0; i < polyline.dimension; i++) {
    float const diff = polyline[point1, i] - polyline[point2, i];
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

float unnormalized_minkowski_distance(Polyline const &polyline, size_t point1, size_t point2, int p) {
  float result = 0;
  for (unsigned int i = 0; i < polyline.dimension; i++) {
    float const diff = polyline[point1, i] - polyline[point2, i];
    result += std::abs(integer_exponentiation(diff, p));
  }
  return result;
}

float minkowski_distance(Polyline const &polyline, size_t point1, size_t point2, int p) {
  return std::pow(unnormalized_minkowski_distance(polyline, point1, point2, p), 1.0 / p);
}

} // namespace DataStructures
